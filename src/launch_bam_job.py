#!/usr/bin/env python
"""
Launch the stand-alone BAM R script as a SageMaker Processing job.

Why a Processing job: the script runs in a separate AWS-managed container, NOT
in this notebook/Studio session. Once submitted it keeps running even if you
close the window or your admin-mandated session times out. Results land in S3.

This uses boto3 directly (the sagemaker Python SDK is not fully installed on
this instance). The container image is the rocker/tidyverse image we mirrored
into ECR; the script install.packages() any extras (emmeans, gratia, ...) at
runtime, which works because the job has internet egress.

Usage (from this instance, inside SageMaker):
    python launch_bam_job.py                          # confirmation screen (default)
    python launch_bam_job.py --screen initial         # the larger initial screen
    python launch_bam_job.py --wait                   # also stream status here
    python launch_bam_job.py --instance-type ml.r7i.4xlarge

Each screen has its own S3 prefix (code + input + output). Upload code/inputs
with sync_to_s3.py (or the aws s3 cp commands it prints) before launching.
"""
import argparse
import os
import sys
import time

import boto3

# ── Account / location constants ─────────────────────────────────────────────
# Account-specific values are read from the environment so nothing account-
# specific is committed. Set these before launching a job:
#   export FINOMENA_AWS_ACCOUNT=<your 12-digit account id>
#   export FINOMENA_AWS_REGION=us-east-1                    # optional, defaults below
#   export FINOMENA_S3_BUCKET=sagemaker-<region>-<account>  # optional, derived below
#   export FINOMENA_ROLE_ARN=arn:aws:iam::<account>:role/<role>  # optional, derived
ACCOUNT = os.environ.get("FINOMENA_AWS_ACCOUNT")
if not ACCOUNT:
    sys.exit("Set FINOMENA_AWS_ACCOUNT to your 12-digit AWS account ID "
             "(see the account/location constants block in launch_bam_job.py).")
REGION = os.environ.get("FINOMENA_AWS_REGION", "us-east-1")
BUCKET = os.environ.get("FINOMENA_S3_BUCKET", f"sagemaker-{REGION}-{ACCOUNT}")
PREFIX = "finomena/bam-run"
IMAGE_URI = os.environ.get(
    "FINOMENA_IMAGE_URI",
    f"{ACCOUNT}.dkr.ecr.{REGION}.amazonaws.com/finomena-r:4.4.2")
ROLE_ARN = os.environ.get(
    "FINOMENA_ROLE_ARN",
    f"arn:aws:iam::{ACCOUNT}:role/sagemaker_execution_role")

# Each screen is a self-contained S3 prefix under PREFIX. The local Data/ dir
# uses "Initial screen" / "Confirmation"; map the CLI flag to both the S3
# sub-prefix and the local input dir (see sync_to_s3.py).
# Per-screen settings. error_rate maps to the multiple-testing correction the
# R script applies (FDR -> Benjamini-Hochberg, FWER -> Holm): the initial screen
# is a discovery screen (FDR), the confirmation screen is confirmatory (FWER).
SCREENS = {
    "confirmation": {"s3": "confirmation", "local": "Data/Confirmation",
                     "error_rate": "FWER",
                     # ablation param parametrization (balanced 2x2 here)
                     "abl_param_mode": "interaction"},
    "initial": {"s3": "initial", "local": "Data/Initial screen",
                "error_rate": "FDR",
                # Genotype*Drug non-estimable on initial -> single-factor param
                "abl_param_mode": "condition"},
}

# Fixed container mount points used by SageMaker Processing.
C_CODE = "/opt/ml/processing/code"
C_INPUT = "/opt/ml/processing/input"
C_OUTPUT = "/opt/ml/processing/output"


C_RESUME = "/opt/ml/processing/resume"


def build_job(name, screen, instance_type, volume_gb, max_runtime_h,
              ablation=False, plate_filter="", abl_only="all", abl_cv="1",
              abl_baseline="", row_filter="", cv_fold_by="", post_only=False,
              skip_bam=False):
    base = f"s3://{BUCKET}/{PREFIX}/{SCREENS[screen]['s3']}"
    s3_code = f"{base}/code"
    s3_input = f"{base}/input"
    # Ablation writes to its own output sub-prefix so it doesn't clobber the
    # production pipeline's output/.
    s3_output = f"{base}/ablation_output" if ablation else f"{base}/output"

    # Common env shared by both modes.
    env = {
        "OPENBLAS_NUM_THREADS": "118",   # backstop under the 128 cap
        "CRAN_MIRROR":
            "https://packagemanager.posit.co/cran/__linux__/noble/latest",
    }

    if ablation:
        # The ablation reads the already-aggregated frame from the production
        # output prefix (output/aggregated_data.csv), mounted at C_INPUT.
        entrypoint = ["Rscript", f"{C_CODE}/ablation_study.R"]
        env.update({
            "ABL_AGG_CSV":   f"{C_INPUT}/aggregated_data.csv",
            "ABL_OUT_DIR":   C_OUTPUT,
            "ABL_PARAM_MODE": SCREENS[screen]["abl_param_mode"],
            "ABL_PLATE_FILTER": plate_filter,
            "ABL_ROW_FILTER": row_filter,
            "ABL_RHO": "0.2",
            "ABL_CV": abl_cv,
            "ABL_ONLY": abl_only,
        })
        # Only set these when overriding — empty would clobber the script's
        # unset-defaults.
        if abl_baseline:
            env["ABL_BASELINE"] = abl_baseline
        if cv_fold_by:
            env["ABL_CV_FOLD_BY"] = cv_fold_by
        # Mount the production aggregated_data.csv as the ablation input.
        s3_input = f"{base}/output"
    elif screen == "initial":
        # Initial screen runs the FULL pipeline in one execution via the wrapper:
        # BAM fit -> consolidate -> robust Glass/SSMD -> movement/Check-A addons -> plots.
        # The wrapper threads C_OUTPUT through every step's env so all artifacts
        # land in the single SageMaker output channel.
        entrypoint = ["bash", f"{C_CODE}/run_pipeline.sh"]
        env.update({
            "FINOMENA_INPUT_DIR": C_INPUT,
            "FINOMENA_OUTPUT_DIR": C_OUTPUT,
            "FINOMENA_HELPERS_DIR": f"{C_CODE}/helpers/scripts",
            "FINOMENA_ERROR_RATE": SCREENS[screen]["error_rate"],
            "FINOMENA_ANALYSIS_MODE": "initial",
            # Post-processing reads the BAM output channel directly.
            "FINOMENA_ADDON_PROD_DIR": C_OUTPUT,
            "FINOMENA_ADDON_AGG_CSV": f"{C_OUTPUT}/aggregated_data.csv",
            "FINOMENA_ADDON_OUTPUT_DIR": C_OUTPUT,
            "FINOMENA_PLOTS_DIR": f"{C_OUTPUT}/Plots",
            "FINOMENA_PLOT_AGG_CSV": f"{C_OUTPUT}/aggregated_data.csv",
        })
        if post_only:
            # Resume: skip the BAM + consolidate + SSMD (STEPs 1-3) and run only
            # the addons + plots, seeding from a prior run's output prefix mounted
            # at C_RESUME.
            env["FINOMENA_POST_ONLY"] = "1"
            env["FINOMENA_RESUME_DIR"] = C_RESUME
        elif skip_bam:
            # Resume: skip ONLY the multi-hour BAM fit (STEP 1); re-run consolidate
            # + SSMD + addons + plots (STEPs 2-5) against the seeded per-plate BAM
            # output. This is the mode that reprocesses EVERY plate the BAM fit
            # (the plate list is derived from the seeded folders by run_pipeline.sh).
            env["FINOMENA_SKIP_BAM"] = "1"
            env["FINOMENA_RESUME_DIR"] = C_RESUME
    else:
        entrypoint = ["Rscript", f"{C_CODE}/stand-alone_script.r"]
        env.update({
            "FINOMENA_INPUT_DIR": C_INPUT,
            "FINOMENA_OUTPUT_DIR": C_OUTPUT,
            "FINOMENA_HELPERS_DIR": f"{C_CODE}/helpers/scripts",
            "FINOMENA_ERROR_RATE": SCREENS[screen]["error_rate"],
        })

    processing_inputs = [
        {
            "InputName": "code",
            "S3Input": {
                "S3Uri": s3_code,
                "LocalPath": C_CODE,
                "S3DataType": "S3Prefix",
                "S3InputMode": "File",
                "S3DataDistributionType": "FullyReplicated",
            },
        },
        {
            "InputName": "input",
            "S3Input": {
                "S3Uri": s3_input,
                "LocalPath": C_INPUT,
                "S3DataType": "S3Prefix",
                "S3InputMode": "File",
                "S3DataDistributionType": "FullyReplicated",
            },
        },
    ]
    if post_only or skip_bam:
        # Mount the prior run's output prefix read-only so the wrapper can seed
        # the BAM artifacts (per-plate folders, aggregated_data, bam_residual_df,
        # master_results, and — for POST_ONLY — the consolidated ALL_PLATES/SSMD
        # tables). POST_ONLY reuses those tables (runs STEPs 4-5); SKIP_BAM rebuilds
        # them from the seeded per-plate folders (runs STEPs 2-5).
        processing_inputs.append({
            "InputName": "resume",
            "S3Input": {
                "S3Uri": s3_output,
                "LocalPath": C_RESUME,
                "S3DataType": "S3Prefix",
                "S3InputMode": "File",
                "S3DataDistributionType": "FullyReplicated",
            },
        })

    return {
        "ProcessingJobName": name,
        "RoleArn": ROLE_ARN,
        "AppSpecification": {
            "ImageUri": IMAGE_URI,
            "ContainerEntrypoint": entrypoint,
        },
        "Environment": env,
        "ProcessingInputs": processing_inputs,
        "ProcessingOutputConfig": {
            "Outputs": [
                {
                    "OutputName": "output",
                    "S3Output": {
                        "S3Uri": s3_output,
                        "LocalPath": C_OUTPUT,
                        "S3UploadMode": "EndOfJob",
                    },
                }
            ]
        },
        "ProcessingResources": {
            "ClusterConfig": {
                "InstanceCount": 1,
                "InstanceType": instance_type,
                "VolumeSizeInGB": volume_gb,
            }
        },
        "StoppingCondition": {
            "MaxRuntimeInSeconds": int(max_runtime_h * 3600)
        },
        "Tags": [
            {"Key": "project", "Value": "finomena-bam"},
            # Optional cost-attribution owner tag; set FINOMENA_OWNER_TAG to add it.
            *([{"Key": "owner", "Value": os.environ["FINOMENA_OWNER_TAG"]}]
              if os.environ.get("FINOMENA_OWNER_TAG") else []),
        ],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen", choices=sorted(SCREENS), default="confirmation",
                   help="Which dataset to run (default: confirmation, the "
                        "smaller test set).")
    # ml.r7i.48xlarge = 192 vCPU / 1536 GB RAM. This is a large, memory-heavy
    # bam() job, so default to the biggest r7i. The r7i family has nonzero
    # processing quota on this account (ml.r5 did not); the 48xlarge size still
    # needs its own quota — request an increase in Service Quotas if it's 0.
    p.add_argument("--instance-type", default="ml.r7i.48xlarge")
    p.add_argument("--volume-gb", type=int, default=100)
    p.add_argument("--max-runtime-h", type=float, default=12.0)
    p.add_argument("--ablation", action="store_true",
                   help="Run the formula ablation study (ablation_study.R) on "
                        "this screen's aggregated_data.csv instead of the main "
                        "pipeline. Writes to <screen>/ablation_output/.")
    p.add_argument("--plate-filter", default="",
                   help="Ablation only: subset to plates whose name contains "
                        "this string, e.g. 'Plate 2' (Phase B).")
    p.add_argument("--abl-only", default="all",
                   help="Ablation only: comma-separated variant keys to run "
                        "(baseline always included), or 'all'. For quick tests.")
    p.add_argument("--abl-cv", default="1",
                   help="Ablation only: '1' to run the 5-fold CV gate, '0' to skip.")
    p.add_argument("--abl-baseline", default="",
                   help="Ablation only: variant key to use as the similarity "
                        "reference (default: full production baseline). Set to a "
                        "tractable cell, e.g. cond30__animalre, when the full "
                        "baseline can't be fit at scale.")
    p.add_argument("--row-filter", default="",
                   help="Ablation only: keep only well-rows in this set, e.g. "
                        "'ABCD' for rows A-D (drops E-H). Tractability lever.")
    p.add_argument("--cv-fold-by", default="",
                   help="Ablation only: 'replicate' (fold by replicate; needed "
                        "for initial screen) or 'animal' (default).")
    p.add_argument("--post-only", action="store_true",
                   help="Initial screen only: SKIP the BAM fit + consolidate + "
                        "SSMD (STEPs 1-3) and run only the addons + plots (STEPs "
                        "4-5), seeding from the prior run's <screen>/output/ "
                        "prefix. Use after a STEP 4/5 code fix to avoid the "
                        "multi-hour BAM refit. Defaults to a smaller instance.")
    p.add_argument("--skip-bam", action="store_true",
                   help="Initial screen only: SKIP only the BAM fit (STEP 1) and "
                        "re-run consolidate + SSMD + addons + plots (STEPs 2-5) "
                        "against the prior run's per-plate BAM output, seeding from "
                        "<screen>/output/. Use to rebuild the consolidated tables "
                        "and figures for EVERY plate the BAM fit (plate list is "
                        "auto-derived) without the multi-hour refit. Defaults to a "
                        "smaller instance.")
    p.add_argument("--wait", action="store_true",
                   help="Stream status here until the job finishes. Safe to "
                        "Ctrl-C — the job keeps running in AWS regardless.")
    args = p.parse_args()

    # Both resume modes are initial-screen-only and mutually exclusive with each
    # other and with ablation.
    for flag, on in (("--post-only", args.post_only), ("--skip-bam", args.skip_bam)):
        if on and args.screen != "initial":
            p.error(f"{flag} is only valid with --screen initial")
        if on and args.ablation:
            p.error(f"{flag} and --ablation are mutually exclusive")
    if args.post_only and args.skip_bam:
        p.error("--post-only and --skip-bam are mutually exclusive")
    # Neither resume mode does a BAM fit, so it doesn't need the giant memory
    # instance; drop to a small one unless --instance-type was overridden. SKIP_BAM
    # re-runs consolidate/SSMD over all plates (a bit heavier than post-only's
    # plots-only) but is still far below the fit's footprint.
    if (args.post_only or args.skip_bam) and args.instance_type == "ml.r7i.48xlarge":
        args.instance_type = "ml.r7i.2xlarge"

    sm = boto3.client("sagemaker", region_name=REGION)
    # Job names must be unique; stamp with epoch seconds.
    kind = ("post" if args.post_only else "skipbam" if args.skip_bam
            else "abl" if args.ablation else "bam")
    name = f"finomena-{kind}-{args.screen}-{int(time.time())}"

    sm.create_processing_job(**build_job(
        name, args.screen, args.instance_type, args.volume_gb,
        args.max_runtime_h, ablation=args.ablation,
        plate_filter=args.plate_filter, abl_only=args.abl_only,
        abl_cv=args.abl_cv, abl_baseline=args.abl_baseline,
        row_filter=args.row_filter, cv_fold_by=args.cv_fold_by,
        post_only=args.post_only, skip_bam=args.skip_bam))

    sub = "ablation_output" if args.ablation else "output"
    out = f"s3://{BUCKET}/{PREFIX}/{SCREENS[args.screen]['s3']}/{sub}/"
    console = (f"https://{REGION}.console.aws.amazon.com/sagemaker/home"
               f"?region={REGION}#/processing-jobs/{name}")
    print(f"Submitted processing job: {name}")
    mode = ("ablation" if args.ablation else
            "resume: post-only (STEPs 4-5)" if args.post_only else
            "resume: skip-BAM (STEPs 2-5, all plates)" if args.skip_bam else
            "main pipeline")
    print(f"  Mode:     {mode}")
    print(f"  Screen:   {args.screen}")
    if args.ablation:
        print(f"  Param:    {SCREENS[args.screen]['abl_param_mode']}"
              + (f"  | plate filter: {args.plate_filter}" if args.plate_filter else ""))
    print(f"  Image:    {IMAGE_URI}")
    print(f"  Instance: {args.instance_type}")
    print(f"  Output:   {out}")
    print(f"  Console:  {console}")
    print("\nThe job now runs independently of this session. You can close the "
          "window or let it time out — it will not stop the job.")

    if args.wait:
        print("\nPolling status (Ctrl-C is safe — job keeps running)...")
        while True:
            d = sm.describe_processing_job(ProcessingJobName=name)
            status = d["ProcessingJobStatus"]
            print(f"  {time.strftime('%H:%M:%S')}  {status}")
            if status in ("Completed", "Failed", "Stopped"):
                if status == "Failed":
                    print("  Reason:", d.get("FailureReason", "(none)"))
                break
            time.sleep(30)


if __name__ == "__main__":
    main()
