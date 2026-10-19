#!/usr/bin/env python
"""
Upload code + helpers + a screen's inputs to the S3 layout launch_bam_job.py
expects. Run this whenever you change the R script, the helpers, or the data.

    python sync_to_s3.py                  # confirmation screen (default)
    python sync_to_s3.py --screen initial # the larger initial screen

S3 layout produced (per screen):
    s3://.../finomena/bam-run/<screen>/code/stand-alone_script.r
    s3://.../finomena/bam-run/<screen>/code/helpers/scripts/*.R
    s3://.../finomena/bam-run/<screen>/input/...   (= the screen's Stand-Alone Inputs)
"""
import argparse
import os
import subprocess
import sys

# Account-specific; read from the environment so nothing account-specific is
# committed. Set FINOMENA_S3_BUCKET (or FINOMENA_AWS_ACCOUNT, from which the
# default bucket name sagemaker-<region>-<account> is derived) before syncing.
# Matches the resolution used by launch_bam_job.py.
BUCKET = os.environ.get("FINOMENA_S3_BUCKET")
if not BUCKET:
    _acct = os.environ.get("FINOMENA_AWS_ACCOUNT")
    _region = os.environ.get("FINOMENA_AWS_REGION", "us-east-1")
    if not _acct:
        sys.exit("Set FINOMENA_S3_BUCKET, or FINOMENA_AWS_ACCOUNT to derive "
                 "sagemaker-<region>-<account>, before syncing.")
    BUCKET = f"sagemaker-{_region}-{_acct}"
PREFIX = "finomena/bam-run"

# Path anchors (robust to the working directory). Code lives alongside this
# script in src/; data lives in Data/ at the repo root (src/'s parent).
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(CODE_DIR)
SCREENS = {
    "confirmation": os.path.join(REPO_ROOT, "Data", "Confirmation"),
    "initial": os.path.join(REPO_ROOT, "Data", "Initial screen"),
}
SCREEN_S3 = {"confirmation": "confirmation", "initial": "initial"}


def run(cmd):
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def code(path):
    """A code file/dir path resolved against src/ (this script's directory)."""
    return os.path.join(CODE_DIR, path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen", choices=sorted(SCREENS), default="confirmation")
    args = p.parse_args()

    local = SCREENS[args.screen]
    base = f"s3://{BUCKET}/{PREFIX}/{SCREEN_S3[args.screen]}"
    inputs = f"{local}/Stand-Alone Inputs"

    # Code shared by every screen: the confirmation/standalone R script, the
    # ablation script, and the helper scripts. Local paths resolve against src/
    # (code()); the S3 layout stays flat under code/ so the container sees every
    # script as a sibling regardless of the local src/ grouping.
    run(["aws", "s3", "cp", code("stand-alone_script.r"), f"{base}/code/stand-alone_script.r"])
    run(["aws", "s3", "cp", code("ablation_study.R"), f"{base}/code/ablation_study.R"])
    run(["aws", "s3", "cp", code("helpers/scripts/"), f"{base}/code/helpers/scripts/",
         "--recursive", "--exclude", "*", "--include", "*.R"])

    # The initial screen runs the FULL pipeline via run_pipeline.sh, so it needs
    # the wrapper + the dedicated initial BAM script + the post-processing steps.
    # (launch_bam_job.py points the initial entrypoint at run_pipeline.sh, which
    # calls these in order.) Without these uploaded, the wrapper would abort at
    # step 1 with "initial_screen_script.r not found" in the container.
    if args.screen == "initial":
        for f in ("run_pipeline.sh",
                  "initial_screen_script.r",
                  "consolidate_results.py",
                  "compute_ssmd_glass.py",
                  "rescue_addons.r",
                  "make_plots.r"):
            run(["aws", "s3", "cp", code(f), f"{base}/code/{f}"])

    # Inputs: this screen's Stand-Alone Inputs (delete so removed files don't linger).
    run(["aws", "s3", "sync", inputs, f"{base}/input/", "--delete"])

    print(f"\nSynced screen '{args.screen}' to {base}/")
    print(f"Launch with:  python launch_bam_job.py --screen {args.screen}")


if __name__ == "__main__":
    sys.exit(main())
