# Finomena zebrafish rescue drug screen — analysis pipeline

Analysis code for a larval-zebrafish behavioural **rescue drug screen**. Larvae
carrying a disease (knockout, KO) allele and wild-type (WT) siblings are recorded
in a video-tracking box; locomotor activity ("pixel difference") is modelled over
time to ask which compounds move KO larvae back toward the WT phenotype.

The statistical core is a **Tweedie generalized additive mixed model (GAMM)** with
an AR(1) error (`mgcv::bam`), fit per plate, from which three per-compound effects
are derived — movement off the disease baseline, landing relative to wild-type, and
a Bayesian posterior probability of equivalence to wild-type. See
[`docs/SAGEMAKER_README.md`](docs/SAGEMAKER_README.md) for running on AWS and
[`docs/BATCH_CONFOUND_WRITEUP.md`](docs/BATCH_CONFOUND_WRITEUP.md) for the design
reasoning.

## Repository layout

```
README.md  LICENSE  CITATION.cff  requirements.txt  R_dependencies.txt   # root

src/                          all analysis + orchestration code
  initial_screen_script.r     Step 1 — the BAM fit (per plate). Main analysis script.
  consolidate_results.py      Step 2 — per-plate results → one wide table.
  compute_ssmd_glass.py       Step 3 — model-free robust effect sizes (Glass's Δ, SSMD).
  rescue_addons.r             Step 4 — movement-confidence + model-vs-model-free check.
  make_plots.r                Step 5 — publication figures (ggplot2).
  run_pipeline.sh             Runs steps 1–5 in order (the in-container entrypoint).
  stand-alone_script.r        Balanced-design BAM script (genotype × drug factorial).
  ablation_study.R            Formula-selection study (why the production formula).
  launch_bam_job.py           Submit the pipeline as an AWS SageMaker Processing job.
  sync_to_s3.py               Upload code + a screen's inputs to S3.
  watch_job.sh                Poll a running job's status.
  make_plots.py               DEPRECATED matplotlib figures (superseded by make_plots.r).
  helpers/scripts/            Shared R: multiple-testing corrections, tree metadata,
                              the original app model script, family selection.

docs/
  MANIFEST.md                 Full file-by-file manifest + data-bundle index.
  SAGEMAKER_README.md         How the AWS job is set up and run.
  BATCH_CONFOUND_WRITEUP.md   The per-plate / batch-confound rationale.
  INITIAL_SCREEN_BUILD_PLAN.md  Design notes for the initial-screen build.
```

## Data (not in this repo)

The raw tracking data, intermediates, and outputs (~2.3 GB) are **not** committed;
`Data/` is git-ignored. The scripts expect a `Data/` tree in the repo root:

```
Data/<Screen>/Stand-Alone Inputs/     design.csv, time_frame.csv, per-plate layouts + .xls
Data/<Screen>/Production Output/       written by the pipeline
```

where `<Screen>` is `Initial screen` or `Confirmation`. Each input plate folder
holds a `plate_layout.csv` and the Zebrabox `.xls` recordings. See
[`docs/MANIFEST.md`](docs/MANIFEST.md) for the archived data bundles (inputs,
intermediates, results) and their contents.
**[Add the download location for the data bundles here — e.g. a Zenodo/figshare DOI
or institutional URL.]**

## Dependencies

- **R ≥ 4.4** with the packages in [`R_dependencies.txt`](R_dependencies.txt). The
  scripts auto-install anything missing at runtime, so no manual step is required;
  install manually only to pre-build an environment or a container image.
- **Python ≥ 3.10** — see [`requirements.txt`](requirements.txt). The two
  post-processing steps use only the standard library; `boto3` + the AWS CLI are
  needed only to run on SageMaker.

## Reproducing the analysis

The model is memory-heavy at full scale, so there are two paths.

### A. One plate / small dataset, locally

If you have R and a single screen's inputs under `Data/`, run the pipeline from the
repo root (the `FINOMENA_*_DIR` paths below are relative to where you invoke it):

```bash
# Initial-screen mode (per-plate rescue design, FDR correction):
FINOMENA_ANALYSIS_MODE=initial FINOMENA_ERROR_RATE=FDR \
FINOMENA_INPUT_DIR="Data/Initial screen/Stand-Alone Inputs" \
FINOMENA_OUTPUT_DIR="Data/Initial screen/Production Output" \
  bash src/run_pipeline.sh
```

`src/run_pipeline.sh` runs steps 1–5 and writes per-plate results, consolidated
tables, and figures under `FINOMENA_OUTPUT_DIR`. The balanced-design
(`confirmation`) script is `src/stand-alone_script.r`; set
`FINOMENA_ANALYSIS_MODE=confirmation` and `FINOMENA_ERROR_RATE=FWER`.

### B. Full screen on AWS SageMaker

The full initial screen (14 plates, ~2,800 animals) exceeds a normal workstation's
memory and is run as a SageMaker Processing job. This requires an AWS account with
SageMaker access, an ECR image (see [`docs/SAGEMAKER_README.md`](docs/SAGEMAKER_README.md)),
and these environment variables (so no account-specific values are committed):

```bash
export FINOMENA_AWS_ACCOUNT=<your 12-digit account id>
export FINOMENA_AWS_REGION=us-east-1            # optional (default)
# optional overrides: FINOMENA_S3_BUCKET, FINOMENA_ROLE_ARN, FINOMENA_IMAGE_URI

python src/sync_to_s3.py       --screen initial     # upload code + inputs to S3
python src/launch_bam_job.py   --screen initial --wait
```

See [`docs/SAGEMAKER_README.md`](docs/SAGEMAKER_README.md) for the one-time
image/role setup, instance sizing, and how to pull results back.

## License and citation

Released under the [MIT License](LICENSE) — free to use, modify, and redistribute,
provided the copyright notice is retained. If you use this code, please cite the
project and the associated paper as described in [`CITATION.cff`](CITATION.cff).
**[Fill in the author, paper, and repository details marked TODO in `CITATION.cff`
and the copyright holder in `LICENSE` before publishing.]**
