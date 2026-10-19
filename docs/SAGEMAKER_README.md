# Running the BAM R script on SageMaker (detached from your session)

This runs `stand-alone_script.r` as a **SageMaker Processing job**: it executes
in a separate AWS-managed container, *not* in your notebook/Studio session. So
**closing the window or hitting the admin session timeout does NOT stop it** —
the job runs to completion on AWS and writes results to S3.

This is the key difference from running `Rscript stand-alone_script.r` in a
terminal: that process dies with your session; a Processing job does not.

## One-time setup (already done on 2026-06-03)

1. **R container image in ECR.** SageMaker pulls the job's image from ECR.
   We mirrored the public `rocker/tidyverse:4.4.2` image (R + tidyverse +
   data.table) into ECR with `crane` (no local Docker needed):
   `<ACCOUNT_ID>.dkr.ecr.us-east-1.amazonaws.com/finomena-r:4.4.2`
   The script `install.packages()` any extras (emmeans, gratia, jsonlite) at
   runtime — the job has internet egress, so CRAN works.

2. **Code + inputs in S3**, under
   `s3://sagemaker-us-east-1-<ACCOUNT_ID>/finomena/bam-run/`:
   - `code/stand-alone_script.r` + `code/helpers/`
   - `input/`  (the contents of `Stand-Alone Inputs/`)
   - `output/` (results land here)

## Run a job

```bash
cd <repo-root>                         # scripts live under src/
python src/launch_bam_job.py            # submit, detached — prints job name + console link
python src/launch_bam_job.py --wait     # also stream status here (Ctrl-C is safe; job keeps running)
python src/launch_bam_job.py --instance-type ml.m7i.4xlarge   # different instance
```

The script's `INPUT_DIR` / `OUTPUT_DIR` are overridden via the
`FINOMENA_INPUT_DIR` / `FINOMENA_OUTPUT_DIR` env vars to the container mount
points. With those unset (running locally), it falls back to the original
`Stand-Alone Inputs` / `Stand-Alone Output` — local behaviour is unchanged.

## If you change the script, helpers, or inputs

Re-upload before launching:

```bash
B=s3://sagemaker-us-east-1-<ACCOUNT_ID>/finomena/bam-run
aws s3 cp src/stand-alone_script.r "$B/code/stand-alone_script.r"
aws s3 cp src/helpers/ "$B/code/helpers/" --recursive
aws s3 cp "Data/Confirmation/Stand-Alone Inputs/" "$B/input/" --recursive   # only if inputs changed
```

## Monitor a running job

- **Console:** the launcher prints a direct link.
- **CLI status:**
  ```bash
  aws sagemaker describe-processing-job --processing-job-name <JOB_NAME> \
    --query '[ProcessingJobStatus, FailureReason]' --output text
  ```
- **Live logs (CloudWatch):** log group `/aws/sagemaker/ProcessingJobs`,
  stream `<JOB_NAME>/...`. In the console: the job page → "Monitoring" →
  "View logs". This is where the R script's `cat()` progress output appears.

## Get the results

```bash
aws s3 cp s3://sagemaker-us-east-1-<ACCOUNT_ID>/finomena/bam-run/output/ \
  "./Stand-Alone Output/" --recursive
```
Produces `master_results.csv`, `summary_statistics_by_group.csv`, and
`aggregated_data.csv`.

## Notes / tuning

- **Instance:** default `ml.r7i.48xlarge` (192 vCPU / 1536 GB) — this is a
  large, memory-heavy job. The script's own estimate is ~35 GB design matrix /
  ~140 GB peak (158K rows × 27.7K cols), so smaller boxes OOM: a 4xlarge
  (128 GB) is far too small. `ml.r5.*` has **0 quota** on this account;
  `m7i` / `r7i` / `c7i` families work. Request quota increases in Service
  Quotas if you need a different size.
- **Package install:** the launcher sets `CRAN_MIRROR` to the Posit binary repo
  so `install.packages()` pulls precompiled binaries (the container has no
  cmake for source builds). `gratia` was removed from the script's required
  list — it is unused in the standalone path and failed to load.
- **Runtime cap:** `--max-runtime-h` (default 12 h). SageMaker stops the job at
  this limit, so set it above your expected runtime.
- **Cost:** you pay only for job wall-clock; the instance is torn down when the
  script exits. No idle notebook cost.
- **Speed-up option:** to skip the per-run `install.packages()` step, build a
  custom image with all packages baked in and push to the same ECR repo (needs
  a Docker daemon or CodeBuild access — neither is available on this instance
  today).
