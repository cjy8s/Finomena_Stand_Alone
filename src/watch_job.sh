#!/usr/bin/env bash
# Check a finomena SageMaker Processing job's STATUS only — done or not, and
# whether it succeeded. No CloudWatch logs.
#
#   ./watch_job.sh                 # newest finomena-*-initial-* job, follow until done
#   ./watch_job.sh <job-name>      # a specific job, follow until done
#   ./watch_job.sh --once          # print status once and exit (no waiting)
#   ./watch_job.sh <job-name> --once
#
# Follow mode prints a single status line every 30s until the job ends, then a
# final DONE / FAILED line. Ctrl-C is safe — the job keeps running in AWS.
set -uo pipefail
REGION="us-east-1"

JOB=""; ONCE=0
for a in "$@"; do
  case "$a" in
    --once) ONCE=1 ;;
    *) JOB="$a" ;;
  esac
done

# Auto-pick the newest initial-screen job (matches both finomena-bam-initial and
# finomena-post-initial) if none given.
if [ -z "$JOB" ]; then
  JOB=$(aws sagemaker list-processing-jobs --region "$REGION" \
        --name-contains "initial" --sort-by CreationTime --sort-order Descending \
        --max-results 1 --query 'ProcessingJobSummaries[0].ProcessingJobName' \
        --output text 2>/dev/null)
  if [ -z "$JOB" ] || [ "$JOB" = "None" ]; then
    echo "No *initial* job found. Pass a job name explicitly."; exit 1
  fi
  echo "Job: $JOB"
fi

get_status() {
  aws sagemaker describe-processing-job --region "$REGION" \
    --processing-job-name "$JOB" --query 'ProcessingJobStatus' --output text 2>/dev/null
}

report_final() {
  local st="$1"
  if [ "$st" = "Completed" ]; then
    echo "✅ DONE — $JOB completed successfully."
  else
    echo "❌ $st — $JOB did not complete."
    # Only on failure, surface the one-line reason (still no full log).
    aws sagemaker describe-processing-job --region "$REGION" \
      --processing-job-name "$JOB" \
      --query 'FailureReason' --output text 2>/dev/null | sed 's/^/   reason: /'
  fi
}

ST=$(get_status)
if [ "$ONCE" = "1" ]; then
  case "$ST" in
    Completed|Failed|Stopped) report_final "$ST" ;;
    *) echo "⏳ $ST — still running." ;;
  esac
  exit 0
fi

# Follow: one status line per poll until terminal.
while true; do
  case "$ST" in
    Completed|Failed|Stopped) report_final "$ST"; break ;;
    *) echo "$(date '+%H:%M:%S')  ⏳ $ST" ;;
  esac
  sleep 30
  ST=$(get_status)
done
