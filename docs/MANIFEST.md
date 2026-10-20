# Finomena rescue-screen — file manifest

What's in this project, where it goes, and how to reproduce. Code goes to GitHub;
data (2.3 GB) is archived as separate tarballs (see **Download bundles** below).

All code lives under `src/` (with shared R in `src/helpers/scripts/`); the narrative
docs live under `docs/`. Script names below are relative to `src/`.

## Pipeline (order of execution — see `src/run_pipeline.sh`)

| Step | Script | Role |
|------|--------|------|
| 1 | `initial_screen_script.r` | BAM fit per plate → per-plate `master_results`, posterior draws, `aggregated_data.csv` |
| 2 | `consolidate_results.py`  | per-plate long → `ALL_PLATES_full_results.csv` |
| 3 | `compute_ssmd_glass.py`   | robust Glass's δ / SSMD → `SSMD_vs_model_ALL.csv` |
| 4 | `rescue_addons.r`         | BAM movement confidence + Check A → `augmented_per_plate.csv`, `check_A_glass_vs_bam.csv` |
| 5 | `make_plots.r`            | 10 figures + `power_diagnostic_volcano.csv` |

`stand-alone_script.r` is the confirmation-screen BAM script; `ablation_study.R`
is the formula-selection study. `helpers/scripts/` holds shared R (`corrections.R`,
`tree_metadata.R`, `TweedieAR1_BAM.R`, `FamilySelection.R`).

## Orchestration (SageMaker Processing job — R is not run locally)

- `launch_bam_job.py` — submits the job (boto3 `create_processing_job`).
- `sync_to_s3.py` — uploads code + a screen's inputs to S3.
- `watch_job.sh` — polls job status.
- `run_pipeline.sh` — the in-container entrypoint that runs steps 1–5.

Regenerate figures only (no multi-hour refit), then pull results back (run from the
repo root; scripts live in `src/`):
```bash
python src/sync_to_s3.py --screen initial
python src/launch_bam_job.py --screen initial --post-only --wait
aws s3 sync s3://sagemaker-us-east-1-<ACCOUNT_ID>/finomena/bam-run/initial/output/Plots/ "Data/Initial screen/Plots/"
```

## Docs
`README.md` (repo landing page / quick start), `SAGEMAKER_README.md`,
`BATCH_CONFOUND_WRITEUP.md`, `INITIAL_SCREEN_BUILD_PLAN.md`. Dependencies:
`Python_dependencies.txt` (Python), `R_dependencies.txt` (R). Licensing / citation:
`LICENSE` (MIT), `CITATION.cff`.

## Deprecated
`make_plots.py` — the old matplotlib figure script, superseded by `make_plots.r`
(the SageMaker R image has no matplotlib). Kept for reference only; not in the pipeline.

## Download bundles (data — NOT in git)

| Tarball | Contents | ~Size |
|---------|----------|-------|
| `finomena_code.tgz` | All scripts, helpers, docs, `.gitignore`, this manifest (= the GitHub repo) | <1 MB |
| `finomena_benchmarks.tgz` | Start-of-project BAM formula ablation: `ablation_study.R` + per-variant metrics/CV + `BENCHMARKS_README.md` | ~25 KB |
| `finomena_results.tgz` | Analysis outputs + figures, both screens (small, high-value; ablation lives in the benchmarks bundle) | ~23 MB |
| `finomena_inputs.tgz` | Raw `.xls` tracking data + layouts, both screens (reproduce from scratch) | ~113 MB |
| `finomena_intermediates.tgz` | `aggregated_data.csv` + per-plate posterior draws (skip the BAM refit) | ~218 MB |

The **benchmarks** bundle is the formula-selection study run at the very start (why the
production formula was chosen — it OOMs the initial screen, so a cheaper functionally-
equivalent formula was benchmarked). See `BENCHMARKS_README.md` inside it for details.
