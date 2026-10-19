#!/usr/bin/env bash
# =============================================================================
# Initial-screen pipeline — single in-container sequencer.
# -----------------------------------------------------------------------------
# Runs the whole rescue screen end-to-end in ONE execution, after a single BAM
# fit, in the order each step's inputs require:
#
#   1. initial_screen_script.r   BAM fit (per plate) -> per-plate master_results,
#                                posterior_equivalence_summary, aggregated_data,
#                                bam_residual_df.csv
#   2. consolidate_results.py    per-plate long  -> ALL_PLATES_full_results.csv
#   3. compute_ssmd_glass.py     aggregated_data -> SSMD_vs_model_ALL.csv (robust
#                                Glass's delta; read-only effect sizes)
#   4. rescue_addons.r           BAM movement confidence + Check A (Glass vs BAM)
#                                -> augmented_per_plate.csv, check_A_glass_vs_bam.csv
#   5. make_plots.r              figures (PNGs) + power_diagnostic_volcano.csv (ggplot2)
#
# Steps 2 and 3 are independent of each other (both read step 1's output) but
# both must follow step 1; steps 4 and 5 read 2+3. The chain is strictly
# sequential here for simplicity and legible logs — the BAM fit dominates wall
# time, so serializing the fast post-processing costs nothing.
#
# Designed to run as the ContainerEntrypoint of the SageMaker Processing job
# (launch_bam_job.py), where SageMaker mounts code at /opt/ml/processing/code,
# inputs at /opt/ml/processing/input, and output at /opt/ml/processing/output.
# All paths are env-driven with local fallbacks so it also runs by hand:
#
#   ./run_pipeline.sh            # local: uses repo-relative defaults
#   (container)                  # env vars set by launch_bam_job.py
#
# `set -e` aborts the chain on the first failing step — addons never run on a
# stale or partial BAM result.
# =============================================================================
set -euo pipefail

# ── Locate the code dir (this script's location) so step calls are absolute ──
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Resolve I/O dirs (container mounts -> local fallback) ────────────────────
# initial_screen_script.r reads FINOMENA_INPUT_DIR / FINOMENA_OUTPUT_DIR; the
# post-processing scripts read FINOMENA_OUTPUT_DIR (per-plate folders live under
# it) and write the consolidated tables + plots there.
: "${FINOMENA_INPUT_DIR:=Data/Initial screen/Stand-Alone Inputs}"
: "${FINOMENA_OUTPUT_DIR:=Data/Initial screen/Production Output}"
: "${FINOMENA_HELPERS_DIR:=helpers/scripts}"
: "${FINOMENA_ERROR_RATE:=FDR}"          # initial screen = BH/FDR
: "${FINOMENA_ANALYSIS_MODE:=initial}"
# FINOMENA_PLATES is intentionally NOT defaulted here. The BAM auto-discovers
# plates from the input data (one output folder per plate it fits), so the plate
# list is DERIVED from that output below (derive_plates) — the single source of
# truth — right after the folders exist. This makes post-processing cover every
# plate the BAM ran, automatically, with no hardcoded count to drift out of date.
# Set FINOMENA_PLATES explicitly only to force a subset; a non-empty value here
# is honoured verbatim and skips derivation.
# Post-processing scripts read these (default to the BAM output dir).
: "${FINOMENA_ADDON_PROD_DIR:=$FINOMENA_OUTPUT_DIR}"
: "${FINOMENA_ADDON_AGG_CSV:=$FINOMENA_OUTPUT_DIR/aggregated_data.csv}"
: "${FINOMENA_ADDON_OUTPUT_DIR:=$FINOMENA_OUTPUT_DIR}"
: "${FINOMENA_PLOTS_DIR:=$FINOMENA_OUTPUT_DIR/Plots}"
# Trajectory panel reads the fresh aggregated_data.csv directly (it now carries a
# bio_plate column), so no /tmp/agg_full.csv precursor is needed.
: "${FINOMENA_PLOT_AGG_CSV:=$FINOMENA_OUTPUT_DIR/aggregated_data.csv}"

# FINOMENA_PLATES is exported separately: it may be derived from the BAM output
# below (see ensure_plates), so it is exported at that point rather than here.
export FINOMENA_INPUT_DIR FINOMENA_OUTPUT_DIR FINOMENA_HELPERS_DIR \
       FINOMENA_ERROR_RATE FINOMENA_ANALYSIS_MODE \
       FINOMENA_ADDON_PROD_DIR FINOMENA_ADDON_AGG_CSV FINOMENA_ADDON_OUTPUT_DIR \
       FINOMENA_PLOTS_DIR FINOMENA_PLOT_AGG_CSV

# Pick interpreters (container may expose python3 only; SageMaker R image has Rscript).
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON="${PYTHON:-python3}"

step() { echo; echo "========== [$(date '+%H:%M:%S')] STEP $1: $2 =========="; }

# Derive the plate list from the per-plate BAM output folders — each plate the
# BAM fit produces holds a master_results.csv, so the folders ARE the ground
# truth for which plates exist. `sort -V` (natural sort) keeps "Plate 10" after
# "Plate 9", not after "Plate 1". Mirrors the same "read plates from the data"
# idiom rescue_addons.r already uses (sort(unique(full_df$Plate))).
derive_plates() {
  local d out=""
  while IFS= read -r d; do
    out="${out:+$out,}$(basename "$d")"
  done < <(find "$FINOMENA_OUTPUT_DIR" -mindepth 2 -maxdepth 2 \
                -name master_results.csv -printf '%h\n' 2>/dev/null | sort -V)
  printf '%s' "$out"
}

# ── RESUME MODES ─────────────────────────────────────────────────────────────
# Two ways to reuse a prior run's committed BAM output instead of re-fitting the
# multi-hour model. Both seed the prior <plate>/ folders + aggregated_data.csv +
# bam_residual_df.csv from FINOMENA_RESUME_DIR (mounted read-only) into
# FINOMENA_OUTPUT_DIR, differing only in which downstream steps re-run:
#
#   FINOMENA_POST_ONLY=1  skip STEPs 1-3 (BAM + consolidate + SSMD), run 4-5.
#                         Use when ONLY the addon/plot code (STEP 4/5) changed;
#                         the seeded consolidate/SSMD tables are reused as-is.
#   FINOMENA_SKIP_BAM=1   skip ONLY the BAM fit (STEP 1), re-run STEPs 2-5.
#                         Use when the fit is fine but the consolidated tables
#                         must be rebuilt — e.g. plate coverage changed, or the
#                         consolidate/SSMD code changed. This is the correct mode
#                         to (re)process every plate the BAM actually fit.
POST_ONLY="${FINOMENA_POST_ONLY:-0}"
SKIP_BAM="${FINOMENA_SKIP_BAM:-0}"
if [ "$POST_ONLY" = "1" ] || [ "$SKIP_BAM" = "1" ]; then
  : "${FINOMENA_RESUME_DIR:?resume modes (POST_ONLY/SKIP_BAM) require FINOMENA_RESUME_DIR}"
  mkdir -p "$FINOMENA_OUTPUT_DIR"
  echo "RESUME: seeding prior BAM outputs from $FINOMENA_RESUME_DIR -> $FINOMENA_OUTPUT_DIR"
  cp -r "$FINOMENA_RESUME_DIR"/. "$FINOMENA_OUTPUT_DIR"/
fi

# ── PARSE-GATE: fail in seconds on a syntax error, not after the hours-long fit ─
# parse() each R script (and py_compile each Python step) WITHOUT executing, so a
# typo aborts the job immediately instead of wasting the BAM fit. R is not
# available in every dev sandbox, so this is the first place the R files are
# actually parsed by an interpreter.
echo "========== [$(date '+%H:%M:%S')] STEP 0: parse-gate (syntax check all steps) =========="
"$RSCRIPT" -e 'invisible(lapply(c("initial_screen_script.r","rescue_addons.r","make_plots.r"), function(f) { parse(f); cat("  R parse OK:", f, "\n") }))'
for f in consolidate_results.py compute_ssmd_glass.py; do
  "$PYTHON" -m py_compile "$f" && echo "  py compile OK: $f"
done

# ── STEP 1: BAM fit (skipped in both resume modes) ───────────────────────────
if [ "$POST_ONLY" = "1" ] || [ "$SKIP_BAM" = "1" ]; then
  echo; echo "RESUME: skipping STEP 1 (BAM fit); per-plate folders come from the seeded output."
else
  step 1 "BAM fit (initial_screen_script.r)"
  "$RSCRIPT" initial_screen_script.r
fi

# ── DERIVE PLATE LIST (now that the per-plate folders exist) ──────────────────
# The BAM fit (or the seeded resume output) has produced one <plate>/ folder per
# plate. Derive FINOMENA_PLATES from those folders unless the caller pinned it to
# force a subset. Every downstream step (consolidate, SSMD, addons, plots) reads
# FINOMENA_PLATES, so setting it once here makes full coverage automatic.
if [ -n "${FINOMENA_PLATES:-}" ]; then
  echo "Plate list PINNED by caller: FINOMENA_PLATES=$FINOMENA_PLATES"
else
  FINOMENA_PLATES="$(derive_plates)"
  [ -n "$FINOMENA_PLATES" ] || { echo "ERROR: no per-plate master_results.csv found under $FINOMENA_OUTPUT_DIR — cannot derive plate list." >&2; exit 1; }
  echo "Derived plate list from BAM output: FINOMENA_PLATES=$FINOMENA_PLATES"
fi
export FINOMENA_PLATES

# ── STEPs 2-3: consolidate + SSMD (skipped only in POST_ONLY) ─────────────────
if [ "$POST_ONLY" = "1" ]; then
  echo; echo "POST-ONLY: skipping STEPs 2-3 (consolidate + SSMD); using seeded tables."
else
  step 2 "Consolidate per-plate -> ALL_PLATES_full_results.csv"
  "$PYTHON" consolidate_results.py

  step 3 "Robust Glass's delta / SSMD -> SSMD_vs_model_ALL.csv"
  "$PYTHON" compute_ssmd_glass.py
fi

step 4 "BAM movement confidence + Check A (rescue_addons.r)"
"$RSCRIPT" rescue_addons.r

step 5 "Plots + power diagnostic (make_plots.r)"
# Plots are rendered in R/ggplot2 (the R image has no pip/matplotlib). Best-effort:
# the analytical outputs (augmented table, Check A/B) are already written by STEP
# 4, so a plot failure must NOT fail the whole job — warn and continue.
"$RSCRIPT" make_plots.r || echo "  WARNING: make_plots.r failed; analytical outputs from STEP 4 are unaffected."

echo
echo "========== [$(date '+%H:%M:%S')] PIPELINE COMPLETE =========="
echo "Outputs under: $FINOMENA_OUTPUT_DIR"
