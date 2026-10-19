# =============================================================================
# Rescue-Screen Addons  —  BAM movement confidence + Check A (Glass vs BAM)
# -----------------------------------------------------------------------------
# A standalone POST-PROCESSING reader. It does NOT re-fit the BAM, recompute the
# robust Glass's delta, or touch the Bayesian functional posterior-equivalence
# test. It reads the frozen outputs of initial_screen_script.r + the
# consolidation/SSMD steps and adds, PER PLATE (never pooled across plates):
#
#   Addition 1 — movement-confidence on the BAM KO+drug vs KO+DMSO contrast.
#                Default: shrinkage-free one-sided sign confidence from the BAM's
#                own (estimate, SE, residual df) — sound because the Tweedie
#                dispersion makes the SE correct. Optional ashr shrinkage via the
#                FINOMENA_USE_ASHR toggle (see note below).
#   Check A    — Spearman of robust Glass's delta vs BAM movement, plus a 2x2
#                set-agreement (counts, b+c, Jaccard, SMC) over a threshold grid.
#
# NOTE: an independent limma engine + "Check B" were trialled and REMOVED — on
# this 2-3 animal/drug screen limma only re-confirmed the BAM (same suppressor
# hits, no independent rescue signal), so it added no decision value. See
# [[finomena-ashr-decision]] for the related ashr call.
#
# RESCUE LOGIC (read downstream, not decided here): a hit = confident MOVEMENT
# toward WT AND high LANDING equivalence to WT (the EXISTING P(equivalent to WT),
# read as-is). The landing contrast is NEVER shrunk (shrinkage biases the
# equivalence call).
#
# Inputs (read-only):
#   ALL_PLATES_full_results.csv  movement DrugEffect_vsKO_*, landing Rescue_vsWT_*,
#                                P(equivalent to WT) = equiv_Pr_vsWT
#   SSMD_vs_model_ALL.csv        robust Glass's delta (Glass_vsKO_log = movement axis)
#   bam_residual_df.csv          BAM residual df per plate (additive side-output of
#                                initial_screen_script.r; sets the t df)
#
# Outputs (written to OUTPUT_DIR):
#   augmented_per_plate.csv      one row per drug per plate
#   check_A_glass_vs_bam.csv     Spearman + 2x2 grid
# =============================================================================


# ── USER CONFIG ──────────────────────────────────────────────────────────────
INPUT_PROD_DIR <- Sys.getenv(
  "FINOMENA_ADDON_PROD_DIR",
  unset = "Data/Initial screen/Production Output")
OUTPUT_DIR <- Sys.getenv(
  "FINOMENA_ADDON_OUTPUT_DIR",
  unset = INPUT_PROD_DIR)                                    # write alongside inputs

# Control conditions (must match initial_screen_script.r's POS/NEG controls).
POS_CONTROL <- Sys.getenv("FINOMENA_POS_CONTROL", unset = "KO+DMSO")  # KO+DMSO
NEG_CONTROL <- Sys.getenv("FINOMENA_NEG_CONTROL", unset = "WT+DMSO")   # WT+DMSO

# ashr toggle (default OFF). ashr shrinks movement effects across the ~80 drugs
# and yields an s-value, but in this 2-3 animal/drug screen it flags NO confident
# positive movers and barely changes the ranking — see [[finomena-ashr-decision]].
# When OFF we instead report a SHRINKAGE-FREE sign-confidence from the BAM's own
# (estimate, SE, residual df): the one-sided tail  P(sign wrong) = pt(-|z|, df)
# on the toward-WT side. Set to 1/true/yes to re-enable ashr shrinkage.
USE_ASHR <- tolower(Sys.getenv("FINOMENA_USE_ASHR", unset = "0")) %in%
            c("1", "true", "yes", "on")

# Gap-anchored grid for the MOVEMENT binary "moved toward WT" call (% recovery
# toward WT). Used ONLY for the supplementary 2x2 set-agreement; the Spearman
# correlation is the primary, threshold-free metric. cutoff = f * disease_gap.
MOVE_FRACTION_GRID <- as.numeric(strsplit(
  Sys.getenv("FINOMENA_MOVE_FRAC_GRID", unset = "0.15,0.20,0.25,0.33"),
  ",", fixed = TRUE)[[1]])

# Glass's delta is in robust-SD units, NOT log2FC, so it cannot share the model's
# gap-anchored cutoff. Binarize it by its own effect-size grid (Cohen-style).
GLASS_CUTOFF_GRID <- as.numeric(strsplit(
  Sys.getenv("FINOMENA_GLASS_CUTOFF_GRID", unset = "0.2,0.5,0.8"),
  ",", fixed = TRUE)[[1]])


# ── PACKAGE LIBRARY LOCATION (mirror initial_screen_script.r) ────────────────
.conda_prefix <- Sys.getenv("CONDA_PREFIX")
.lib <- if (nzchar(.conda_prefix)) {
  file.path(.conda_prefix, "lib", "R", "library")
} else {
  file.path(Sys.getenv("HOME", unset = "/tmp"), "R", "library")
}
dir.create(.lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(.lib, .libPaths()))

cran_mirror <- Sys.getenv("CRAN_MIRROR", unset = "https://cloud.r-project.org")

# CRAN packages. ashr only needed when the shrinkage toggle is on. (limma was
# removed along with the independent engine.)
.cran_pkgs <- c("data.table", "dplyr", "readr")
if (USE_ASHR) .cran_pkgs <- c(.cran_pkgs, "ashr")
.missing <- .cran_pkgs[!vapply(.cran_pkgs, requireNamespace, logical(1),
                               quietly = TRUE)]
if (length(.missing) > 0) {
  cat("Installing CRAN packages: ", paste(.missing, collapse = ", "), "\n",
      sep = "")
  install.packages(.missing, repos = cran_mirror)
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  if (USE_ASHR) library(ashr)
})

options(width = 320)
options(warn = 1)
`%||%` <- function(a, b) if (is.null(a)) b else a


# ── ANCHOR RELATIVE PATHS TO THIS SCRIPT'S LOCATION ──────────────────────────
.get_script_dir <- function() {
  args0 <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args0, value = TRUE)
  if (length(fa) > 0) return(normalizePath(dirname(sub("^--file=", "", fa[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (nzchar(p)) return(normalizePath(dirname(p)))
  }
  normalizePath(getwd())
}
SCRIPT_DIR <- .get_script_dir()
.is_abs <- function(p) grepl("^([A-Za-z]:)?[/\\\\]", p)
.resolve <- function(p) {
  if (length(p) == 0 || !nzchar(p)) return(p)
  if (.is_abs(p)) p else normalizePath(file.path(SCRIPT_DIR, p),
                                       mustWork = FALSE)
}
INPUT_PROD_DIR <- .resolve(INPUT_PROD_DIR)
OUTPUT_DIR     <- .resolve(OUTPUT_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Production input dir: ", INPUT_PROD_DIR, "\n")
cat("Output dir:           ", OUTPUT_DIR, "\n")
cat("ashr shrinkage:       ", if (USE_ASHR) "ON" else "OFF (one-sided sign confidence)", "\n")
cat("Movement frac grid:   ", paste(MOVE_FRACTION_GRID, collapse = ", "), "\n")
cat("Glass cutoff grid:    ", paste(GLASS_CUTOFF_GRID, collapse = ", "), "\n\n")


# ── READ FROZEN INPUTS ───────────────────────────────────────────────────────
.req <- function(p, what) {
  if (!file.exists(p)) stop("Missing ", what, " at: ", p)
  p
}
full_path  <- .req(file.path(INPUT_PROD_DIR, "ALL_PLATES_full_results.csv"),
                   "ALL_PLATES_full_results.csv")
ssmd_path  <- .req(file.path(INPUT_PROD_DIR, "SSMD_vs_model_ALL.csv"),
                   "SSMD_vs_model_ALL.csv")

full_df <- as.data.table(read_csv(full_path, show_col_types = FALSE))
ssmd_df <- as.data.table(read_csv(ssmd_path, show_col_types = FALSE))

# bam_residual_df.csv is the additive side-output of initial_screen_script.r.
# Sets the t df for the sign-confidence (or ashr) likelihood. If absent, fall
# back to the normal approximation (df = Inf) and warn.
resid_df_path <- file.path(INPUT_PROD_DIR, "bam_residual_df.csv")
if (file.exists(resid_df_path)) {
  bam_df_tab <- as.data.table(read_csv(resid_df_path, show_col_types = FALSE))
  cat("BAM residual df table:\n"); print(bam_df_tab); cat("\n")
} else {
  bam_df_tab <- NULL
  warning("bam_residual_df.csv not found — movement sign confidence will use the ",
          "normal approximation (df = Inf). Re-run initial_screen_script.r to ",
          "emit the exact residual df.")
}

plate_list <- sort(unique(full_df$Plate))
cat("Plates: ", paste(plate_list, collapse = ", "), "\n\n")


# ── SMALL METRIC HELPERS ─────────────────────────────────────────────────────
# Spearman correlation, NA-safe (pairwise complete). Returns NA if <3 pairs.
.spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

# 2x2 set agreement between two binary "moved" calls. a = both moved, b = only
# the second moved, c = only the first moved, d = neither. Reports raw counts,
# disagreement b+c, Jaccard = a/(a+b+c), SMC = (a+d)/(a+b+c+d). No chance
# correction (kappa breaks under the heavy imbalance of a triage screen).
.agreement_2x2 <- function(call1, call2) {
  ok <- !is.na(call1) & !is.na(call2)
  c1 <- call1[ok]; c2 <- call2[ok]
  a <- sum( c1 &  c2); b <- sum(!c1 &  c2)
  c <- sum( c1 & !c2); d <- sum(!c1 & !c2)
  denom_j <- a + b + c
  data.frame(
    a = a, b = b, c = c, d = d,
    disagreement = b + c,
    jaccard = if (denom_j > 0) a / denom_j else NA_real_,
    smc     = (a + d) / max(a + b + c + d, 1L)
  )
}


# =============================================================================
# PER-PLATE LOOP
# =============================================================================
augmented_rows <- list()
checkA_rows    <- list()

for (.plate in plate_list) {
  cat(strrep("=", 70), "\n", "PLATE: ", .plate, "\n", strrep("=", 70), "\n",
      sep = "")

  fp <- full_df[Plate == .plate]
  sp <- ssmd_df[Plate == .plate]
  if (nrow(fp) == 0L) { cat("  no rows in ALL_PLATES for this plate; skip.\n"); next }

  # ── DISEASE GAP (read-only) ────────────────────────────────────────────────
  # KO->WT magnitude. For every drug, DrugEffect_vsKO - Rescue_vsWT algebraically
  # equals (WT+DMSO) - (KO+DMSO); the median across drugs denoises it. Positive
  # because WT moves more than KO. The gap-anchored movement cutoff is f * gap.
  gap <- median(fp$DrugEffect_vsKO_est_log2 - fp$Rescue_vsWT_est_log2,
                na.rm = TRUE)
  cat(sprintf("  Disease gap (KO->WT, log2FC): %.3f\n", gap))

  # ── ADDITION 1: BAM MOVEMENT confidence ─────────────────────────────────────
  bam_df_val <- if (!is.null(bam_df_tab)) {
    v <- bam_df_tab$residual_df[bam_df_tab$bio_plate == .plate]
    if (length(v) >= 1L && is.finite(v[1])) v[1] else NULL
  } else NULL

  if (USE_ASHR) {
    # ashr shrinkage across drugs (default mixcompdist = "uniform"; NOT "normal"
    # because ashr lacks a normal-mixture + student-t likelihood combination).
    cat("  ashr(BAM movement) df = ",
        if (is.null(bam_df_val)) "NULL (normal approx)" else sprintf("%.0f", bam_df_val),
        "\n", sep = "")
    ash_bam <- ash(betahat = fp$DrugEffect_vsKO_est_log2,
                   sebetahat = fp$DrugEffect_vsKO_SE, df = bam_df_val)
    fp[, `:=`(
      bam_move_conf_log2   = ash_bam$result$PosteriorMean,
      bam_move_svalue      = ash_bam$result$svalue,
      bam_move_lfsr        = ash_bam$result$lfsr
    )]
  } else {
    # Shrinkage-free sign confidence from the BAM's own (estimate, SE, df): the
    # one-sided tail toward WT. svalue = P(true sign is NOT toward WT). With no
    # shrinkage the "confidence" point estimate is just the raw effect.
    df_t <- if (is.null(bam_df_val)) Inf else bam_df_val
    z <- fp$DrugEffect_vsKO_est_log2 / fp$DrugEffect_vsKO_SE
    # toward WT = positive movement; tail probability the sign is wrong:
    sval_flat <- ifelse(z > 0, pt(-z, df = df_t), pt(z, df = df_t))
    cat("  ashr OFF — shrinkage-free one-sided sign confidence (df = ",
        if (is.infinite(df_t)) "normal" else sprintf("%.0f", df_t), ")\n", sep = "")
    fp[, `:=`(
      bam_move_conf_log2   = DrugEffect_vsKO_est_log2,   # raw (no shrinkage)
      bam_move_svalue      = sval_flat,
      bam_move_lfsr        = sval_flat
    )]
  }

  # ── ATTACH ROBUST GLASS'S DELTA (read-only, movement axis = vs KO) ──────────
  glass_cols <- intersect(c("drug", "Glass_vsKO_log", "Glass_vsWT_log"),
                          names(sp))
  fp <- merge(fp, sp[, ..glass_cols], by = "drug", all.x = TRUE, sort = FALSE)

  fp[, `:=`(Plate = .plate, disease_gap_log2 = gap)]
  augmented_rows[[.plate]] <- fp

  # ── CHECK A: Glass vs BAM movement ──────────────────────────────────────────
  # Primary: Spearman of Glass_vsKO_log against the BAM movement effect. (When
  # ashr is off, bam_move_conf_log2 == the raw effect, so the two Spearman
  # columns coincide; the column pair is kept so the ashr-on case still reports
  # both.)
  rA_raw  <- .spearman(fp$Glass_vsKO_log, fp$DrugEffect_vsKO_est_log2)
  rA_conf <- .spearman(fp$Glass_vsKO_log, fp$bam_move_conf_log2)
  cat(sprintf("  Check A Spearman(Glass, BAM) = %.3f%s\n", rA_raw,
              if (USE_ASHR) sprintf(" ; (Glass, BAM ashr) = %.3f", rA_conf) else ""))
  # Supplementary 2x2: model "moved" = movement >= f*gap (gap-anchored grid);
  # Glass "moved" = Glass_vsKO_log >= g (its own effect-size grid, different unit).
  for (f in MOVE_FRACTION_GRID) {
    model_call <- fp$DrugEffect_vsKO_est_log2 >= f * gap
    for (g in GLASS_CUTOFF_GRID) {
      glass_call <- fp$Glass_vsKO_log >= g
      ag <- .agreement_2x2(model_call, glass_call)
      checkA_rows[[length(checkA_rows) + 1L]] <- cbind(
        data.frame(Plate = .plate,
                   spearman = rA_raw,
                   move_fraction = f, move_cutoff_log2 = f * gap,
                   glass_cutoff = g, stringsAsFactors = FALSE),
        ag)
    }
  }
  cat("\n")
}


# =============================================================================
# ASSEMBLE + WRITE OUTPUTS
# =============================================================================
# Augmented per-plate table — movement (BAM raw + confidence; Glass), landing
# (BAM raw + existing P(equiv to WT)), and the disease gap.
augmented <- rbindlist(augmented_rows, use.names = TRUE, fill = TRUE)

.col_order <- c(
  # identity
  "Plate", "drug", "disease_gap_log2",
  # MOVEMENT — BAM effect + confidence
  "DrugEffect_vsKO_est_log2", "DrugEffect_vsKO_SE", "DrugEffect_vsKO_raw_p",
  "DrugEffect_vsKO_adj_p_FDR",
  "bam_move_conf_log2", "bam_move_svalue", "bam_move_lfsr",
  # MOVEMENT — robust Glass's delta (read-only, vs KO)
  "Glass_vsKO_log",
  # LANDING — BAM raw + EXISTING posterior equivalence (read as-is, NOT shrunk)
  "Rescue_vsWT_est_log2", "Rescue_vsWT_SE", "Rescue_vsWT_raw_p",
  "Rescue_vsWT_adj_p_FDR", "equiv_Pr_vsWT", "equiv_M_median_vsWT",
  # landing robust Glass (vs WT) for completeness
  "Glass_vsWT_log"
)
.present <- intersect(.col_order, names(augmented))
.rest    <- setdiff(names(augmented), .present)
setcolorder(augmented, c(.present, .rest))

aug_out <- file.path(OUTPUT_DIR, "augmented_per_plate.csv")
write_csv(augmented, aug_out)
cat("Saved: ", aug_out, "  (", nrow(augmented), " rows)\n", sep = "")

checkA <- rbindlist(checkA_rows, use.names = TRUE, fill = TRUE)
write_csv(checkA, file.path(OUTPUT_DIR, "check_A_glass_vs_bam.csv"))
cat("Saved: check_A_glass_vs_bam.csv  (", nrow(checkA), " rows)\n", sep = "")

# Compact Spearman summary to the console (the headline of Check A).
cat("\n── Check A Spearman summary (per plate) ──\n")
print(unique(checkA[, .(Plate, Glass_vs_BAM = spearman)]), row.names = FALSE)

cat("\nDone.\n")
