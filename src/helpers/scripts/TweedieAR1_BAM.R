# ── PACKAGE CHECK ─────────────────────────────────────────────────────────────
# Packages are managed by the Python app (installed into App/R/library/).
# R_LIBS is set by the app before launching this script, so R finds them there.
cat("R library path:", paste(.libPaths(), collapse = "\n               "), "\n")
required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel", "emmeans", "jsonlite",
                   "gratia", "MASS", "itsadug")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them with:  Rscript App/R/install_packages.R"
  )
}

library(tidyverse)
library(data.table)
library(mgcv)
library(parallel)
library(emmeans)
library(jsonlite)
library(gratia)
# MASS::mvrnorm used directly without library() to avoid masking dplyr::select
# itsadug is loaded conditionally — only the ablation metric block uses it,
# and only acf_resid() / start_value_rho() get called.

# Console width — R wraps printed data.frames whose column widths exceed
# `width`. The default is ~80 chars which splits the master table into two
# chunks. The app's R log pane is much wider than that, so we set width to
# a generous value so master_results prints in one contiguous block.
# 320 is a safe-but-generous default: any terminal/widget narrower than
# that will soft-wrap, which is no worse visually than R's own split would
# have been, and wider surfaces get clean single-line output.
options(width = 320)

# Print warnings immediately (warn=1) rather than queueing them for end-of-
# script flush. The deferred flush has been observed to produce a parse-
# error trailer ("unexpected ')' in ' key)'") on some Rscript / mgcv combos,
# which then surfaces as "Execution halted" AFTER the script's own last
# cat() — confusing both the BAM widget's error dialog (which shows the
# tail of the log) and the user. Immediate warnings appear in the log in
# the order they fire, and the script exits cleanly after the final cat.
options(warn = 1)

# Null-coalescing operator (used when reading optional sidecar fields)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Multiple-testing correction algorithms + tree metadata helpers
.script_dir <- tryCatch({
  args0 <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args0, value = TRUE)
  if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else getwd()
}, error = function(e) getwd())
source(file.path(.script_dir, "corrections.R"))
source(file.path(.script_dir, "tree_metadata.R"))

# ── ARGUMENT PARSING ──────────────────────────────────────────────────────────
# Called from the app as:
#  Rscript "TweedieAR1_BAM.R" <input_csv> <output_dir> <var_names> <ref_values> <ref_condition> <global_corr> <contrast_adj> <roles> <family_name> <shift_val>
# var_names and ref_values are comma-separated strings, e.g. "Genotype,Drug" and "WT,DMSO"
# family_name: "Tweedie" | "Gamma" | "NegBinomial"  (default: "Tweedie")
# shift_val:   numeric constant added to pxl_diff before fitting (default: 0)
args <- commandArgs(trailingOnly = TRUE)
input_csv  <- args[1]
output_dir <- if (length(args) >= 2 && !is.na(args[2])) args[2] else dirname(input_csv)

if (is.na(input_csv) || !file.exists(input_csv)) {
  stop("ERROR: input_csv argument missing or file not found: ", input_csv)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ── DYNAMIC VARIABLE PARSING ─────────────────────────────────────────────────
var_names_raw <- if (length(args) >= 3 && !is.na(args[3]) && nchar(args[3]) > 0) args[3] else "Genotype,Drug"
ref_values_raw <- if (length(args) >= 4 && !is.na(args[4]) && nchar(args[4]) > 0) args[4] else "WT,DMSO"
ref_condition  <- if (length(args) >= 5 && !is.na(args[5]) && nchar(args[5]) > 0) args[5] else "WT+DMSO"
# Legacy positional args — the Python BAM widget still passes "BH"/"none" here
# for back-compat with older R scripts, but the new correction architecture
# drives every step from correction.json (read further down) instead. These
# values aren't used by the analysis any more.
.legacy_global_correction <- if (length(args) >= 6 && !is.na(args[6])) args[6] else "BH"
.legacy_contrast_adjust   <- if (length(args) >= 7 && !is.na(args[7])) args[7] else "none"
roles_raw                 <- if (length(args) >= 8 && !is.na(args[8]) && nchar(args[8]) > 0) args[8] else ""
family_name       <- if (length(args) >= 9 && !is.na(args[9]) && nchar(args[9]) > 0) args[9] else "Tweedie"
shift_val         <- if (length(args) >= 10 && !is.na(args[10]) && nchar(args[10]) > 0) as.numeric(args[10]) else 0

# Build the family object from the name
chosen_family <- switch(family_name,
  Tweedie     = tw(),
  Gamma       = Gamma(link = "log"),
  NegBinomial = nb(),
  {
    warning("Unknown family '", family_name, "' — defaulting to Tweedie.")
    tw()
  }
)

var_names_raw_vec <- str_split(var_names_raw, ",")[[1]]
ref_values        <- str_split(ref_values_raw, ",")[[1]]
n_vars            <- length(var_names_raw_vec)

# Make variable names R-safe (replace spaces, special chars with dots)
# Keep original names for display (plot titles, labels)
var_names_display <- var_names_raw_vec
var_names         <- make.names(var_names_raw_vec, unique = TRUE)

# Build a named reference vector: c(Genotype = "WT", Drug = "DMSO", ...)
ref_map <- setNames(ref_values[1:n_vars], var_names)

# Build display-name map for plot labels
var_display_map <- setNames(var_names_display, var_names)

# Parse roles: "WT+DMSO=RC,KO+DMSO=EXP,..." → named vector c("WT+DMSO"="RC", ...)
role_map <- c()
if (nchar(roles_raw) > 0) {
  role_pairs <- str_split(roles_raw, ",")[[1]]
  for (pair in role_pairs) {
    parts <- str_split(pair, "=")[[1]]
    if (length(parts) == 2) {
      role_map[parts[1]] <- parts[2]
    }
  }
}

# Helper: append role abbreviation to a condition name for plot labels
label_with_role <- function(cond) {
  role <- role_map[cond]
  if (!is.na(role)) paste0(cond, "\n[", role, "]") else cond
}

# Smooth-basis dimension (k) used for every time-series smooth in the formula.
# Defined here (rather than just-in-time) so the run log can echo it back to
# the user alongside the other run parameters.
k_start <- 30

# Random-effect basis for animal_id and plate. Picked via env var so the
# Architecture × Method Ablation tab can flip between options without
# changing the CLI. Values (increasing flexibility):
#   none     — no animal/plate random effect (baseline; AR(1) only)
#   re       — pure random intercept per level (no time argument)
#   re_slope — random intercept + random slope over time_in_group
#   fs       — factor smooth over time_in_group (unconstrained per-animal trajectory)
#   sz       — sum-to-zero constrained smooth over time_in_group (deviation-from-mean)
#
# Default = "sz" — selected as the production winner via the ablation
# diagnostic (see ABLATION_METHODOLOGY.md at the project root: Arch C HGAM
# with sz random basis won by ~12% lower 5-fold CV deviance vs Arch A).
random_basis <- Sys.getenv("RANDOM_BASIS", "sz")
.valid_bases <- c("none", "re", "re_slope", "fs", "sz")
if (!(random_basis %in% .valid_bases)) {
  warning("RANDOM_BASIS='", random_basis,
          "' not in {", paste(.valid_bases, collapse=", "),
          "} — defaulting to 're'.")
  random_basis <- "re"
}

# ── ABLATION ENV VARS ───────────────────────────────────────────────────────
# Set by the Architecture × Method Ablation tab. When ABL_ARCH is non-empty
# the script runs in ablation mode: it follows the architecture-specific
# formula construction, optionally subsets the data (Arch A group), bypasses
# Pass-1 rho estimation when ABL_RHO_FIXED is given, skips Bayesian posterior
# equivalence, and writes a run_metrics.json sidecar at the end. CV mode
# (ABL_CV_FOLD non-empty) further restricts the fit to train animals and
# writes cv_metrics.json after computing holdout deviance.
abl_arch         <- Sys.getenv("ABL_ARCH",          "")  # "A", "B", "C", or ""
abl_by_variable  <- Sys.getenv("ABL_BY_VARIABLE",   "")  # cond_combo|group_factor|interaction
abl_ar_start     <- Sys.getenv("ABL_AR_START_MODE", "")  # animal|phase
abl_group_index  <- Sys.getenv("ABL_GROUP_INDEX",   "")  # Arch A: which Group to fit
abl_rho_fixed    <- Sys.getenv("ABL_RHO_FIXED",     "")  # numeric string or ""
abl_cv_fold      <- Sys.getenv("ABL_CV_FOLD",       "")  # 0..4 or ""
abl_skip_post    <- toupper(Sys.getenv("ABL_SKIP_POSTERIOR", "")) %in% c("1","TRUE","T","YES")
abl_mode         <- nzchar(abl_arch)
if (abl_mode) {
  if (!(abl_arch %in% c("A","B","C"))) {
    stop("ABL_ARCH must be one of {A, B, C} (got '", abl_arch, "').")
  }
  if (abl_arch %in% c("B","C") && !(abl_by_variable %in% c("cond_combo","group_factor","interaction"))) {
    stop("ABL_BY_VARIABLE required for Arch B/C: cond_combo|group_factor|interaction")
  }
  if (!(abl_ar_start %in% c("animal","phase"))) {
    stop("ABL_AR_START_MODE must be 'animal' or 'phase'.")
  }
  cat("\n── Ablation mode ───────────────────────────────────────\n")
  cat("  Architecture:        ", abl_arch, "\n")
  cat("  By-variable:         ", abl_by_variable, "\n")
  cat("  AR.start mode:       ", abl_ar_start, "\n")
  if (nzchar(abl_group_index)) cat("  Group (Arch A):      ", abl_group_index, "\n")
  if (nzchar(abl_rho_fixed))   cat("  Rho fixed:           ", abl_rho_fixed, "\n")
  if (nzchar(abl_cv_fold))     cat("  CV holdout fold:     ", abl_cv_fold, "\n")
  cat("  Skip posterior eq.:  ", abl_skip_post, "\n")
}

cat("Input CSV:        ", input_csv, "\n")
cat("Output dir:       ", output_dir, "\n")
cat("Variables:        ", paste(var_names_display, collapse = ", "), "\n")
# Only show the make.names()-cleaned form if it actually differs from the original
if (!identical(var_names_display, var_names)) {
  cat("Variables (R):    ", paste(var_names, collapse = ", "), "\n")
}
cat("Reference values: ", paste(ref_map, collapse = ", "), "\n")
cat("Ref Condition:    ", ref_condition, "\n")
# emmeans contrasts now always run with adjust = "none"; multiple-testing
# correction (flat BH/Holm or tree TreeBH/Holm-gatekeeping/graphicalMCP) is
# applied post-hoc against the raw p-values via the correction.json sidecar.
cat("Roles:            ", roles_raw, "\n")
cat("Family:           ", family_name, "\n")
cat("Shift applied:    ", shift_val, "\n")
cat("k (smooth basis): ", k_start, "\n")
cat("Random-eff basis: ", random_basis, "\n")

# Find how many cores the computer has, and leave 1 free so the computer doesn't freeze
usable_cores <- max(1, detectCores() - 1)

# ── DATA LOADING ──────────────────────────────────────────────────────────────
# The Python app exports a CSV with columns:
#   time_sec, location, loc_coord, pixel_diff, Condition, Phase, Group, animal_id, plate
#   + individual variable columns (Genotype, Drug, etc.)
full_df <- read_csv(input_csv, show_col_types = FALSE) %>%
  rename(pxl_diff = pixel_diff) %>%
  mutate(
    Condition = as.character(Condition),
    Phase     = as.character(Phase),
    Group     = as.integer(Group),
    time_sec  = as.integer(time_sec),
    pxl_diff  = as.numeric(pxl_diff),
    plate     = as.character(plate),
    location  = as.character(location),
    loc_coord = as.character(loc_coord)
  ) %>%
  filter(!is.na(Condition), Condition != "")

# If individual variable columns are already in the CSV, use them;
# otherwise split Condition on "+"
# Check using the original (display) names since CSV columns come from Python
if (!all(var_names_display %in% colnames(full_df))) {
  full_df <- full_df %>%
    separate(Condition, into = var_names, sep = "\\+", remove = FALSE)
} else {
  # Rename columns from display names to R-safe names
  for (i in seq_along(var_names)) {
    if (var_names_display[i] != var_names[i] && var_names_display[i] %in% colnames(full_df)) {
      full_df <- full_df %>% rename(!!var_names[i] := !!var_names_display[i])
    }
  }
}

# ---------------------------------------------------------
# 1. FINAL DATA PREP FOR MGCV
# ---------------------------------------------------------
# Unified-model design:
#   - One bam() fit across all phase groups, not one per group.
#   - time_in_group resets at each phase boundary so per-(condition × phase)
#     smooths live entirely inside their own phase.
#   - time_sec is the global continuous-time variable used by the per-well
#     and per-plate random-effect smooths so they span the whole recording.
#   - Group_factor encodes phase identity as a categorical predictor for the
#     parametric expansion and the by-interaction smooth.
#   - start_event resets the AR(1) chain at BOTH well starts AND phase
#     boundaries within wells — otherwise mgcv would attempt to correlate
#     residuals across a discontinuous phase transition.
# AR.start mode controls where the AR(1) chain restarts. Default ("phase")
# resets at both animal-start and phase-boundary within an animal. The
# ablation tab can switch to "animal" via ABL_AR_START_MODE — chains restart
# only at animal boundaries, treating each animal as one continuous time
# series across phases.
.ar_phase_resets <- if (abl_mode) (abl_ar_start == "phase") else TRUE
gam_df <- full_df %>%
  mutate(
    animal_id       = as.factor(paste(plate, location, sep = "_")),
    plate           = as.factor(plate),
    Condition_Combo = as.factor(Condition),
    Group_factor    = as.factor(Group)
  ) %>%
  group_by(Group) %>%
  mutate(time_in_group = time_sec - min(time_sec)) %>%
  ungroup() %>%
  arrange(animal_id, Group, time_in_group) %>%
  group_by(animal_id) %>%
  mutate(start_event = if (.ar_phase_resets) {
                         (row_number() == 1) | (Group != lag(Group, default = -1L))
                       } else {
                         row_number() == 1
                       }) %>%
  ungroup()

# ── Arch A: subset to a single phase Group ─────────────────────────────────
if (abl_mode && abl_arch == "A") {
  if (!nzchar(abl_group_index)) {
    stop("ABL_GROUP_INDEX required when ABL_ARCH=A.")
  }
  g_id <- as.integer(abl_group_index)
  gam_df <- gam_df[gam_df$Group == g_id, , drop = FALSE]
  if (nrow(gam_df) == 0) stop("No data rows for Group ", g_id, " in Arch A subset.")
  # Drop now-empty Group_factor levels so Arch A's per-phase model doesn't
  # carry phantom parametric columns from other groups.
  gam_df$Group_factor <- droplevels(gam_df$Group_factor)
  # Re-baseline time_in_group (it's already per-Group, but recompute for safety)
  gam_df$time_in_group <- as.integer(gam_df$time_sec - min(gam_df$time_sec))
  cat("Arch A: subset to Group ", g_id, " (", nrow(gam_df), " rows).\n", sep = "")
}

# ── CV mode: split animals into 5 folds, drop the holdout from training ────
.cv_holdout_df <- NULL
if (abl_mode && nzchar(abl_cv_fold)) {
  fold_id <- as.integer(abl_cv_fold)
  if (!(fold_id %in% 0:4)) stop("ABL_CV_FOLD must be 0..4 (got ", abl_cv_fold, ").")
  # Deterministic seeded animal-to-fold assignment so all variants see the
  # same split. Sort first to make the order deterministic across machines.
  uniq_animals <- sort(levels(gam_df$animal_id))
  set.seed(42L)
  fold_assign  <- sample(rep(0:4, length.out = length(uniq_animals)))
  names(fold_assign) <- uniq_animals
  holdout_animals <- uniq_animals[fold_assign == fold_id]
  cat("CV fold ", fold_id, ": holding out ", length(holdout_animals),
      " animal(s) for prediction.\n", sep = "")
  .cv_holdout_df <- gam_df[gam_df$animal_id %in% holdout_animals, , drop = FALSE]
  gam_df         <- gam_df[!(gam_df$animal_id %in% holdout_animals), , drop = FALSE]
  if (nrow(gam_df) == 0) stop("CV training set is empty.")
}

# Convert each design variable to factor and set reference level globally
for (v in var_names) {
  gam_df[[v]] <- as.factor(gam_df[[v]])
  gam_df[[v]] <- relevel(gam_df[[v]], ref = ref_map[[v]])
}
gam_df$Condition_Combo <- relevel(gam_df$Condition_Combo, ref = ref_condition)

unique_groups <- sort(unique(gam_df$Group))
n_groups      <- length(unique_groups)
cat("Phase groups detected:", n_groups, "  →  ", paste(unique_groups, collapse = ", "), "\n")

# Apply shift to pxl_diff if requested (required when using Gamma family)
if (shift_val != 0) {
  gam_df <- gam_df %>% mutate(pxl_diff = pxl_diff + shift_val)
  cat(sprintf("\nShift of %.4f applied to pxl_diff (family: %s).\n", shift_val, family_name))
}

# ---------------------------------------------------------
# 2. BUILD UNIFIED BAM FORMULA
# ---------------------------------------------------------
# Parametric:   Var1 * Var2 * ... * VarN * Group_factor
#               → all design-variable main effects + all interactions
#                 INCLUDING phase-group interactions; captures per-(condition,
#                 phase) intercept shifts and any cross-phase modulation.
# Smooths:
#   s(time_in_group, k = k_start)                                  — shared within-phase shape baseline
#   s(time_in_group, by = interaction(Condition_Combo, Group_factor), k = k_start)
#                                                                  — per-(condition × phase) trajectory deviation
#   s(time_sec, plate,     bs = 'fs', m = 1)                       — plate random smooth across all time
#   s(time_sec, animal_id, bs = 'fs', m = 1)                       — well random smooth across all time
#
# Transition spikes between phases are inferred from the early-time shape of
# the next phase's per-(condition, phase) smooth plus the parametric phase
# intercepts. The random-effect smooths over time_sec span the entire
# recording so per-well / per-plate variance is estimated once jointly.
# Random-effect smooth string builder. Returns the formula fragment for
# ONE factor (animal_id or plate) given the chosen basis and time variable.
# Returns "" for "none" so the caller can omit the term entirely.
#   re       → intercept only (no time argument)
#   re_slope → intercept + random slope on `time_var`
#   fs       → unconstrained factor smooth on `time_var`
#   sz       → sum-to-zero constrained smooth on `time_var`
build_re_smooth <- function(factor_var, basis, time_var = "time_in_group") {
  switch(basis,
    none     = "",
    re       = sprintf("s(%s, bs='re')", factor_var),
    re_slope = sprintf("s(%s, bs='re') + s(%s, %s, bs='re')",
                       factor_var, time_var, factor_var),
    fs       = sprintf("s(%s, %s, bs='fs', m=1)", time_var, factor_var),
    sz       = sprintf("s(%s, %s, bs='sz')",      time_var, factor_var)
  )
}

# Prepend " + " only when the term is non-empty so a "none" choice doesn't
# leak a dangling plus into the formula.
.prefix_re <- function(term) if (nzchar(term)) paste0(" + ", term) else ""

# ── Pick the time variable + by-variable for the chosen architecture ───────
#   Arch A — per-phase: fixed smooths on time_in_group, by = Condition_Combo
#            (Group_factor is single-valued in the subset, omitted from parametric)
#   Arch B — fully global: smooths on time_sec, by = user-selected (cond_combo /
#            group_factor / interaction); parametric includes Group_factor
#   Arch C — HGAM: smooths on time_in_group, by = user-selected; parametric
#            includes Group_factor (current architecture)
re_time_var <- "time_in_group"
fixed_time_var <- "time_in_group"
if (abl_mode && abl_arch == "B") {
  re_time_var    <- "time_sec"
  fixed_time_var <- "time_sec"
}

build_by_smooth <- function(by_choice, time_var) {
  by_expr <- switch(by_choice,
    cond_combo    = "Condition_Combo",
    group_factor  = "Group_factor",
    interaction   = "interaction(Condition_Combo, Group_factor)"
  )
  sprintf("s(%s, by = %s, k = k_start)", time_var, by_expr)
}

n_plates   <- length(unique(gam_df$plate))
plate_term <- build_re_smooth("plate", random_basis, re_time_var)
if (n_plates < 2 && nzchar(plate_term)) {
  cat("Only ", n_plates, " plate(s) detected — omitting the plate factor smooth.\n",
      sep = "")
  plate_term <- ""
}
plate_smooth_str  <- .prefix_re(plate_term)
animal_smooth_str <- .prefix_re(build_re_smooth("animal_id", random_basis, re_time_var))

# Parametric expansion — Arch A drops Group_factor because the subset has
# a single Group level. B and C include it (interaction term spans phases).
if (abl_mode && abl_arch == "A") {
  parametric_str <- paste(var_names, collapse = " * ")
} else {
  parametric_str <- paste(c(var_names, "Group_factor"), collapse = " * ")
}

# Trajectory smooth: global shape + by-smooth deviations.
# Default by-variable matches the legacy unified-model choice when not in
# ablation mode (interaction across Condition_Combo and Group_factor).
default_by_choice <- "interaction"
if (abl_mode) {
  if (abl_arch == "A") {
    # Arch A: only Condition_Combo varies within a single phase
    chosen_by_choice <- "cond_combo"
  } else {
    chosen_by_choice <- abl_by_variable
  }
} else {
  chosen_by_choice <- default_by_choice
}
by_smooth_term <- build_by_smooth(chosen_by_choice, fixed_time_var)

formula_str <- paste0(
  "pxl_diff ~ ", parametric_str,
  " + s(", fixed_time_var, ", k = k_start)",
  " + ", by_smooth_term,
  plate_smooth_str,
  animal_smooth_str
)
if (abl_mode) {
  cat("\nBAM formula (Arch ", abl_arch, ", by=", chosen_by_choice,
      ", re=", random_basis, ", AR.start=", abl_ar_start, "):\n  ",
      formula_str, "\n\n", sep = "")
} else {
  cat("\nBAM formula (unified across all phase groups):\n  ", formula_str, "\n\n")
}

bam_formula <- as.formula(formula_str)

# ---------------------------------------------------------
# B. PASS 1 — ESTIMATE AR(1) RHO FROM RESIDUALS
# ---------------------------------------------------------
# Procedure: standard itsadug start_value_rho equivalent — fit a no-AR model,
# then take lag-1 autocorrelation of the overall residual vector. Returns a
# single ρ that Pass 2 plugs in. (itsadug start_value_rho() body is
# acf(resid(model), plot=FALSE)$acf[2]; we inline the call to avoid the extra
# package dependency at this site — itsadug is still loaded for diagnostics.)
# Skipped entirely when ABL_RHO_FIXED is set — useful for ablation runs where
# we want ρ identical across variants.
fallback_rho <- 0.20
if (abl_mode && nzchar(abl_rho_fixed)) {
  optimal_rho <- as.numeric(abl_rho_fixed)
  optimal_rho <- max(min(optimal_rho, 0.95), -0.95)
  cat(sprintf("\nUsing fixed rho = %.4f (ablation mode — Pass 1 skipped).\n",
              optimal_rho))
} else {
  cat("\n========================================\n")
  cat("Pass 1: fitting no-AR model to estimate autocorrelation (rho)...\n")
  optimal_rho <- tryCatch({
    model_no_ar <- bam(
      formula  = bam_formula,
      data     = gam_df,
      family   = chosen_family,
      select   = TRUE,
      method   = "fREML",
      discrete = TRUE,
      nthreads = usable_cores
    )
    rho_est <- acf(resid(model_no_ar), plot = FALSE)$acf[2]
    rho_est <- max(min(rho_est, 0.95), -0.95)
    rm(model_no_ar)
    rho_est
  }, error = function(e) {
    warning("Pass 1 failed (", conditionMessage(e),
            ") — falling back to fixed rho = ", fallback_rho)
    fallback_rho
  })
  cat(sprintf("Pass 1 estimated rho: %.4f\n", optimal_rho))
}

# ---------------------------------------------------------
# C. PASS 2 — FIT THE MODEL WITH AR(1) CORRECTION
# ---------------------------------------------------------
cat(sprintf("\nPass 2: fitting unified model with AR1 correction (rho = %.4f)...\n",
            optimal_rho))

final_model <- bam(
  formula  = bam_formula,
  data     = gam_df,
  family   = chosen_family,
  rho      = optimal_rho,
  AR.start = gam_df$start_event,
  select   = TRUE,
  method   = "fREML",
  discrete = TRUE,
  nthreads = usable_cores
)

cat("Unified model fitted across all", n_groups, "phase groups.\n")

# ---------------------------------------------------------
# CV MODE — predict on holdout and exit before contrast work
# ---------------------------------------------------------
# When ABL_CV_FOLD is set, the script runs in CV mode: fit on training
# animals only, predict on the held-out animals at the population level
# (random-effect smooths excluded), and write the holdout Tweedie deviance.
# The contrast/correction/posterior blocks below are skipped — they would
# only waste compute since the widget cares about deviance only.
if (abl_mode && !is.null(.cv_holdout_df)) {
  cat("\n========================================\n")
  cat("CV mode: predicting on holdout animals (population-level)\n")

  # Random-effect smooth term labels to exclude from the linear predictor —
  # the held-out animals don't have estimated random coefficients, so we
  # predict their trajectory at the population mean.
  excl_terms <- switch(random_basis,
    none     = character(0),
    re       = c("s(animal_id)", "s(plate)"),
    re_slope = c("s(animal_id)", "s(plate)",
                 sprintf("s(%s,animal_id)", re_time_var),
                 sprintf("s(%s,plate)",     re_time_var)),
    fs       = c(sprintf("s(%s,animal_id)", re_time_var),
                 sprintf("s(%s,plate)",     re_time_var)),
    sz       = c(sprintf("s(%s,animal_id)", re_time_var),
                 sprintf("s(%s,plate)",     re_time_var))
  )

  # animal_id and plate factor levels in the holdout aren't in the training
  # set; mgcv predict() needs levels it knows. Since we exclude the random
  # smooths above, the actual level value doesn't enter the linear predictor
  # — assign a training-set placeholder so predict() doesn't error.
  #
  # IMPORTANT: gam_df and .cv_holdout_df come from a tidyverse pipeline so
  # they are tibbles. mgcv's predict() can fail with "cannot xtfrm data
  # frames" when any newdata column has tibble-like structure that tripsorder()/match() during the model-frame build. Coerce to plain data.frame
  # first to strip that.
  cv_pred_df <- as.data.frame(.cv_holdout_df)
  # Use the FIRST training-set animal/plate as the placeholder. Both factor
  # columns get re-factored against the training-set levels (droplevels'd to
  # contain only those actually present) so mgcv has unambiguous coefs to
  # ignore via `exclude` below.
  train_animals <- droplevels(factor(as.character(as.data.frame(gam_df)$animal_id)))
  train_plates  <- droplevels(factor(as.character(as.data.frame(gam_df)$plate)))
  placeholder_animal <- levels(train_animals)[1]
  placeholder_plate  <- levels(train_plates)[1]
  cv_pred_df$animal_id <- factor(placeholder_animal, levels = levels(train_animals))
  cv_pred_df$plate     <- factor(placeholder_plate,  levels = levels(train_plates))

  # Predict in two stages so we can pinpoint whether prediction or deviance
  # computation is the failure mode if the holdout deviance comes out NA.
  cv_predict_error  <- NA_character_
  cv_devresid_error <- NA_character_
  cv_n_invalid_pred <- NA_integer_
  cv_deviance       <- NA_real_

  preds <- tryCatch({
    predict(final_model, newdata = cv_pred_df,
            type = "response", exclude = excl_terms)
  }, error = function(e) {
    cv_predict_error <<- conditionMessage(e)
    warning("CV predict() failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(preds)) {
    n_bad <- sum(is.na(preds) | !is.finite(preds))
    cv_n_invalid_pred <- as.integer(n_bad)
    cat(sprintf("CV preds: %d total, %d invalid (NA/Inf)\n",
                length(preds), n_bad))
    if (n_bad == length(preds)) {
      cv_predict_error <- "all predictions invalid (NA/Inf)"
    } else {
      cv_deviance <- tryCatch({
        fam <- final_model$family
        obs <- .cv_holdout_df$pxl_diff
        wts <- rep(1, length(obs))
        if (is.null(fam$dev.resids)) {
          cv_devresid_error <<- "family$dev.resids is NULL"
          NA_real_
        } else {
          dev_per_obs <- fam$dev.resids(obs, preds, wts)
          n_bad_dev <- sum(is.na(dev_per_obs) | !is.finite(dev_per_obs))
          if (n_bad_dev > 0) {
            cat(sprintf("dev.resids: %d/%d invalid (NA/Inf)\n",
                        n_bad_dev, length(dev_per_obs)))
          }
          total <- sum(dev_per_obs, na.rm = TRUE)
          if (is.na(total) || !is.finite(total)) {
            cv_devresid_error <<- "sum(dev_per_obs) is NA/Inf"
            NA_real_
          } else {
            total
          }
        }
      }, error = function(e) {
        cv_devresid_error <<- conditionMessage(e)
        warning("CV dev.resids failed: ", conditionMessage(e))
        NA_real_
      })
    }
  }

  cv_metrics <- list(
    fold              = as.integer(abl_cv_fold),
    arch              = abl_arch,
    random_basis      = random_basis,
    by_variable       = chosen_by_choice,
    ar_start_mode     = abl_ar_start,
    group_index       = if (nzchar(abl_group_index)) as.integer(abl_group_index) else NA_integer_,
    n_train_obs       = nrow(gam_df),
    n_holdout_obs     = nrow(.cv_holdout_df),
    n_holdout_animals = length(unique(as.character(.cv_holdout_df$animal_id))),
    holdout_deviance  = cv_deviance,
    n_invalid_pred    = cv_n_invalid_pred,
    predict_error     = cv_predict_error,
    devresid_error    = cv_devresid_error,
    excl_terms_used   = excl_terms
  )
  jsonlite::write_json(cv_metrics, file.path(output_dir, "cv_metrics.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
  cat(sprintf("CV holdout deviance: %s\n",
              if (is.na(cv_deviance)) "NA" else sprintf("%.4f", cv_deviance)))
  cat("Saved: cv_metrics.json. Exiting.\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------
# 3. GLOBAL MULTIPLE COMPARISONS: ACROSS ALL PHASE GROUPS
# ---------------------------------------------------------
# Dynamically generate per-variable contrasts and the full interaction contrast

# ── Optional contrast-selection sidecar ────────────────────────────────────
# The Python app may write contrast_selection.json into the output dir to
# narrow the Full_Interaction family to a user-chosen subset of pairs. The
# filter is applied BEFORE building the emmeans contrast set so the local
# (Sidak / Dunnett / Tukey) and global (FDR / Holm) corrections both operate
# on the trimmed family.
#   - file present, kept_pairs non-empty → build only those pairs
#   - file present, kept_pairs empty     → skip Full_Interaction entirely
#   - file absent                        → fall back to method = "pairwise"
contrast_sidecar <- file.path(output_dir, "contrast_selection.json")
if (file.exists(contrast_sidecar)) {
  contrast_spec <- tryCatch(
    fromJSON(contrast_sidecar, simplifyVector = FALSE),
    error = function(e) {
      warning("Could not parse contrast_selection.json: ", conditionMessage(e),
              " — falling back to full pairwise.")
      NULL
    }
  )
  if (!is.null(contrast_spec)) {
    cat("Contrast selection sidecar found: ",
        length(contrast_spec$kept_pairs), " pair(s) requested.\n", sep = "")
  }
} else {
  contrast_spec <- NULL
}

all_contrasts_list <- list()

# ── Reference-grid sizing ───────────────────────────────────────────────────
# Single unified model: the grid now expands across (var levels × Group_factor).
# Bump rg.limit to accommodate.
random_factor_names <- character(0)
for (sm in final_model$smooth) {
  if (inherits(sm, "fs.interaction")) {
    random_factor_names <- c(random_factor_names, sm$fterm)
  } else if (inherits(sm, "random.effect")) {
    random_factor_names <- c(random_factor_names, sm$term)
  }
}
random_factor_names <- unique(random_factor_names)

candidate_factor_vars <- intersect(
  c(var_names, "Condition_Combo", "Group_factor", "animal_id"),
  colnames(gam_df)
)
factor_sizes <- vapply(candidate_factor_vars, function(v) {
  if (is.factor(gam_df[[v]])) nlevels(gam_df[[v]]) else 1L
}, integer(1))
fixed_factor_vars  <- setdiff(candidate_factor_vars, random_factor_names)
fixed_factor_sizes <- factor_sizes[fixed_factor_vars]
grid_rows <- if (length(fixed_factor_sizes) > 0)
               prod(as.numeric(fixed_factor_sizes)) else 1
cat(sprintf("\nRef-grid fixed factors: %s\n",
            paste(sprintf("%s(%d)", fixed_factor_vars, fixed_factor_sizes),
                  collapse = " × ")))
cat(sprintf("Estimated emmeans grid: %d rows\n", as.integer(grid_rows)))
emm_options(rg.limit = max(50000, as.integer(grid_rows * 2)))

# ── Per-variable effects ────────────────────────────────────────────────────
# Group_factor joins the conditioning so the same trt-vs-ref contrast is
# produced once per (other-vars × Group) cell from the single unified model.
# In Arch A the dataset is subset to a single Group, so Group_factor has
# only one level and emmeans strips it from the reference grid. Skip the
# Group_factor conditioning in that case and fill Group from abl_group_index.
.use_group_in_emm <- !(abl_mode && abl_arch == "A")
.arch_a_group_id  <- if (abl_mode && abl_arch == "A") as.integer(abl_group_index) else NA_integer_

for (vi in seq_along(var_names)) {
  focal_var   <- var_names[vi]
  other_vars  <- var_names[-vi]
  ref_val     <- ref_map[[focal_var]]

  cond_vars <- if (.use_group_in_emm) c(other_vars, "Group_factor") else other_vars
  if (length(cond_vars) > 0) {
    spec_str <- paste0("~ ", focal_var, " | ", paste(cond_vars, collapse = " + "))
  } else {
    spec_str <- paste0("~ ", focal_var)
  }
  spec_formula <- as.formula(spec_str)

  emm_base <- emmeans(final_model, specs = spec_formula)
  emm_contrasts <- contrast(emm_base, method = "trt.vs.ctrl",
                            ref = ref_val, adjust = "none")

  df_contrasts <- as.data.frame(emm_contrasts)
  if (.use_group_in_emm) {
    df_contrasts <- df_contrasts %>%
      mutate(
        Test_Family  = paste0(focal_var, "_Effect"),
        Passed_Contrasts = as.character(contrast),
        Group        = as.integer(as.character(Group_factor))
      ) %>%
      select(-Group_factor)
  } else {
    df_contrasts <- df_contrasts %>%
      mutate(
        Test_Family  = paste0(focal_var, "_Effect"),
        Passed_Contrasts = as.character(contrast),
        Group        = .arch_a_group_id
      )
  }

  # Collapse remaining other-variable columns into a single Split_By string
  split_cols <- intersect(other_vars, colnames(df_contrasts))
  if (length(split_cols) > 0) {
    df_contrasts$Split_By <- apply(df_contrasts[, split_cols, drop = FALSE], 1,
                                    function(x) paste(x, collapse = " + "))
    df_contrasts <- df_contrasts %>% select(-all_of(split_cols))
  } else {
    df_contrasts$Split_By <- "All"
  }

  all_contrasts_list[[focal_var]] <- df_contrasts
}

# ── Full interaction: pairwise condition comparisons by Group ───────────────
# Use the DESIGN variables (Var1 * Var2 * ...) — NOT Condition_Combo — because
# emmeans detects Condition_Combo's levels as nested in the parametric
# Genotype*Drug*Group_factor expansion and strips it out of the reference
# grid. Asking for Condition_Combo by name then errors with
#   "No variable named Condition_Combo in the reference grid".
# Building the spec from the design variables sidesteps the nesting issue
# entirely; we reconstruct the "+"-joined condition string from the grid
# columns for matching against the sidecar's kept_pairs entries.
#
# In Arch A the dataset is subset to one Group so Group_factor is single-
# valued — drop it from the spec and the by-variable; the contrast block
# fills Group from abl_group_index.
inter_spec <- if (.use_group_in_emm) {
  as.formula(paste0("~ ", paste(var_names, collapse = " * "), " | Group_factor"))
} else {
  as.formula(paste0("~ ", paste(var_names, collapse = " * ")))
}
emm_inter_base <- emmeans(final_model, specs = inter_spec)

# Pull a single by-level's sub-grid to learn the within-group row ordering;
# the weight vectors built below have length = nrow(that sub-grid) and are
# replayed by emmeans across every Group_factor level via by="Group_factor".
full_inter_grid <- emm_inter_base@grid
if (.use_group_in_emm) {
  first_group <- levels(gam_df$Group_factor)[1]
  within_grid <- full_inter_grid[
    as.character(full_inter_grid$Group_factor) == first_group,
    var_names, drop = FALSE
  ]
} else {
  within_grid <- full_inter_grid[, var_names, drop = FALSE]
}
cond_strs <- apply(within_grid, 1,
                   function(r) paste(as.character(r), collapse = "+"))

emm_inter_contrasts <- NULL
if (!is.null(contrast_spec)) {
  kept_pairs <- contrast_spec$kept_pairs
  if (length(kept_pairs) == 0) {
    cat("Full_Interaction skipped (no pairs selected).\n")
  } else {
    custom_contrasts <- list()
    missing <- c()
    for (pair in kept_pairs) {
      lhs <- pair[[1]]; rhs <- pair[[2]]
      lhs_i <- which(cond_strs == lhs)
      rhs_i <- which(cond_strs == rhs)
      if (length(lhs_i) == 1L && length(rhs_i) == 1L) {
        v <- numeric(length(cond_strs))
        v[lhs_i] <-  1
        v[rhs_i] <- -1
        custom_contrasts[[paste(lhs, "-", rhs)]] <- v
      } else {
        missing <- c(missing, paste(lhs, "vs", rhs))
      }
    }
    if (length(missing) > 0) {
      warning("Couldn't locate pair(s) in the design-variable grid: ",
              paste(missing, collapse = "; "))
    }
    if (length(custom_contrasts) > 0) {
      emm_inter_contrasts <- if (.use_group_in_emm) {
        contrast(emm_inter_base, method = custom_contrasts,
                 by = "Group_factor", adjust = "none")
      } else {
        contrast(emm_inter_base, method = custom_contrasts, adjust = "none")
      }
    }
  }
} else {
  # Legacy fallback (no contrast_selection.json sidecar): all-pairs.
  emm_inter_contrasts <- if (.use_group_in_emm) {
    contrast(emm_inter_base, method = "pairwise",
             by = "Group_factor", adjust = "none")
  } else {
    contrast(emm_inter_base, method = "pairwise", adjust = "none")
  }
}

if (!is.null(emm_inter_contrasts)) {
  df_inter <- as.data.frame(emm_inter_contrasts)
  if (.use_group_in_emm) {
    df_inter <- df_inter %>%
      mutate(
        Test_Family  = "Full_Interaction",
        Passed_Contrasts = as.character(contrast),
        Split_By     = "None",
        Group        = as.integer(as.character(Group_factor))
      ) %>%
      select(-Group_factor)
  } else {
    df_inter <- df_inter %>%
      mutate(
        Test_Family  = "Full_Interaction",
        Passed_Contrasts = as.character(contrast),
        Split_By     = "None",
        Group        = .arch_a_group_id
      )
  }
  all_contrasts_list[["interaction"]] <- df_inter
}

# ---------------------------------------------------------
# COMBINE + APPLY MULTIPLE-TESTING CORRECTION
# ---------------------------------------------------------
# The contrast loop above produced RAW p-values (emmeans `adjust = "none"`).
# Now apply the user's chosen correction method to those raw values per the
# correction.json sidecar:
#   - strategy   = "flat" | "tree"
#   - error_rate = "FDR"  | "FWER"
#   - contrast_set = "all_pairs" | "ref_only"   (affects graphicalMCP vs Holm)
#   - threshold  = q (FDR) or alpha (FWER), default 0.05
#   - tree       = list of level specs (only used when strategy=="tree")
#
# If no sidecar is present, default to flat BH per Test_Family at q = 0.05
# (matches the historical behaviour from before the new architecture).
master_results_df <- bind_rows(all_contrasts_list)

final_master_table <- master_results_df %>%
  mutate(
    raw_pvalue = p.value,
    estimate   = estimate / log(2),
    SE         = SE       / log(2),
    Group      = as.integer(Group)
  ) %>%
  select(Test_Family, Group, Split_By, Passed_Contrasts, estimate, SE, raw_pvalue) %>%
  arrange(Test_Family, Group)

# Read correction.json if present
correction_sidecar <- file.path(output_dir, "correction.json")
correction_spec <- if (file.exists(correction_sidecar)) {
  tryCatch(
    fromJSON(correction_sidecar, simplifyVector = FALSE),
    error = function(e) {
      warning("Could not parse correction.json: ", conditionMessage(e),
              " — falling back to flat BH at q=0.05.")
      NULL
    }
  )
} else NULL

if (is.null(correction_spec)) {
  # Default: flat BH per-family at q=0.05
  correction_spec <- list(
    strategy   = "flat",
    error_rate = "FDR",
    threshold  = 0.05,
    scope      = "per_family"
  )
  cat("\nNo correction.json sidecar — using default: flat BH per family at q=0.05.\n")
} else {
  cat("\nCorrection.json found: strategy=", correction_spec$strategy %||% "flat",
      " error_rate=", correction_spec$error_rate %||% "FDR",
      " threshold=", correction_spec$threshold %||% 0.05, "\n", sep = "")
}

# If a tree strategy, annotate rows with tree metadata BEFORE correction.
if (identical(correction_spec$strategy, "tree")) {
  tree_spec <- correction_spec$tree
  if (is.null(tree_spec) || length(tree_spec) == 0) {
    cat("Tree strategy chosen but no tree spec provided — using default rescue tree.\n")
    tree_spec <- default_rescue_tree_spec(var_names, ref_map)
  }
  # Pass var_names and the Python-supplied condition_lookup so design-variable
  # linkage decodes Full_Interaction's multi-variable contrast strings via
  # dict access (no R-side string parsing of "+"-joined condition names).
  final_master_table <- assign_tree_metadata(
    final_master_table, tree_spec, var_names,
    condition_lookup = correction_spec$condition_lookup
  )
}

final_master_table <- apply_correction(final_master_table, correction_spec)

# Tidy up column order — tree columns may not have been added by flat methods
# but they're always present after apply_correction (NAs are fine).
final_master_table <- final_master_table %>%
  select(Test_Family, Group, Split_By, Passed_Contrasts,
         estimate, SE,
         raw_pvalue, adjusted_pvalue, correction_method,
         tidyselect::any_of(c("tree_id", "tree_level",
                              "tree_parent_id", "tree_status")))

cat("\n--- CLEAN MASTER TABLE ---\n")
# Print via as.data.frame to avoid a partial-arg-matching crash when the
# table is a plain data.frame (mgcv select() loses the tibble class).
print(as.data.frame(final_master_table), row.names = FALSE)

# Persist the per-contrast results for downstream visualisation.
write_csv(final_master_table, "master_results.csv")
cat("Saved: master_results.csv\n")

# ---------------------------------------------------------
# SUMMARY TABLE FOR QUICK REFERENCE
# ---------------------------------------------------------

summary_table <- final_master_table %>%
  group_by(Test_Family, Group) %>%
  summarise(
    n_total = n(),
    n_sig = sum(adjusted_pvalue < 0.05, na.rm = TRUE),
    mean_effect = mean(abs(estimate)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(n_sig, mean_effect),
    names_glue = "G{Group}_{.value}"
  )

write_csv(summary_table, "summary_statistics_by_group.csv")

# =============================================================================
# ABLATION METRICS — run_metrics.json
# =============================================================================
# Computed only when ABL_ARCH is set. The widget reads these per variant to
# populate the metrics matrix, the Pareto front, and the random-eats-fixed
# bar chart. Everything is wrapped in tryCatch so an individual diagnostic
# failure doesn't kill the whole run — failed metrics come back as NA.
if (abl_mode) {
  cat("\n========================================\n")
  cat("Ablation metrics: computing run_metrics.json\n")

  safe <- function(expr) tryCatch(eval(expr), error = function(e) NA_real_)

  # Captured at fit time — collect any mgcv warnings that surfaced for this
  # model. warnings() reads the global warning buffer.
  conv_warnings <- tryCatch({
    w <- names(warnings())
    if (is.null(w)) character(0) else as.character(w)
  }, error = function(e) character(0))

  smry <- tryCatch(summary(final_model), error = function(e) NULL)

  # ── Fit-quality block ────────────────────────────────────────────────────
  fit_quality <- list(
    AIC                = safe(quote(AIC(final_model))),
    BIC                = safe(quote(BIC(final_model))),
    logLik             = safe(quote(as.numeric(logLik(final_model)))),
    fREML              = safe(quote(as.numeric(final_model$gcv.ubre))),
    dev_explained      = if (!is.null(smry)) smry$dev.expl %||% NA_real_ else NA_real_,
    r_sq_adj           = if (!is.null(smry)) smry$r.sq    %||% NA_real_ else NA_real_,
    tweedie_phi        = if (!is.null(smry)) smry$scale   %||% NA_real_ else NA_real_,
    total_edf          = safe(quote(sum(smry$edf %||% NA_real_, na.rm = TRUE)))
  )

  # ── Adequacy block ───────────────────────────────────────────────────────
  # k.check returns the same k-index table that gam.check prints.
  k_check_tbl <- tryCatch(mgcv::k.check(final_model), error = function(e) NULL)
  if (!is.null(k_check_tbl)) {
    k_index_p_vals <- k_check_tbl[, "p-value"]
    smooth_pass    <- !is.na(k_index_p_vals) & k_index_p_vals >= 0.05
    k_index_summary <- list(
      n_smooths         = nrow(k_check_tbl),
      n_pass            = sum(smooth_pass, na.rm = TRUE),
      n_fail            = sum(!smooth_pass & !is.na(k_index_p_vals)),
      pass_rate         = if (length(smooth_pass)) mean(smooth_pass) else NA_real_,
      per_smooth_p      = setNames(as.list(k_index_p_vals), rownames(k_check_tbl)),
      per_smooth_edf    = setNames(as.list(k_check_tbl[, "edf"]), rownames(k_check_tbl)),
      per_smooth_k      = setNames(as.list(k_check_tbl[, "k'"]), rownames(k_check_tbl))
    )
  } else {
    k_index_summary <- list(error = "k.check failed")
  }

  # Lag-1 ACF of AR-corrected residuals — should be ~0 if AR(1) captured
  # within-animal autocorrelation. final_model$std.rsd holds the AR-corrected
  # standardized residuals when rho != 0.
  ar1_residual_lag1 <- tryCatch({
    rsd <- if (!is.null(final_model$std.rsd)) final_model$std.rsd
           else resid(final_model)
    acf(rsd, plot = FALSE)$acf[2]
  }, error = function(e) NA_real_)

  adequacy <- list(
    convergence_warnings = conv_warnings,
    k_index              = k_index_summary,
    ar1_residual_lag1    = ar1_residual_lag1
  )

  # ── Identifiability + Random-eats-fixed block ────────────────────────────
  # concurvity(model, full=FALSE)$worst is a matrix; we extract submatrices
  # to separate "random-vs-fixed" from overall maximum concurvity.
  conc <- tryCatch(mgcv::concurvity(final_model, full = FALSE),
                   error = function(e) NULL)
  smooth_labels <- if (!is.null(conc)) rownames(conc$worst) else character(0)

  # Identify random-effect smooth labels directly from the model object
  # rather than by regex on the label text. The label strings differ
  # between summary(), concurvity(), and gam.vcomp() — gam.vcomp in
  # particular appends a numeric suffix to fs.interaction smooths (one
  # entry per smoothing parameter), so a strict label regex misses them.
  .random_factor_set <- c("animal_id", "plate")
  .random_smooth_labels <- character(0)
  for (sm in final_model$smooth) {
    cls <- class(sm)
    # Smooth classes that imply a random-effect-like structure
    is_re_class <- any(c("random.effect", "fs.interaction") %in% cls) ||
                   any(grepl("^sz", cls, ignore.case = TRUE))
    # Grouping variables — random.effect uses $term, factor-smooths use $fterm
    grouping_vars <- character(0)
    if (!is.null(sm$term))  grouping_vars <- c(grouping_vars, sm$term)
    if (!is.null(sm$fterm)) grouping_vars <- c(grouping_vars, sm$fterm)
    if (is_re_class && any(grouping_vars %in% .random_factor_set)) {
      .random_smooth_labels <- c(.random_smooth_labels, sm$label)
    }
  }
  .random_smooth_labels <- unique(.random_smooth_labels)
  cat("Random smooth labels in model: ",
      paste(.random_smooth_labels, collapse = " | "), "\n", sep = "")

  # Label-classifier: exact match against the model's random labels first,
  # then prefix match (catches gam.vcomp's "<label>1", "<label>2" suffixes
  # for fs.interaction smooths), then a final substring fallback for any
  # label-format we haven't anticipated.
  is_random_label <- function(lbl) {
    if (is.null(lbl) || !is.character(lbl) || nchar(lbl) == 0) return(FALSE)
    if (lbl %in% .random_smooth_labels) return(TRUE)
    for (rl in .random_smooth_labels) {
      if (startsWith(lbl, rl)) return(TRUE)
    }
    grepl("animal_id|plate", lbl)
  }
  rand_idx  <- which(sapply(smooth_labels, is_random_label))
  fixed_idx <- which(!sapply(smooth_labels, is_random_label))

  conc_summary <- if (!is.null(conc)) {
    worst   <- conc$worst
    diag(worst) <- NA   # ignore self-concurvity
    max_overall <- suppressWarnings(max(worst, na.rm = TRUE))
    max_overall <- if (is.finite(max_overall)) max_overall else NA_real_
    if (length(rand_idx) && length(fixed_idx)) {
      sub_rf <- worst[rand_idx, fixed_idx, drop = FALSE]
      max_rf <- suppressWarnings(max(sub_rf, na.rm = TRUE))
      max_rf <- if (is.finite(max_rf)) max_rf else NA_real_
    } else {
      max_rf <- NA_real_
    }
    list(max_overall = max_overall, max_random_vs_fixed = max_rf,
         n_random_terms = length(rand_idx), n_fixed_terms = length(fixed_idx))
  } else {
    list(error = "concurvity failed")
  }

  # gam.vcomp: standard deviations of each smooth's random effect (where
  # applicable). The return shape depends on the smoothness selection
  # method (mgcv source, mgcv.r ~line 4314-4338):
  #   * fREML/REML/ML with Hessian → "gam.vcomp" S3 list:
  #       $vc        matrix with cols (std.dev, lower, upper),
  #                  rownames = variance-component names (smoothing-param
  #                  labels + optional "scale" row)
  #       $all       optional named numeric vector
  #       $rank,$rank.hess,$conf.lev
  #   * Without Hessian → list($vc=<sqrt sps>, $all=<sqrt full.sp>) OR a
  #     bare numeric vector
  # NOTE: for fs.interaction smooths gam.vcomp returns multiple rows per
  # smooth (one per smoothing parameter); rownames are decorated with the
  # smoothing-param index. The prefix-match in is_random_label() picks
  # these up.
  vcomp_summary <- tryCatch({
    vc_out <- mgcv::gam.vcomp(final_model, rescale = FALSE)
    if (is.list(vc_out) && !is.null(vc_out$vc)) {
      vc_inner <- vc_out$vc
      if (is.matrix(vc_inner)) {
        sd_vec   <- vc_inner[, 1]
        sd_names <- rownames(vc_inner)
      } else {
        sd_vec   <- as.numeric(vc_inner)
        sd_names <- names(vc_inner)
      }
    } else if (is.matrix(vc_out)) {
      sd_vec   <- vc_out[, 1]
      sd_names <- rownames(vc_out)
    } else {
      sd_vec   <- as.numeric(vc_out)
      sd_names <- names(vc_out)
    }
    cat("gam.vcomp labels: ",
        paste(sd_names, collapse = " | "), "\n", sep = "")
    rand_sd_idx  <- which(sapply(sd_names, is_random_label))
    fixed_sd_idx <- which(!sapply(sd_names, is_random_label))
    resid_scale <- fit_quality$tweedie_phi
    list(
      per_smooth_sd        = setNames(as.list(sd_vec), sd_names),
      n_random_components  = length(rand_sd_idx),
      mean_random_sd       = if (length(rand_sd_idx))
                               mean(sd_vec[rand_sd_idx], na.rm = TRUE) else NA_real_,
      mean_fixed_smooth_sd = if (length(fixed_sd_idx))
                               mean(sd_vec[fixed_sd_idx], na.rm = TRUE) else NA_real_,
      random_sd_to_phi     = if (length(rand_sd_idx) && !is.na(resid_scale) && resid_scale > 0)
                               mean(sd_vec[rand_sd_idx], na.rm = TRUE) / sqrt(resid_scale)
                             else NA_real_
    )
  }, error = function(e) list(error = conditionMessage(e)))

  # EDF ratio: random / fixed-smooth. The fixed-effect smooths exclude the
  # global s(time_*) and the by-smooth (which is what we want to compare
  # the random part against).
  edf_summary <- tryCatch({
    per_smooth_edf <- smry$s.table[, "edf"]
    s_names <- rownames(smry$s.table)
    rand_e_idx <- which(sapply(s_names, is_random_label))
    edf_random <- sum(per_smooth_edf[rand_e_idx], na.rm = TRUE)
    edf_fixed  <- sum(per_smooth_edf[-rand_e_idx], na.rm = TRUE)
    list(
      edf_random        = edf_random,
      edf_fixed_smooth  = edf_fixed,
      random_to_fixed   = if (edf_fixed > 0) edf_random / edf_fixed else NA_real_,
      per_smooth_edf    = setNames(as.list(per_smooth_edf), s_names)
    )
  }, error = function(e) list(error = conditionMessage(e)))

  # ── Contrast variance summary (per Test_Family) ──────────────────────────
  # Uses raw p-values and SE from the master table built above.
  contrast_summary <- tryCatch({
    families <- unique(as.character(final_master_table$Test_Family))
    setNames(lapply(families, function(fam) {
      rows <- final_master_table[final_master_table$Test_Family == fam, , drop = FALSE]
      list(
        n_contrasts = nrow(rows),
        mean_SE     = mean(rows$SE,                 na.rm = TRUE),
        median_SE   = median(rows$SE,               na.rm = TRUE),
        max_SE      = max(rows$SE,                  na.rm = TRUE),
        mean_abs_est = mean(abs(rows$estimate),     na.rm = TRUE)
      )
    }), families)
  }, error = function(e) list(error = conditionMessage(e)))

  # ── vcov diagnostic on one focal emmeans grid (Full interaction) ─────────
  vcov_diag <- tryCatch({
    em_design <- emmeans(final_model,
                         specs = as.formula(paste0(
                           "~ ", paste(var_names, collapse = " * "),
                           if (abl_arch != "A") " | Group_factor" else ""
                         )))
    V <- vcov(em_design)
    diag(V) <- pmax(diag(V), .Machine$double.eps)  # guard against zero
    Vcor <- cov2cor(V)
    off  <- Vcor; diag(off) <- NA
    list(
      n_emm         = nrow(V),
      condition_num = kappa(V),
      max_abs_offdiag_cor = suppressWarnings(max(abs(off), na.rm = TRUE))
    )
  }, error = function(e) list(error = conditionMessage(e)))

  # ── Assemble + write ─────────────────────────────────────────────────────
  run_metrics <- list(
    ablation = list(
      arch              = abl_arch,
      random_basis      = random_basis,
      by_variable       = chosen_by_choice,
      ar_start_mode     = abl_ar_start,
      group_index       = if (nzchar(abl_group_index)) as.integer(abl_group_index) else NA_integer_,
      rho_fixed         = if (nzchar(abl_rho_fixed)) as.numeric(abl_rho_fixed) else NA_real_,
      n_obs             = nrow(gam_df),
      n_animals         = length(unique(as.character(gam_df$animal_id))),
      n_plates          = n_plates
    ),
    fit_quality        = fit_quality,
    adequacy           = adequacy,
    concurvity         = conc_summary,
    variance_components = vcomp_summary,
    edf                = edf_summary,
    contrast_variance  = contrast_summary,
    vcov_diagnostic    = vcov_diag
  )
  jsonlite::write_json(run_metrics,
                       file.path(output_dir, "run_metrics.json"),
                       auto_unbox = TRUE, pretty = TRUE,
                       na = "null")
  cat("Saved: run_metrics.json\n")
}

# =============================================================================
# 4. POSTERIOR EQUIVALENCE TESTING (Bayesian)
# =============================================================================
# Optional. Runs when contrast_selection.json sets posterior_equivalence.enabled.
# For each kept pair (lhs, rhs) and each phase group, draws from the BAM
# posterior to compute
#     M = max_t |η_lhs(t) − η_rhs(t)|     (log_e fold change at each draw)
# and converts to log_2 by dividing by log(2). Reports the posterior of M and
# Pr(M < δ) — a direct Bayesian credibility statement of equivalence within ±δ.
# Implementation uses mgcv's lpmatrix posterior + MASS::mvrnorm coefficient
# draws. (gratia is loaded for users who want to do additional posterior work.)

# Helper: for one (model, full_df, lhs, rhs, group_id), return the posterior of M.
# Unified model produces per-(pair × group) draws by varying Condition_Combo
# AND Group_factor in the newdata frame.
.posterior_equivalence_pair <- function(model, full_df, lhs, rhs, var_names,
                                        group_id, delta_log2, n_draws,
                                        n_time = 200, seed = 42) {
  cond_levels  <- levels(full_df$Condition_Combo)
  group_levels <- levels(full_df$Group_factor)
  group_str    <- as.character(group_id)
  if (!(lhs %in% cond_levels) || !(rhs %in% cond_levels)) return(NULL)
  if (!(group_str %in% group_levels))                     return(NULL)

  parse_cond <- function(cond_str) {
    parts <- strsplit(cond_str, "+", fixed = TRUE)[[1]]
    if (length(parts) != length(var_names)) {
      stop("Cannot parse '", cond_str, "': expected ", length(var_names),
           " parts (one per variable), got ", length(parts))
    }
    setNames(as.list(parts), var_names)
  }
  lhs_vars <- parse_cond(lhs)
  rhs_vars <- parse_cond(rhs)

  # Use the time_in_group range for THIS phase group so the trajectory is
  # built only over the phase the comparison lives in.
  group_rows <- full_df[full_df$Group == group_id, ]
  if (nrow(group_rows) == 0) return(NULL)
  time_grid <- seq(0, max(group_rows$time_in_group), length.out = n_time)
  # time_sec is needed for the random-effect smooths even though those are
  # excluded from the linear predictor below — predict() still requires the
  # column to exist. Anchor it at the midpoint of the group's range.
  time_sec_anchor <- mean(range(group_rows$time_sec))

  build_nd <- function(cond_str, vars_list) {
    nd <- data.frame(
      time_in_group = time_grid,
      time_sec      = time_sec_anchor
    )
    for (v in var_names) {
      nd[[v]] <- factor(vars_list[[v]], levels = levels(full_df[[v]]))
    }
    nd$Condition_Combo <- factor(cond_str,  levels = cond_levels)
    nd$Group_factor    <- factor(group_str, levels = group_levels)
    # animal_id and plate are required by predict() but excluded from the
    # linear predictor (population-level trajectories).
    nd$animal_id <- factor(full_df$animal_id[1],
                           levels = levels(full_df$animal_id))
    nd$plate <- factor(full_df$plate[1],
                       levels = levels(full_df$plate))
    nd
  }
  nd_lhs <- build_nd(lhs, lhs_vars)
  nd_rhs <- build_nd(rhs, rhs_vars)

  # Linear-predictor design matrices, excluding the random-effect smooths so
  # we get population-level (not animal- or plate-specific) trajectories.
  # The exact term labels depend on the chosen random_basis:
  #   none     → no random-effect terms to exclude
  #   re       → s(animal_id), s(plate)
  #   re_slope → both the intercept terms AND the time slope terms
  #   fs / sz  → only the time-in-group factor-smooth terms
  excl <- switch(random_basis,
    none     = character(0),
    re       = c("s(animal_id)", "s(plate)"),
    re_slope = c("s(animal_id)", "s(plate)",
                 "s(time_in_group,animal_id)", "s(time_in_group,plate)"),
    fs       = c("s(time_in_group,animal_id)", "s(time_in_group,plate)"),
    sz       = c("s(time_in_group,animal_id)", "s(time_in_group,plate)")
  )
  X_lhs <- predict(model, newdata = nd_lhs, type = "lpmatrix", exclude = excl)
  X_rhs <- predict(model, newdata = nd_rhs, type = "lpmatrix", exclude = excl)
  X_diff <- X_lhs - X_rhs                      # n_time × n_coef

  # Bayesian (large-sample Gaussian) posterior of β
  set.seed(seed)
  Cv <- vcov(model, unconditional = TRUE)
  beta_samples <- MASS::mvrnorm(n_draws, mu = coef(model), Sigma = Cv)  # n_draws × n_coef

  # Posterior draws of η_lhs(t) − η_rhs(t) on log_e scale → divide by log(2)
  diff_samples <- X_diff %*% t(beta_samples)   # n_time × n_draws
  M_log2 <- apply(abs(diff_samples), 2, max) / log(2)

  list(
    M_log2  = M_log2,
    median  = median(M_log2),
    lower95 = unname(quantile(M_log2, 0.025)),
    upper95 = unname(quantile(M_log2, 0.975)),
    pr_equiv = mean(M_log2 < delta_log2)
  )
}

if (abl_mode && abl_skip_post) {
  cat("\nAblation mode: skipping Bayesian posterior equivalence testing.\n")
} else if (!is.null(contrast_spec) &&
    !is.null(contrast_spec$posterior_equivalence) &&
    isTRUE(contrast_spec$posterior_equivalence$enabled) &&
    length(contrast_spec$kept_pairs) > 0) {

  delta_log2 <- as.numeric(contrast_spec$posterior_equivalence$delta %||% 1.0)
  n_draws    <- as.integer(contrast_spec$posterior_equivalence$n_draws %||% 10000)

  cat("\n========================================\n")
  cat("Posterior Equivalence Testing (Bayesian)\n")
  cat(sprintf("  delta = %.3f (log_2 fold change)\n", delta_log2))
  cat(sprintf("  posterior draws = %d\n", n_draws))
  cat("========================================\n")

  pe_summary_rows <- list()
  pe_draws_named  <- list()

  # Unified-model PE: iterate over (group, pair) directly against the single
  # final_model. Prediction grids carry both Condition_Combo and Group_factor.
  for (g in unique_groups) {
    g_str <- as.character(g)
    cat(sprintf("Group %s: posterior equivalence...\n", g_str))

    for (i in seq_along(contrast_spec$kept_pairs)) {
      pair <- contrast_spec$kept_pairs[[i]]
      lhs <- pair[[1]]; rhs <- pair[[2]]

      pe <- tryCatch(
        .posterior_equivalence_pair(final_model, gam_df, lhs, rhs, var_names,
                                    group_id = g, delta_log2 = delta_log2,
                                    n_draws = n_draws, seed = 42L + i),
        error = function(e) {
          warning(sprintf("Group %s, %s vs %s: %s", g_str, lhs, rhs,
                          conditionMessage(e)))
          NULL
        }
      )
      if (is.null(pe)) next

      pe_summary_rows[[length(pe_summary_rows) + 1]] <- data.frame(
        Group         = as.integer(g_str),
        lhs           = lhs,
        rhs           = rhs,
        pair          = paste(lhs, "vs", rhs),
        M_median_log2 = pe$median,
        M_lower95     = pe$lower95,
        M_upper95     = pe$upper95,
        Pr_equiv      = pe$pr_equiv,
        delta_log2    = delta_log2,
        stringsAsFactors = FALSE
      )
      pe_draws_named[[paste0("g", g_str, "__", lhs, "__", rhs)]] <- pe$M_log2
    }
  }

  if (length(pe_summary_rows) > 0) {
    pe_summary_df <- bind_rows(pe_summary_rows) %>%
      arrange(Group, pair) %>%
      as_tibble()
    cat("\n--- POSTERIOR EQUIVALENCE SUMMARY ---\n")
    # Print via as.data.frame to avoid a partial-arg-matching crash when the
    # table is a plain data.frame instead of a tibble.
    print(as.data.frame(pe_summary_df), row.names = FALSE)
    write_csv(pe_summary_df, "posterior_equivalence_summary.csv")

    # Long form of the raw posterior draws — persisted so the Python
    # Visualisations tab can render the PE density and heatmap interactively.
    draws_df <- bind_rows(lapply(seq_along(pe_draws_named), function(j) {
      key <- names(pe_draws_named)[j]
      m <- regmatches(key, regexec("^g(.+?)__(.+?)__(.+)$", key))[[1]]
      data.frame(
        Group  = as.integer(m[2]),
        pair   = paste(m[3], "vs", m[4]),
        M_log2 = pe_draws_named[[j]],
        stringsAsFactors = FALSE
      )
    }))
    write_csv(draws_df, "posterior_equivalence_draws.csv")
    cat("Saved: posterior_equivalence_draws.csv (",
        nrow(draws_df), " draws total)\n", sep = "")
  } else {
    cat("No posterior-equivalence results were produced (all pairs skipped).\n")
  }
}

cat("\nAnalysis complete. All outputs saved to:", output_dir, "\n")
