# =============================================================================
# Ablation Study — formula variants for the BAM movement model
# -----------------------------------------------------------------------------
# Purpose: find a formula that is FUNCTIONALLY similar to the production formula
# (same fitted movement trajectories) but far cheaper to fit — so the initial
# screen (242 conditions × 658 animals, which OOMs at 1.5 TB with the production
# formula) becomes tractable.
#
# Design (mirrors the confirmation-study ablation methodology):
#   Stage 1  per-variant metric panel  -> run_metrics_<key>.json
#            (AIC/BIC/fREML/logLik/dev_explained/r_sq_adj/phi/edf, k.check
#             pass-rate, AR1 residual lag-1, concurvity, gam.vcomp,
#             edf.random_to_fixed) + wall-clock + peak RSS
#   Stage 1  adequacy filters          -> eligible flag (catastrophic fails out)
#   Stage 2  Pareto front              -> on_front flag (within this run)
#   Stage 3  5-fold animal-stratified CV Tweedie deviance  -> cv_deviance (GATE)
#   FINAL    winner = best (similarity_to_baseline + speedup) among CV-passers
#
#   similarity_to_baseline (FUNCTIONAL): predict every variant AND the baseline
#   on one shared (condition × time) grid (population level, random smooths
#   excluded), then compare fitted values on the response scale:
#     fit_corr (Pearson), fit_rmse, fit_ccc (Lin's concordance),
#     plus per-condition Δη(t) effect-curve distance vs baseline.
#
# Reads the already-aggregated long frame (aggregated_data.csv) so it does NOT
# re-parse the raw Zebrabox .xls files — aggregation is identical across
# variants, so we do it once (in the main pipeline) and reuse it here.
#
# rho is FIXED (not the two-pass optimisation) to keep each fit cheap; this is
# the same compute-saving deviation FamilySelection.R uses. Set via env.
#
# Usage (env-driven, container-friendly):
#   AGG_CSV=/path/aggregated_data.csv  OUT_DIR=/path/out \
#   VAR_NAMES=Genotype,Drug  REF_VALUES=WT,DMSO  REF_CONDITION=WT+DMSO \
#   ABL_RHO=0.2  ABL_KCHECK=1  Rscript ablation_study.R
# =============================================================================


# ── CONFIG (all env-overridable; local-friendly defaults) ────────────────────
AGG_CSV   <- Sys.getenv("ABL_AGG_CSV",
                        unset = "Data/Confirmation/Stand-Alone Output/aggregated_data.csv")
OUT_DIR   <- Sys.getenv("ABL_OUT_DIR", unset = "Data/Confirmation/Ablation Output")
VAR_NAMES_RAW  <- Sys.getenv("ABL_VAR_NAMES",     unset = "Genotype,Drug")
REF_VALUES_RAW <- Sys.getenv("ABL_REF_VALUES",    unset = "WT,DMSO")
REF_CONDITION  <- Sys.getenv("ABL_REF_CONDITION", unset = "WT+DMSO")
# Fixed AR(1) rho — skip the two-pass optimisation to halve compute.
ABL_RHO   <- as.numeric(Sys.getenv("ABL_RHO", unset = "0.2"))
# How the `param` condition level is parametrized (the design differs per
# dataset — see plan): "interaction" -> Genotype*Drug (balanced confirmation
# data); "condition" -> Condition_Combo as a single factor (initial screen,
# where Genotype*Drug is non-estimable because only DMSO is shared genotypes).
ABL_PARAM_MODE <- Sys.getenv("ABL_PARAM_MODE", unset = "interaction")
# Optional: subset the aggregated frame to plates whose name contains this
# string (e.g. "Plate 2") before fitting. Empty = use all plates.
ABL_PLATE_FILTER <- Sys.getenv("ABL_PLATE_FILTER", unset = "")
# Optional: keep only wells whose loc_coord row-letter is in this set (e.g.
# "ABCD" to use rows A-D and drop E-H). Empty = all rows. A tractability lever
# for the initial screen (fewer conditions -> smaller smooth -> faster).
ABL_ROW_FILTER <- Sys.getenv("ABL_ROW_FILTER", unset = "")
# Which variants to run (comma list of keys, or "all").
ABL_ONLY  <- Sys.getenv("ABL_ONLY", unset = "all")
# Run CV gate (expensive). 1=yes, 0=skip.
ABL_CV    <- Sys.getenv("ABL_CV", unset = "1") %in% c("1", "TRUE", "T", "yes")
# CV fold unit: "animal" (default, balanced data like confirmation) or
# "replicate" (initial screen — each condition is ~1 animal/replicate, so
# animal-folding drops whole conditions from training; fold by replicate keeps
# every condition represented, since each condition appears in >=2 replicates).
CV_FOLD_BY <- Sys.getenv("ABL_CV_FOLD_BY", unset = "animal")
N_CV_FOLDS <- 5L
CV_SEED    <- 42L
# Common prediction grid resolution for functional similarity.
SIM_N_TIME <- 200L


# ── PACKAGE LIBRARY + INSTALL (matches main script's container logic) ────────
.conda_prefix <- Sys.getenv("CONDA_PREFIX")
.lib <- if (nzchar(.conda_prefix)) {
  file.path(.conda_prefix, "lib", "R", "library")
} else {
  file.path(Sys.getenv("HOME", unset = "/tmp"), "R", "library")
}
dir.create(.lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(.lib, .libPaths()))

required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel", "jsonlite",
                   "RhpcBLASctl")
cran_mirror <- Sys.getenv("CRAN_MIRROR", unset = "https://cloud.r-project.org")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1),
                                 quietly = TRUE)]
if (length(missing) > 0) {
  installable <- setdiff(missing, "parallel")
  if (length(installable) > 0) {
    cat("Installing: ", paste(installable, collapse = ", "), "\n", sep = "")
    install.packages(installable, repos = cran_mirror)
  }
}
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(mgcv); library(jsonlite)
})
options(width = 320); options(warn = 1)
`%||%` <- function(a, b) if (is.null(a)) b else a

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── THREADING (same OpenBLAS-cap logic as the main script) ───────────────────
OPENBLAS_BUILD_CAP <- 128L
n_cores        <- max(1L, parallel::detectCores())
thread_ceiling <- max(1L, OPENBLAS_BUILD_CAP - 8L)
FIT_PERF_CAP   <- 32L
usable_cores   <- max(1L, min(n_cores - 1L, thread_ceiling, FIT_PERF_CAP))
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
}


# ── PEAK-RSS HELPER (Linux /proc; matches main script) ───────────────────────
peak_rss_mb <- function() {
  if (file.exists("/proc/self/status")) {
    hwm <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
    if (length(hwm) == 0) return(NA_real_)
    as.numeric(sub(".*?(\\d+).*", "\\1", hwm)) / 1024
  } else NA_real_
}

# Reset the process peak-RSS high-water mark (VmHWM) so the NEXT peak_rss_mb()
# reading reflects only the memory used AFTER this call — not the cumulative
# maximum since process start. Without this, every cell inherits the run's
# running maximum and per-cell memory is unmeasurable (all cells read the same
# inflated peak). Linux-only (/proc/self/clear_refs); no-op elsewhere.
reset_peak_rss <- function() {
  gc(full = TRUE)  # release freed memory first so the reset baseline is clean
  if (file.exists("/proc/self/clear_refs")) {
    tryCatch(writeLines("5", "/proc/self/clear_refs"), error = function(e) NULL)
  }
  invisible(NULL)
}


# ── BUILD gam_df FROM THE AGGREGATED FRAME (identical to main script) ─────────
var_names_display <- str_split(VAR_NAMES_RAW, ",")[[1]]
ref_values        <- str_split(REF_VALUES_RAW, ",")[[1]]
var_names         <- make.names(var_names_display, unique = TRUE)
ref_map           <- setNames(ref_values, var_names)

cat("Reading aggregated frame: ", AGG_CSV, "\n", sep = "")
if (!file.exists(AGG_CSV)) stop("aggregated_data.csv not found: ", AGG_CSV)
full_df <- fread(AGG_CSV)
setDT(full_df)

# Optional plate subset (Phase B: e.g. "Plate 2"). Done BEFORE factor build so
# dropped plates don't leave empty factor levels. Match on the plate string,
# which in the aggregated frame looks like "Plate 2_Replicate 1".
if (nzchar(ABL_PLATE_FILTER)) {
  before <- nrow(full_df)
  full_df <- full_df[grepl(ABL_PLATE_FILTER, plate, fixed = TRUE)]
  cat("Plate filter '", ABL_PLATE_FILTER, "': kept ", nrow(full_df), " / ",
      before, " rows.\n", sep = "")
  if (nrow(full_df) == 0L) stop("Plate filter matched no rows: ", ABL_PLATE_FILTER)
}

# Optional well-row subset (e.g. "ABCD" keeps rows A-D, drops E-H). loc_coord
# looks like "A1","B7"; its first char is the row letter. Done before factor
# build so dropped conditions don't leave empty Condition_Combo levels.
if (nzchar(ABL_ROW_FILTER)) {
  before <- nrow(full_df)
  keep_rows <- strsplit(ABL_ROW_FILTER, "", fixed = TRUE)[[1]]
  full_df <- full_df[substr(as.character(loc_coord), 1L, 1L) %in% keep_rows]
  cat("Row filter '", ABL_ROW_FILTER, "': kept ", nrow(full_df), " / ",
      before, " rows.\n", sep = "")
  if (nrow(full_df) == 0L) stop("Row filter matched no rows: ", ABL_ROW_FILTER)
}

# The aggregated frame already has Condition, plate, animal_id, Group, time_sec,
# pixel_diff and the split var columns (Genotype, Drug). Rebuild factors +
# time_in_group + start_event exactly as the production script does.
full_df[, `:=`(
  plate           = as.factor(plate),
  animal_id       = as.factor(animal_id),
  Condition_Combo = as.factor(Condition),
  Group_factor    = as.factor(Group)
)]
full_df[, time_in_group := as.integer(time_sec - min(time_sec)), by = Group]
setorder(full_df, animal_id, Group, time_in_group)
full_df[, start_event := (seq_len(.N) == 1L) |
                        (Group != shift(Group, fill = -1L)), by = animal_id]
gam_df <- full_df

for (v in var_names) {
  gam_df[[v]] <- as.factor(gam_df[[v]])
  if (ref_map[[v]] %in% levels(gam_df[[v]])) {
    gam_df[[v]] <- relevel(gam_df[[v]], ref = ref_map[[v]])
  }
}
if (REF_CONDITION %in% levels(gam_df$Condition_Combo)) {
  gam_df$Condition_Combo <- relevel(gam_df$Condition_Combo, ref = REF_CONDITION)
}

n_groups <- length(unique(gam_df$Group))
n_plates <- length(unique(gam_df$plate))
use_group_factor <- n_groups >= 2L
k_prod <- 30L  # production k

cat("Rows: ", nrow(gam_df), " | animals: ",
    nlevels(gam_df$animal_id), " | conditions: ",
    nlevels(gam_df$Condition_Combo), " | plates: ", n_plates,
    " | phase groups: ", n_groups, "\n", sep = "")


# ── FORMULA BUILDERS ─────────────────────────────────────────────────────────
# A variant = (cond_mode, cond_k, animal_mode, animal_k). The shared time smooth
# is constant across variants. The parametric part is constant EXCEPT for the
# `param` condition level, whose form depends on ABL_PARAM_MODE (see below).
parametric_str <- if (use_group_factor) {
  paste(c(var_names, "Group_factor"), collapse = " * ")
} else {
  paste(var_names, collapse = " * ")
}
by_term <- if (use_group_factor) {
  "interaction(Condition_Combo, Group_factor)"
} else "Condition_Combo"
plate_smooth <- if (n_plates >= 2) " + s(time_in_group, plate, bs = 'sz')" else ""

# Parametric expansion used when the condition effect is dropped to a parametric
# term (cond_mode == "param"). The estimable form differs by dataset:
#   interaction -> Genotype*Drug (balanced confirmation: full 2x2)
#   condition   -> Condition_Combo as ONE factor (initial: Genotype*Drug aliased)
param_cond_str <- switch(ABL_PARAM_MODE,
  interaction = parametric_str,
  condition   = if (use_group_factor) "Condition_Combo * Group_factor" else "Condition_Combo",
  stop("ABL_PARAM_MODE must be 'interaction' or 'condition'"))

build_formula <- function(cond_mode, cond_k, animal_mode, animal_k,
                           base_k = k_prod) {
  # base_k = k of the SHARED baseline time smooth s(time_in_group, k=base_k) —
  # the population-average trajectory common to all conditions. Defaults to
  # k_prod (30) so existing callers are unchanged; the baseline-k scan overrides
  # it. (Distinct from cond_k, which is the per-condition deviation smooth.)
  param_str <- if (cond_mode == "param") param_cond_str else parametric_str
  # condition effect on the time-course
  cond_term <- switch(cond_mode,
    smooth = sprintf(" + s(time_in_group, by = %s, k = %d)", by_term, cond_k),
    param  = "",                       # condition shape dropped; parametric only
    stop("bad cond_mode"))
  # per-animal effect
  animal_term <- switch(animal_mode,
    sz = sprintf(" + s(time_in_group, animal_id, bs = 'sz', k = %d)", animal_k),
    re = " + s(animal_id, bs = 're')",  # random intercept per animal
    stop("bad animal_mode"))
  as.formula(paste0(
    "pixel_diff ~ ", param_str,
    " + s(time_in_group, k = ", base_k, ")",
    cond_term, plate_smooth, animal_term))
}

# ── VARIANT GRID (Phase A, confirmation) — MATCHED 2×N re-vs-sz10 ─────────────
# Prior run findings drove this design:
#   * condparam dropped — fit_ccc ~0.64 (parametric destroys trajectory shape).
#   * animal sz axis is flat (sz10≈sz5≈sz3); sz10 was the ONLY significantly
#     better-than-baseline cell (paired t=-3.9), so sz10 is the sz representative.
#   * re cells were big-speedup but "within noise" on CV (high se_diff).
# To cleanly separate the ANIMAL-TERM effect (re vs sz10) from the CONDITION-k
# effect, run a MATCHED 2×N factorial: both animal terms at EACH cond-k. Then a
# paired per-fold test compares re vs sz10 AT each cond level (same folds), with
# no confound. cond{30,20,10,5} × {re, sz10} + baseline = 9 cells.
# Similarity reference cell. Defaults to the production baseline, but is
# env-overridable: at initial-screen scale the full-k baseline is compute-
# intractable (select=TRUE serial penalty loop over 82 conditions), so it can't
# serve as the reference there. Setting ABL_BASELINE to a tractable cell (e.g.
# cond30__animalre) makes that cell the reference; all similarity (fit_ccc,
# est_corr) is then measured relative to it.
BASELINE_KEY <- Sys.getenv("ABL_BASELINE", unset = "cond30__animalsz30")
# BASELINE-k × COND-k FACTORIAL: the two time-smooths — shared baseline
# s(time_in_group, k=base_k) and per-condition s(time_in_group, by=Cond, k=cond_k)
# — may INTERACT. base_k 40≈50 was equivalent at cond_k=30; this 3×3 grid checks
# whether that holds at cond_k=40 and 50, isolating any cross-term. animal=re.
VARIANTS <- list()
for (bk in c(50, 40, 30)) for (ck in c(50, 40, 30)) {
  VARIANTS[[sprintf("base%d__cond%d__re", bk, ck)]] <-
    list(cond_mode="smooth", cond_k=ck, animal_mode="re", animal_k=0, base_k=bk)
}
# Put the chosen similarity reference first (default ABL_BASELINE = cond30__animalre
# for this scan — see launcher). Falls back gracefully if not in the grid.
if (BASELINE_KEY %in% names(VARIANTS))
  VARIANTS <- VARIANTS[c(BASELINE_KEY, setdiff(names(VARIANTS), BASELINE_KEY))]

if (!identical(ABL_ONLY, "all")) {
  keep <- str_split(ABL_ONLY, ",")[[1]] |> trimws()
  VARIANTS <- VARIANTS[intersect(c(BASELINE_KEY, keep), names(VARIANTS))]
}
cat("Variant grid: ", length(VARIANTS), " cells (param_mode=", ABL_PARAM_MODE,
    ").\n", sep = "")


# ── A SINGLE FIT WITH FIXED RHO ──────────────────────────────────────────────
fit_variant <- function(spec) {
  f <- build_formula(spec$cond_mode, spec$cond_k, spec$animal_mode, spec$animal_k, base_k = spec$base_k %||% k_prod)
  bam(formula = f, data = gam_df, family = tw(),
      rho = ABL_RHO, AR.start = gam_df$start_event,
      select = TRUE, method = "fREML", discrete = TRUE,
      nthreads = usable_cores, chunk.size = 10000)
}


# ── STAGE 1 METRIC PANEL ─────────────────────────────────────────────────────
metric_panel <- function(model) {
  kc <- tryCatch(mgcv::k.check(model), error = function(e) NULL)
  kc_pass <- if (!is.null(kc)) mean(kc[, "p-value"] > 0.05, na.rm = TRUE) else NA_real_
  conc <- tryCatch(mgcv::concurvity(model, full = FALSE), error = function(e) NULL)
  conc_max <- if (!is.null(conc)) max(conc$worst, na.rm = TRUE) else NA_real_
  vc <- tryCatch(mgcv::gam.vcomp(model, rescale = FALSE), error = function(e) NULL)
  ar1 <- tryCatch(as.numeric(acf(residuals(model), plot = FALSE)$acf[2]),
                  error = function(e) NA_real_)
  sm <- summary(model)
  list(
    aic            = tryCatch(AIC(model),  error = function(e) NA_real_),
    bic            = tryCatch(BIC(model),  error = function(e) NA_real_),
    freml          = tryCatch(model$gcv.ubre[[1]], error = function(e) NA_real_),
    loglik         = tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_),
    dev_explained  = tryCatch(sm$dev.expl, error = function(e) NA_real_),
    r_sq_adj       = tryCatch(sm$r.sq,    error = function(e) NA_real_),
    tweedie_phi    = tryCatch(sm$dispersion, error = function(e) NA_real_),
    total_edf      = tryCatch(sum(model$edf), error = function(e) NA_real_),
    n_coef         = length(coef(model)),
    k_index_pass   = kc_pass,
    ar1_resid_lag1 = ar1,
    concurvity_max = conc_max,
    converged      = isTRUE(model$converged),       # surfaced in the summary
    conv_warn      = isTRUE(!model$converged)
  )
}


# ── FUNCTIONAL SIMILARITY: shared (condition × time) prediction grid ─────────
# Population-level predictions (random smooths excluded), compared on the
# response scale. The grid is built once; every model predicts on it.
build_sim_grid <- function() {
  conds <- levels(gam_df$Condition_Combo)
  tg <- seq(0, max(gam_df$time_in_group), length.out = SIM_N_TIME)
  tsec <- mean(range(gam_df$time_sec))
  grid <- CJ(Condition_Combo = factor(conds, levels = conds),
             time_in_group = tg, sorted = FALSE)
  grid[, time_sec := tsec]
  # split Condition back into the design vars (so parametric terms resolve)
  parts <- tstrsplit(as.character(grid$Condition_Combo), "+", fixed = TRUE)
  for (i in seq_along(var_names)) {
    grid[[var_names[i]]] <- factor(parts[[i]], levels = levels(gam_df[[var_names[i]]]))
  }
  if (use_group_factor) {
    grid[, Group_factor := factor(levels(gam_df$Group_factor)[1],
                                  levels = levels(gam_df$Group_factor))]
  }
  grid[, animal_id := factor(levels(gam_df$animal_id)[1],
                             levels = levels(gam_df$animal_id))]
  grid[, plate := factor(levels(gam_df$plate)[1], levels = levels(gam_df$plate))]
  grid[]
}
# terms to exclude so predictions are population-level (no random smooths)
excl_terms <- function(spec) {
  e <- if (spec$animal_mode == "sz") "s(time_in_group,animal_id)" else "s(animal_id)"
  if (n_plates >= 2) e <- c(e, "s(time_in_group,plate)")
  e
}
ccc <- function(x, y) {            # Lin's concordance correlation coefficient
  mx <- mean(x); my <- mean(y); vx <- var(x); vy <- var(y); cxy <- cov(x, y)
  2 * cxy / (vx + vy + (mx - my)^2)
}

# ── INFERENTIAL CONTRASTS: each condition vs the reference (DMSO), per model ──
# Reproduces the *conclusion* each model reaches, in a way that works for EVERY
# variant (smooth or parametric) — unlike emmeans, which breaks when the
# by-condition smooth is dropped. Uses the same lpmatrix + vcov machinery as the
# posterior-equivalence block: for each condition C, contrast its time-averaged
# population trajectory against the reference condition's, on the link scale.
#   estimate = mean_t[ η_C(t) − η_ref(t) ]   (β-linear, so = mean rows of X_diff %*% β)
#   SE       = sqrt(w' V w)                  where w = colMeans(X_C − X_ref)
#   z, p     = Wald test of estimate = 0
# Returns a data.frame keyed by condition: estimate, p_value, significant.
SIG_THRESHOLD <- 0.05
contrast_vs_ref <- function(model, spec) {
  conds <- levels(gam_df$Condition_Combo)
  non_ref <- setdiff(conds, REF_CONDITION)
  if (length(non_ref) == 0) return(NULL)
  excl <- excl_terms(spec)
  beta <- coef(model)
  V    <- vcov(model, unconditional = TRUE)
  tg   <- seq(0, max(gam_df$time_in_group), length.out = SIM_N_TIME)
  tsec <- mean(range(gam_df$time_sec))
  # lpmatrix for one condition's trajectory (population level).
  Xrow <- function(cond) {
    nd <- data.table(Condition_Combo = factor(cond, levels = conds),
                     time_in_group = tg, time_sec = tsec)
    parts <- tstrsplit(cond, "+", fixed = TRUE)
    for (i in seq_along(var_names))
      nd[[var_names[i]]] <- factor(parts[[i]], levels = levels(gam_df[[var_names[i]]]))
    if (use_group_factor)
      nd[, Group_factor := factor(levels(gam_df$Group_factor)[1],
                                  levels = levels(gam_df$Group_factor))]
    nd[, animal_id := factor(levels(gam_df$animal_id)[1],
                             levels = levels(gam_df$animal_id))]
    nd[, plate := factor(levels(gam_df$plate)[1], levels = levels(gam_df$plate))]
    colMeans(predict(model, newdata = nd, type = "lpmatrix", exclude = excl))
  }
  x_ref <- Xrow(REF_CONDITION)
  out <- lapply(non_ref, function(cond) {
    w   <- Xrow(cond) - x_ref          # contrast weight vector (time-averaged)
    est <- as.numeric(w %*% beta)
    se  <- sqrt(as.numeric(t(w) %*% V %*% w))
    z   <- est / se
    p   <- 2 * pnorm(-abs(z))
    data.frame(condition = cond, estimate = est, p_value = p,
               significant = p < SIG_THRESHOLD, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


# ── RUN ALL VARIANTS ─────────────────────────────────────────────────────────
sim_grid <- build_sim_grid()
baseline_pred <- NULL
baseline_contr <- NULL          # baseline's per-condition contrasts vs DMSO
contr_all <- list()             # per-variant per-condition contrast tables (saved)
records <- list()

for (key in names(VARIANTS)) {
  spec <- VARIANTS[[key]]
  cat("\n──", key, "──\n")
  # Reset the peak-RSS high-water mark so rss1 below captures THIS fit's own peak
  # memory, not the run's cumulative maximum (see reset_peak_rss).
  reset_peak_rss()
  t0 <- Sys.time(); rss0 <- peak_rss_mb()
  model <- tryCatch(fit_variant(spec), error = function(e) {
    cat("  FIT FAILED: ", conditionMessage(e), "\n"); NULL
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rss1 <- peak_rss_mb()

  if (is.null(model)) {
    records[[key]] <- list(variant = key, spec = spec, fit_ok = FALSE,
                           wall_secs = elapsed)
    next
  }

  pred <- as.numeric(predict(model, newdata = sim_grid, type = "response",
                             exclude = excl_terms(spec)))
  contr <- tryCatch(contrast_vs_ref(model, spec), error = function(e) {
    cat("  contrast failed: ", conditionMessage(e), "\n"); NULL })
  if (key == BASELINE_KEY) { baseline_pred <- pred; baseline_contr <- contr }
  # Persist this variant's full per-condition contrast table (tagged by variant)
  # so we can later see exactly WHICH conditions flip significance between
  # formulas — essential for auditing the 82-condition initial screen.
  if (!is.null(contr)) { contr$variant <- key; contr_all[[key]] <- contr }

  rec <- c(list(variant = key, spec = spec, fit_ok = TRUE,
                wall_secs = elapsed, peak_rss_mb = rss1,
                rss_delta_mb = if (is.na(rss1) || is.na(rss0)) NA_real_ else rss1 - rss0),
           metric_panel(model))

  # (1) FUNCTIONAL similarity — do the models predict the same trajectories?
  if (!is.null(baseline_pred) && key != BASELINE_KEY) {
    ok <- is.finite(pred) & is.finite(baseline_pred)
    rec$fit_corr <- suppressWarnings(cor(pred[ok], baseline_pred[ok]))
    rec$fit_rmse <- sqrt(mean((pred[ok] - baseline_pred[ok])^2))
    rec$fit_ccc  <- ccc(pred[ok], baseline_pred[ok])
  }
  # (2) INFERENTIAL similarity — do the models reach the same CONCLUSIONS?
  # Compare per-condition contrast estimates + significance calls vs baseline.
  #   est_corr  = PRIMARY inferential metric. Correlation of the continuous
  #               effect ESTIMATES — threshold-free, so it does NOT flip at the
  #               p=0.05 boundary. This is the robust "same conclusions?" signal.
  #   sig_agree = SECONDARY, decision-level sanity check. Fraction of conditions
  #               where the binary significant/not call matches. It is BRITTLE at
  #               the boundary by construction (baseline p=0.049 vs variant
  #               p=0.051 counts as a disagreement even though the models agree),
  #               so it UNDERSTATES similarity for borderline conditions. Read it
  #               together with est_corr: high est_corr + imperfect sig_agree =
  #               same effects with boundary jitter (fine); low est_corr = the
  #               conclusions genuinely diverge (not fine).
  #   sig_recall= of baseline's significant hits, fraction the variant also flags.
  rec$n_sig <- if (!is.null(contr)) sum(contr$significant, na.rm = TRUE) else NA_integer_
  if (!is.null(baseline_contr) && !is.null(contr) && key != BASELINE_KEY) {
    m <- merge(baseline_contr, contr, by = "condition",
               suffixes = c("_base", "_var"))
    if (nrow(m) > 0) {
      rec$est_corr     <- suppressWarnings(cor(m$estimate_base, m$estimate_var))
      rec$sig_agree    <- mean(m$significant_base == m$significant_var)
      bs <- m$significant_base
      rec$sig_recall   <- if (any(bs)) mean(m$significant_var[bs]) else NA_real_
    }
  }

  records[[key]] <- rec
  write_json(rec, file.path(OUT_DIR, paste0("run_metrics_", key, ".json")),
             auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
  cat(sprintf("  wall=%.1fs  peakRSS=%s MB  nCoef=%d  dev=%.4f\n",
              elapsed, format(round(rec$peak_rss_mb)), rec$n_coef,
              rec$dev_explained %||% NA_real_))
}


# ── STAGE 1 ADEQUACY FILTERS (catastrophic-fail screen) ──────────────────────
phi_min <- min(vapply(records, function(r) r$tweedie_phi %||% NA_real_, numeric(1)),
               na.rm = TRUE)
eligible <- function(r) {
  # LOOSE early-stage screen — meant to drop only catastrophic failures, not to
  # be hard to pass (matches the spirit of the original confirmation-study
  # ablation). The decisive metrics at the end of Phase A are convergence,
  # functional similarity, and speedup — not these filters.
  #
  # NOTE on k.check: we are DELIBERATELY shrinking k (k=3,5,10), which fails
  # k.check by design — so the k-index filter is intentionally NOT a gate here
  # (it would disqualify exactly the cheap variants we want to evaluate). It is
  # still reported in the summary for inspection.
  if (!isTRUE(r$fit_ok)) return(FALSE)
  ok <- TRUE
  if (isTRUE(r$conv_warn)) ok <- FALSE                          # must converge (the key gate)
  if (!is.na(r$ar1_resid_lag1 %||% NA) && abs(r$ar1_resid_lag1 %||% 0) > 0.9) ok <- FALSE # only severe residual autocorr
  if (!is.na(r$tweedie_phi %||% NA) && is.finite(phi_min) &&
      (r$tweedie_phi %||% Inf) > 1000 * phi_min) ok <- FALSE     # only runaway dispersion
  ok
  # NOTE: concurvity is NOT a gate. With bs='sz' factor-smooth random effects,
  # concurvity ≈ 1 is expected and benign (the per-animal smooths overlap the
  # population smooth by construction). It is reported in the summary but does
  # not disqualify a variant.
}


# ── STAGE 3 GATE: 5-FOLD ANIMAL-STRATIFIED CV TWEEDIE DEVIANCE ────────────────
# Mirrors the proven CV logic in helpers/scripts/TweedieAR1_BAM.R (the original
# ablation). Three lessons borrowed from it that fixed silent NA/NaN deviance:
#   (a) predict() newdata MUST be a plain data.frame — mgcv's predict() fails
#       ("cannot xtfrm data frames") on tibble/data.table newdata, and the
#       failure is otherwise swallowed. coerce with as.data.frame().
#   (b) held-out animals/plates aren't in the training fit; since we exclude the
#       random smooths anyway, set them to a TRAINING placeholder level
#       (droplevels'd) so predict() has known coefs to ignore.
#   (c) stage predict vs dev.resids in separate trycatch and RECORD why a fold
#       failed — never swallow into a bare NULL.
# We do NOT recompute start_event on the fold: holdout is by whole animals, so
# each kept animal's contiguity (and thus AR.start) is preserved.
cv_deviance <- function(spec) {
  # FOLD BY REPLICATE, not by animal. On the initial screen each condition is a
  # unique well with ~1 animal per replicate, so animal-stratified CV drops whole
  # CONDITIONS from training ("factor levels not in original fit") — ill-posed.
  # Each condition appears once per replicate, so holding out a whole replicate
  # keeps every condition represented in the remaining replicates. Folds = the
  # distinct `plate` values (e.g. "Plate 1_Replicate 1/2/3"). If CV_FOLD_BY is
  # "animal" (default for balanced data like confirmation), fall back to the
  # original animal-stratified scheme.
  excl <- excl_terms(spec)
  if (identical(CV_FOLD_BY, "replicate")) {
    fold_units <- sort(unique(as.character(gam_df$plate)))
    fold_col   <- as.character(gam_df$plate)
  } else {
    set.seed(CV_SEED)
    animals <- sort(levels(gam_df$animal_id))
    fa <- setNames(sample(rep_len(seq_len(N_CV_FOLDS), length(animals))), animals)
    fold_units <- as.character(seq_len(N_CV_FOLDS))
    fold_col   <- as.character(fa[as.character(gam_df$animal_id)])
  }
  n_folds <- length(fold_units)
  devs <- numeric(n_folds); errs <- character(0)
  for (k in seq_len(n_folds)) {
    ho <- fold_units[k]
    tr <- as.data.frame(gam_df[fold_col != ho])
    te <- gam_df[fold_col == ho]
    f <- build_formula(spec$cond_mode, spec$cond_k, spec$animal_mode, spec$animal_k, base_k = spec$base_k %||% k_prod)
    # bam(discrete=TRUE) captures its call and RE-EVALUATES data / AR.start
    # lazily IN THE FORMULA'S ENVIRONMENT, where this function's local `tr` /
    # AR.start aren't visible ("object 'tr' not found"). Fix: build an env that
    # holds tr + ar + the needed globals, attach it AS the formula's environment,
    # and eval the bam() call there. (Verified: environment(f)<-e is the part
    # that actually resolves it — eval env alone is not enough.)
    fit_env <- new.env(parent = globalenv())
    fit_env$tr <- tr
    fit_env$ar <- tr$start_event
    fit_env$ABL_RHO <- ABL_RHO
    fit_env$usable_cores <- usable_cores
    environment(f) <- fit_env
    fit_env$f <- f
    m <- tryCatch(eval(quote(bam(formula = f, data = tr, family = tw(),
                      rho = ABL_RHO, AR.start = ar, select = TRUE,
                      method = "fREML", discrete = TRUE,
                      nthreads = usable_cores, chunk.size = 10000)),
                  envir = fit_env),
                  error = function(e) { errs[[length(errs)+1]] <<-
                    sprintf("fold %d fit: %s", k, conditionMessage(e)); NULL })
    if (is.null(m)) { devs[k] <- NA_real_; next }

    # (a)+(b): plain data.frame newdata, placeholder train levels for the
    # excluded random smooths.
    cv_df <- as.data.frame(te)
    tr_animals <- droplevels(factor(as.character(tr$animal_id)))
    cv_df$animal_id <- factor(levels(tr_animals)[1], levels = levels(tr_animals))
    if (n_plates >= 2) {
      tr_plates <- droplevels(factor(as.character(tr$plate)))
      cv_df$plate <- factor(levels(tr_plates)[1], levels = levels(tr_plates))
    }
    # (c): staged predict
    mu <- tryCatch(predict(m, newdata = cv_df, type = "response", exclude = excl),
                   error = function(e) { errs[[length(errs)+1]] <<-
                     sprintf("fold %d predict: %s", k, conditionMessage(e)); NULL })
    if (is.null(mu) || all(!is.finite(mu))) {
      if (is.null(mu)) {} else errs[[length(errs)+1]] <- sprintf("fold %d: all preds invalid", k)
      devs[k] <- NA_real_; next
    }
    # (c): staged deviance. mgcv's Tweedie dev.resids signature is
    # (y, mu, wt, theta) — calling it with only 3 positional args leaves `theta`
    # to recycle and triggers "longer object length is not a multiple of
    # shorter" on some mgcv builds. Pass theta explicitly (the fitted power
    # parameter), with a 3-arg fallback for non-Tweedie families.
    y  <- te$pixel_diff
    wt <- rep(1, length(y))
    dr <- tryCatch({
      th <- tryCatch(m$family$getTheta(TRUE), error = function(e) NULL)
      if (!is.null(th) && "theta" %in% names(formals(m$family$dev.resids))) {
        m$family$dev.resids(y, mu, wt, theta = th)
      } else {
        m$family$dev.resids(y, mu, wt)
      }
    }, error = function(e) { errs[[length(errs)+1]] <<-
        sprintf("fold %d dev.resids: %s", k, conditionMessage(e)); NULL })
    devs[k] <- if (is.null(dr)) NA_real_ else sum(dr, na.rm = TRUE)
  }
  if (length(errs)) cat("    CV issues: ", paste(errs, collapse = " | "), "\n", sep = "")
  n_ok <- sum(is.finite(devs))
  # Return the PER-FOLD deviance vector too. The fold assignment is identical
  # across variants (fixed CV_SEED + sorted animals), so per-fold deviances are
  # PAIRED — enabling a paired comparison vs baseline (fold-to-fold variation
  # cancels, much tighter than the between-fold SE).
  list(cv_dev_mean  = if (n_ok > 0) mean(devs, na.rm = TRUE) else NA_real_,
       cv_dev_se    = if (n_ok > 1) sd(devs, na.rm = TRUE) / sqrt(n_ok) else NA_real_,
       cv_dev_folds = devs)
}


# ── STAGE 3 GATE: run CV ONLY on eligible finalists ──────────────────────────
# Stage 1 (one fit/cell) is done above for every cell. Now filter by adequacy,
# then run the expensive 5-fold CV gate only on survivors.
base_wall <- records[[BASELINE_KEY]]$wall_secs %||% NA_real_
base_rss  <- records[[BASELINE_KEY]]$peak_rss_mb %||% NA_real_
elig_keys <- names(records)[vapply(records, eligible, logical(1))]
cv_keys   <- if (ABL_CV) elig_keys else character(0)
cat("\nEligible after Stage-1 filters: ", length(elig_keys), "/",
    length(records), "  -> CV on: ",
    if (length(cv_keys)) paste(cv_keys, collapse = ", ") else "(none)", "\n", sep = "")
cv_results <- list()
for (key in cv_keys) {
  cat("  CV: ", key, " ...\n", sep = "")
  cv_results[[key]] <- tryCatch(cv_deviance(records[[key]]$spec),
                                error = function(e) list(cv_dev_mean = NA_real_,
                                                         cv_dev_se = NA_real_))
}

rows <- lapply(records, function(r) {
  elig <- eligible(r)
  cv <- cv_results[[r$variant]] %||% list(cv_dev_mean = NA_real_, cv_dev_se = NA_real_)
  data.frame(
    variant       = r$variant,
    param_mode    = ABL_PARAM_MODE,
    fit_ok        = isTRUE(r$fit_ok),
    converged     = if (isTRUE(r$fit_ok)) isTRUE(r$converged) else NA,  # key end-of-Phase-A metric
    wall_secs     = r$wall_secs %||% NA_real_,
    speedup_x     = if (is.finite(base_wall) && is.finite(r$wall_secs %||% NA)) base_wall / r$wall_secs else NA_real_,
    peak_rss_mb   = r$peak_rss_mb %||% NA_real_,
    mem_ratio     = if (is.finite(base_rss) && is.finite(r$peak_rss_mb %||% NA)) (r$peak_rss_mb %||% NA) / base_rss else NA_real_,
    n_coef        = r$n_coef %||% NA_integer_,
    dev_explained = r$dev_explained %||% NA_real_,
    aic           = r$aic %||% NA_real_,
    concurvity_max= r$concurvity_max %||% NA_real_,
    k_index_pass  = r$k_index_pass %||% NA_real_,
    fit_corr      = r$fit_corr %||% NA_real_,
    fit_rmse      = r$fit_rmse %||% NA_real_,
    fit_ccc       = r$fit_ccc %||% NA_real_,
    # inferential similarity (same conclusions?)
    n_sig         = r$n_sig %||% NA_integer_,
    est_corr      = r$est_corr %||% NA_real_,
    sig_agree     = r$sig_agree %||% NA_real_,
    sig_recall    = r$sig_recall %||% NA_real_,
    cv_dev_mean   = cv$cv_dev_mean %||% NA_real_,
    cv_dev_se     = cv$cv_dev_se %||% NA_real_,
    eligible      = elig,
    stringsAsFactors = FALSE)
})
summary_df <- bind_rows(rows)
# Sort: eligible first, then by functional similarity (fit_corr desc).
summary_df <- summary_df[order(-summary_df$eligible,
                               -(summary_df$fit_corr %||% 0), na.last = TRUE), ]

cat("\n================ ABLATION SUMMARY ================\n")
print(summary_df, row.names = FALSE, digits = 5)
write_csv(summary_df, file.path(OUT_DIR, "ablation_summary.csv"))

# ── PER-CONDITION CONTRASTS (long) + significance-flip vs baseline ───────────
# Every variant's per-condition contrast vs DMSO, tagged by variant, plus a flag
# for whether its significance CALL differs from baseline's for that condition.
# This is what lets us audit exactly which drugs change call between formulas.
if (length(contr_all)) {
  contr_long <- bind_rows(contr_all)
  if (!is.null(baseline_contr)) {
    bsig <- setNames(baseline_contr$significant, baseline_contr$condition)
    contr_long$base_significant <- bsig[contr_long$condition]
    contr_long$flips_vs_base <- contr_long$significant != contr_long$base_significant
  }
  contr_long <- contr_long[order(contr_long$condition, contr_long$variant), ]
  write_csv(contr_long, file.path(OUT_DIR, "per_condition_contrasts.csv"))
  cat("\nSaved: per_condition_contrasts.csv  (", nrow(contr_long),
      " rows)\n", sep = "")
  # Quick report: which (variant, condition) pairs flip significance vs baseline
  if ("flips_vs_base" %in% names(contr_long)) {
    flips <- contr_long[isTRUE(contr_long$flips_vs_base) |
                        (!is.na(contr_long$flips_vs_base) & contr_long$flips_vs_base), ]
    if (nrow(flips)) {
      cat("\n── Significance flips vs baseline ──\n")
      print(flips[, c("variant", "condition", "estimate", "p_value",
                      "significant", "base_significant")],
            row.names = FALSE, digits = 4)
    } else cat("\n(No significance flips vs baseline.)\n")
  }
}
cat("\nSaved: ", file.path(OUT_DIR, "ablation_summary.csv"), "\n", sep = "")

# ── PAIRED CV COMPARISON vs BASELINE ─────────────────────────────────────────
# The CV fold assignment is identical across variants (fixed seed + sorted
# animals), so per-fold deviances are PAIRED with the baseline. A paired test
# cancels fold-to-fold variation (the dominant noise in the between-fold SE),
# giving a far more sensitive read on whether a variant truly differs from
# baseline out-of-sample. Negative diff = LOWER deviance = BETTER than baseline.
base_folds <- cv_results[[BASELINE_KEY]]$cv_dev_folds
if (!is.null(base_folds) && sum(is.finite(base_folds)) >= 2) {
  paired_rows <- list()
  for (key in names(cv_results)) {
    if (key == BASELINE_KEY) next
    vf <- cv_results[[key]]$cv_dev_folds
    if (is.null(vf)) next
    d  <- vf - base_folds                      # per-fold paired difference
    ok <- is.finite(d)
    if (sum(ok) < 2) next
    md <- mean(d[ok]); sed <- sd(d[ok]) / sqrt(sum(ok))
    paired_rows[[key]] <- data.frame(
      variant       = key,
      mean_diff     = md,                       # vs baseline (neg = better)
      se_diff       = sed,
      t_paired      = if (sed > 0) md / sed else NA_real_,
      pct_change    = 100 * md / mean(base_folds[ok], na.rm = TRUE),
      n_folds       = sum(ok),
      stringsAsFactors = FALSE)
  }
  if (length(paired_rows)) {
    paired_df <- bind_rows(paired_rows)
    paired_df <- paired_df[order(paired_df$mean_diff), ]   # best (most negative) first
    cat("\n──────── PAIRED CV vs baseline (per-fold, same folds) ────────\n")
    cat("mean_diff < 0  => variant generalizes BETTER than baseline.\n")
    cat("|t_paired| >~ 2.5 (df=4) => difference is significant at ~0.05.\n\n")
    print(paired_df, row.names = FALSE, digits = 4)
    write_csv(paired_df, file.path(OUT_DIR, "cv_paired_vs_baseline.csv"))
    cat("\nSaved: cv_paired_vs_baseline.csv\n")
  }
} else {
  cat("\n(Paired CV skipped — baseline per-fold deviances unavailable.)\n")
}

# ── PAIRED CV: pairwise ANIMAL-TERM comparisons at each cond-k ────────────────
# At each cond-k, compare the three animal terms pairwise, per fold (same folds),
# so the animal-term effect is isolated from cond-k. mean_diff = mean(A - B):
# negative => term A generalizes BETTER than term B at that cond-k.
animal_pairs <- list(
  c("re",   "sz10"),   # random intercept vs reduced smooth
  c("re",   "sz30"),   # random intercept vs full smooth
  c("sz10", "sz30")    # reduced vs full smooth
)
avs_rows <- list()
for (ck in c(30, 20, 10, 5)) {
  for (pr in animal_pairs) {
    ka <- sprintf("cond%d__animal%s", ck, pr[1])
    kb <- sprintf("cond%d__animal%s", ck, pr[2])
    fa <- cv_results[[ka]]$cv_dev_folds; fb <- cv_results[[kb]]$cv_dev_folds
    if (is.null(fa) || is.null(fb)) next
    d <- fa - fb; ok <- is.finite(d)
    if (sum(ok) < 2) next
    md <- mean(d[ok]); sed <- sd(d[ok]) / sqrt(sum(ok))
    avs_rows[[paste(ck, pr[1], pr[2])]] <- data.frame(
      cond_k = ck, comparison = sprintf("%s - %s", pr[1], pr[2]),
      mean_diff = md, se_diff = sed,
      t_paired = if (sed > 0) md / sed else NA_real_,
      A_better = md < 0, n_folds = sum(ok), stringsAsFactors = FALSE)
  }
}
if (length(avs_rows)) {
  avs_df <- bind_rows(avs_rows)
  cat("\n──────── PAIRED CV: animal-term comparisons at each cond-k ────────\n")
  cat("mean_diff < 0 => term A (left) generalizes BETTER than term B at that cond-k.\n")
  cat("|t_paired| >~ 2.78 (df=4) => significant at ~0.05.\n\n")
  print(avs_df, row.names = FALSE, digits = 4)
  write_csv(avs_df, file.path(OUT_DIR, "cv_animal_term_pairwise.csv"))
  cat("\nSaved: cv_animal_term_pairwise.csv\n")
}

cat("\nHow to read this:\n",
    "  PRIMARY    = functional similarity (fit_corr/fit_ccc: same trajectories?)\n",
    "             + est_corr (same effect estimates / conclusions?), threshold-free\n",
    "  SPEEDUP    = speedup_x (wall) and mem_ratio (peak RSS) vs baseline\n",
    "  GATE       = converged == TRUE, and CV deviance (cv_dev_mean) among eligible\n",
    "  sig_agree  = decision-level sanity check ONLY; brittle at the p=0.05\n",
    "               boundary, so it understates similarity for borderline\n",
    "               conditions. Trust est_corr over sig_agree.\n",
    "Winner = highest functional + estimate similarity with the best speedup,\n",
    "among converged, CV-passing variants.\n",
    sep = "")
