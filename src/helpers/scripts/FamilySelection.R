# ── FamilySelection.R ──────────────────────────────────────────────────────────
# Optional distribution family selection test for Finomena BAM analysis.
# Compares Tweedie, Gamma(log), and Negative Binomial using the full production
# BAM formula (by-condition smooth, plate factor smooth, per-animal factor
# smooth) fit per phase group. AR1 rho is fixed at 0.25 to halve compute (vs the
# production two-pass rho optimisation); this is the only deviation from the
# downstream model. Family ranking is robust to small rho mis-specification.
#
# Outputs (to output_dir):
#   family_selection_results.csv   — per-family aggregated AIC/BIC/dispersion table
#   family_selection_by_group.csv  — per-group × per-family breakdown
#   family_selection_winner.csv    — winning family name and shift value
#
# Called from the app as:
#   Rscript "FamilySelection.R" <input_csv> <output_dir> <var_names> <ref_values> <ref_condition>
# ──────────────────────────────────────────────────────────────────────────────

# ── Package check ──────────────────────────────────────────────────────────────
cat("R library path:", paste(.libPaths(), collapse = "\n               "), "\n")
required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel")
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

# ── Argument parsing ───────────────────────────────────────────────────────────
args           <- commandArgs(trailingOnly = TRUE)
input_csv      <- args[1]
output_dir     <- if (length(args) >= 2 && !is.na(args[2])) args[2] else dirname(input_csv)
var_names_raw  <- if (length(args) >= 3 && !is.na(args[3]) && nchar(args[3]) > 0) args[3] else "Genotype,Drug"
ref_values_raw <- if (length(args) >= 4 && !is.na(args[4]) && nchar(args[4]) > 0) args[4] else "WT,DMSO"
ref_condition  <- if (length(args) >= 5 && !is.na(args[5]) && nchar(args[5]) > 0) args[5] else "WT+DMSO"

if (is.na(input_csv) || !file.exists(input_csv)) {
  stop("ERROR: input_csv argument missing or file not found: ", input_csv)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

k_start <- 30  # smooth-basis dimension (echoed below alongside the formula)

cat("Input CSV:  ", input_csv,  "\n")
cat("Output dir: ", output_dir, "\n")
cat("k (smooth basis):", k_start, "\n")

# ── Variable parsing ────────────────────────────────────────────────────────────
var_names_raw_vec <- str_split(var_names_raw, ",")[[1]]
ref_values        <- str_split(ref_values_raw, ",")[[1]]
n_vars            <- length(var_names_raw_vec)
var_names         <- make.names(var_names_raw_vec, unique = TRUE)
ref_map           <- setNames(ref_values[1:n_vars], var_names)

cat("Variables:  ", paste(var_names, collapse = ", "), "\n")
cat("References: ", paste(ref_map, collapse = ", "), "\n\n")

# ── Data loading ───────────────────────────────────────────────────────────────
full_df <- read_csv(input_csv, show_col_types = FALSE) %>%
  rename(pxl_diff = pixel_diff) %>%
  mutate(
    Condition = as.character(Condition),
    Phase     = as.character(Phase),
    Group     = as.integer(Group),
    time_sec  = as.integer(time_sec),
    pxl_diff  = as.numeric(pxl_diff),
    plate     = as.character(plate),
    location  = as.character(location)
  ) %>%
  filter(!is.na(Condition), Condition != "")

# ── Data prep ──────────────────────────────────────────────────────────────────
# Mirrors the unified-model setup in TweedieAR1_BAM.R: one fit across all phase
# groups, with time_in_group reset per phase, and start_event resetting AR(1)
# at both well starts AND phase boundaries within wells.
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
  mutate(start_event = (row_number() == 1) | (Group != lag(Group, default = -1L))) %>%
  ungroup()

if (!all(var_names_raw_vec %in% colnames(gam_df))) {
  gam_df <- gam_df %>%
    separate(Condition, into = var_names, sep = "\\+", remove = FALSE)
} else {
  for (i in seq_along(var_names)) {
    if (var_names_raw_vec[i] != var_names[i] && var_names_raw_vec[i] %in% colnames(gam_df)) {
      gam_df <- gam_df %>% rename(!!var_names[i] := !!var_names_raw_vec[i])
    }
  }
}

for (v in var_names) {
  gam_df[[v]] <- relevel(as.factor(gam_df[[v]]), ref = ref_map[[v]])
}
gam_df$Condition_Combo <- relevel(as.factor(gam_df$Condition_Combo), ref = ref_condition)

unique_groups <- sort(unique(gam_df$Group))
n_groups      <- length(unique_groups)

# ── Observed zero proportion ────────────────────────────────────────────────────
zero_prop_obs <- mean(gam_df$pxl_diff == 0, na.rm = TRUE)
cat(sprintf("Observed zero proportion: %.2f%%\n", zero_prop_obs * 100))

# Shift so Gamma can run: shift minimum value up to 1 if any zeros/negatives
min_val   <- min(gam_df$pxl_diff, na.rm = TRUE)
shift_val <- if (min_val <= 0) (1 - min_val) else 0
cat(sprintf("Shift for Gamma family:   %.4f\n\n", shift_val))

# ── Full production-style unified formula ──────────────────────────────────────
# Mirrors the BAM script: ONE fit across all phase groups, with phase encoded
# via Group_factor in the parametric expansion + by-interaction smooth, and
# random-effect smooths over global time_sec. Rho is fixed (not optimised)
# for speed.
n_plates <- length(unique(gam_df$plate))
if (n_plates >= 2) {
  plate_smooth_str <- " + s(time_sec, plate, bs = 'fs', m = 1)"
} else {
  plate_smooth_str <- ""
  cat("Only ", n_plates, " plate(s) detected — omitting the plate factor smooth.\n",
      sep = "")
}

fixed_rho      <- 0.20
parametric_str <- paste(c(var_names, "Group_factor"), collapse = " * ")
test_formula   <- as.formula(paste0(
  "pxl_diff_test ~ ", parametric_str,
  " + s(time_in_group, k = k_start)",
  " + s(time_in_group, by = interaction(Condition_Combo, Group_factor), k = k_start)",
  plate_smooth_str,
  " + s(time_sec, animal_id, bs = 'fs', m = 1)"
))
cat("Test formula:", deparse(test_formula, width.cutoff = 200), "\n")
cat(sprintf("Fixed AR1 rho: %.2f (production optimises this in the unified fit)\n\n", fixed_rho))

usable_cores <- max(1, detectCores() - 1)

# ── Candidate families ──────────────────────────────────────────────────────────
candidates <- list(
  Tweedie     = list(family = tw(),                shift = 0),
  Gamma       = list(family = Gamma(link = "log"), shift = shift_val),
  NegBinomial = list(family = nb(),                shift = 0)
)

# ── One fit per family (no per-group loop in unified design) ───────────────────
total_n      <- nrow(gam_df)
fit_failures <- list()
family_rows  <- list()

for (fname in names(candidates)) {
  cat("========================================\n")
  cat("Family:", fname, "\n")
  cand <- candidates[[fname]]
  gam_df$pxl_diff_test <- gam_df$pxl_diff + cand$shift

  tryCatch({
    m <- bam(
      test_formula,
      data     = gam_df,
      family   = cand$family,
      rho      = fixed_rho,
      AR.start = gam_df$start_event,
      select   = TRUE,
      method   = "fREML",
      discrete = TRUE,
      nthreads = usable_cores
    )

    aic_val  <- AIC(m)
    bic_val  <- BIC(m)
    p_resid  <- residuals(m, type = "pearson")
    disp     <- sum(p_resid^2, na.rm = TRUE) / m$df.residual
    dev_expl <- summary(m)$dev.expl * 100

    pred_zero <- switch(fname,
      Tweedie = {
        p_hat  <- m$family$getTheta(TRUE)
        phi    <- m$sig2
        mu_hat <- fitted(m)
        mean(exp(-mu_hat^(2 - p_hat) / (phi * (2 - p_hat))), na.rm = TRUE)
      },
      NegBinomial = {
        theta  <- m$family$getTheta(TRUE)
        mu_hat <- fitted(m)
        mean((theta / (theta + mu_hat))^theta, na.rm = TRUE)
      },
      Gamma = NA_real_
    )

    cat(sprintf("  AIC=%.1f  BIC=%.1f  disp=%.3f  dev=%.1f%%\n",
                aic_val, bic_val, disp, dev_expl))

    family_rows[[fname]] <- data.frame(
      Family        = fname,
      AIC           = aic_val,
      BIC           = bic_val,
      Dispersion    = disp,
      Dev_Explained = dev_expl,
      Pred_Zero     = pred_zero,
      Obs_Zero      = zero_prop_obs,
      n_obs         = total_n,
      Shift_Applied = cand$shift,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    cat("  FAILED:", msg, "\n")
    fit_failures[[fname]] <<- msg
    family_rows[[fname]] <<- data.frame(
      Family        = fname,
      AIC           = NA_real_,
      BIC           = NA_real_,
      Dispersion    = NA_real_,
      Dev_Explained = NA_real_,
      Pred_Zero     = NA_real_,
      Obs_Zero      = zero_prop_obs,
      n_obs         = total_n,
      Shift_Applied = cand$shift,
      stringsAsFactors = FALSE
    )
  })
}

per_family_df <- bind_rows(family_rows)
# Keep the "by_group" filename for backwards compatibility with the Python
# widget that reads it; the unified fit produces one row per family here.
write_csv(per_family_df, file.path(output_dir, "family_selection_by_group.csv"))

# ── Format final results table (one row per family) ────────────────────────────
results_df <- per_family_df %>%
  mutate(
    Zero_Match = case_when(
      Family == "Gamma"          ~ sprintf("N/A (shift of %.4f applied; Gamma requires >0)", Shift_Applied),
      is.na(Pred_Zero)           ~ "FAILED — see log",
      TRUE                       ~ sprintf("Obs: %.1f%% | Pred: %.1f%%",
                                           zero_prop_obs * 100, Pred_Zero * 100)
    ),
    AIC           = ifelse(is.na(AIC), Inf, round(AIC, 1)),
    BIC           = ifelse(is.na(BIC), Inf, round(BIC, 1)),
    Dispersion    = round(Dispersion, 3),
    Dev_Explained = round(Dev_Explained, 1)
  ) %>%
  select(Family, AIC, BIC, Dispersion, Dev_Explained, Zero_Match, Shift_Applied)

# ── Select winner ────────────────────────────────────────────────────────────────
valid <- results_df %>% filter(is.finite(AIC))

if (nrow(valid) == 0) {
  stop("All family tests failed. Check the log above for details.")
}

winner_row   <- valid %>% slice_min(AIC, n = 1)
winner_name  <- winner_row$Family
winner_shift <- winner_row$Shift_Applied

cat("\n========================================\n")
cat("FAMILY SELECTION RESULTS:\n")
print(results_df, row.names = FALSE)
cat(sprintf("\nRECOMMENDED FAMILY: %s  (lowest AIC = %.1f)\n", winner_name, winner_row$AIC))
cat(sprintf("SHIFT VALUE:        %.4f\n", winner_shift))

# ── Save outputs ─────────────────────────────────────────────────────────────────
write_csv(results_df, file.path(output_dir, "family_selection_results.csv"))
write_csv(
  data.frame(
    family         = winner_name,
    shift          = winner_shift,
    zero_prop_obs  = round(zero_prop_obs * 100, 2)
  ),
  file.path(output_dir, "family_selection_winner.csv")
)

cat("\nFamily selection complete. Results saved to:", output_dir, "\n")
