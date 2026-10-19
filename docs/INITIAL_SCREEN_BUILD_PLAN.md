# Initial-Screen Production Run — Build Plan (for review)

## Context / why this exists

The initial screen (3 plates, ~82 drugs/plate, 658 animals) cannot use the
confirmation production formula: the per-animal smooth is computationally
intractable at scale, Genotype×Drug is non-estimable (no WT+drug cells), and
combining plates confounds drug effect with batch (each drug lives on one
plate). After an extensive ablation we locked an initial-screen formula and a
per-plate analysis design. This plan turns that into a production script.

Work happens in **`initial_screen_script.r`** (a copy of `stand-alone_script.r`).
The confirmation script `stand-alone_script.r` is NOT touched.

## Locked decisions (from the ablation + design discussion)

- **Formula (per plate):**
  `pixel_diff ~ Condition_Combo + s(time_in_group, k=30) + s(time_in_group, by=Condition_Combo, k=30) + s(animal_id, bs="re")`
  - parametric = `Condition_Combo` (single factor; Genotype×Drug non-estimable here)
  - shared baseline smooth k=30 (verified adequate: contrasts invariant above k=30)
  - by-condition smooth k=30 (FDR hit list identical k=30→50; D08 edge case didn't survive correction)
  - animal = random intercept `bs="re"` (sz per-animal smooth was intractable at scale)
- **AR(1):** two-pass estimated rho (production-correct), per plate.
- **Per-plate fitting:** each plate fit separately (own rho, own DMSO reference,
  own error scale) — required by the batch confound. ALL replicates per plate.
- **One job looping all 3 plates**, per-plate output folders, clear per-plate
  log markers.
- **Multiple-testing correction:** FDR (BH), applied WITHIN each family separately
  (not pooled across families).

## Experimental design → contrast structure

Conditions present: `WT+DMSO` (negative control), `MUT+DMSO` (= KO+DMSO, positive
control), `MUT+<drug>` for each drug. No `WT+drug`. This is a rescue screen.

Three contrast families, each scored TWO ways (difference test + equivalence test):

1. **Genotype (QC, computed once):** `MUT+DMSO` vs `WT+DMSO` — the disease window.
   Difference: Wald p, raw (single test, no FDR). Equivalence: Pr(M<δ).
2. **Drug-effect family (per drug):** `MUT+drug` vs `MUT+DMSO` — did the drug move
   the mutant? Difference: Wald p → FDR within the ~82. Equivalence: Pr(M<δ).
3. **Rescue family (per drug):** `MUT+drug` vs `WT+DMSO` — did the drug reach
   wildtype? Difference: Wald p → FDR within the ~82. Equivalence: Pr(M<δ).

FDR is applied separately to family 2 and family 3 (each is "many drugs vs one
control" — the correct multiplicity axis; pooling them would over-correct the
primary rescue endpoint).

**Rescue hit logic (read from the joint pattern, computed downstream / by user):**
rescue = drug-effect difference significant (FDR) AND high Pr(equiv to WT+DMSO).
i.e. moved away from KO+DMSO AND landed within δ of WT+DMSO. The two families are
two coordinates classifying each drug, not two shots at one hit.

## What changes vs the copied confirmation script

### 1. Formula section (small)
- parametric_str → `Condition_Combo` (initial mode)
- animal smooth → `s(animal_id, bs="re")`
- baseline + by-condition k stay 30; single-phase-group path already drops
  Group_factor and plate_smooth when fitting one plate.

### 2. CONTRASTS section (full replacement — the big change)
The current section uses emmeans over `var_names` (Genotype, Drug) for
"per-variable effects" + a Full_Interaction. This is non-estimable under the
Condition_Combo parametrization and must be replaced with the validated
lpmatrix + vcov pairwise-contrast approach (same machinery the ablation's
`contrast_vs_ref` used). New contrasts produce, per family, for each pair:
estimate (time-averaged Δη/log2), SE, z, raw_p. Then `apply_correction()` (the
existing helper, FDR) is run within each family.

### 3. POSTERIOR EQUIVALENCE (generalize, keep the method)
Reuse the existing `.posterior_equivalence_pair` (Bayesian posterior draws;
M = max-over-time |Δη(t)|/log2; Pr(M<δ)) UNCHANGED in method, but:
- run it on BOTH the drug-effect pairs (vs MUT+DMSO) and rescue pairs (vs WT+DMSO),
  plus the genotype pair — i.e. drive it from the same per-family pair lists as
  the difference contrasts (one consistent reference per family).
- the `parse_cond`/`build_nd` internals assume `var_names` parts; under
  Condition_Combo parametrization, build_nd must set Condition_Combo (and the
  split vars for any parametric use) consistently — verify the newdata matches
  what the initial formula needs (Condition_Combo is the parametric term now).
- exclude string: animal term is now `s(animal_id)` (re), not
  `s(time_in_group,animal_id)`. Must update the `excl` accordingly.

### 4. δ (equivalence margin)
- Keep δ = 1.0 (the original/confirmation value) BUT make it env-overridable
  (e.g. FINOMENA_EQUIV_DELTA) rather than hardcoded.
- Calibration deferred to after the first run: the genotype pair's M (peak
  separation of MUT+DMSO vs WT+DMSO) appears directly in
  posterior_equivalence_summary.csv. Verify the positive control is NOT called
  equivalent to WT (Pr_equiv should be low); re-set δ from real output if needed.
  (Note: earlier 0.91 log2 was the TIME-AVERAGED genotype contrast; the
  equivalence test uses MAX-over-time M, which is larger — so δ=1.0 is plausibly
  fine but unverified until the run.)

### 5. Per-plate loop + outputs
- Wrap the analysis (model fit → contrasts → correction → posterior eq → writes)
  in a loop over plates. Each plate: subset to that plate's rows (all replicates),
  refit independently, write to its own folder, e.g. `output/Plate 1/`, `Plate 2/`,
  `Plate 3/`.
- Clear per-plate log markers (e.g. "===== PLATE 2/3: Plate 2 =====") so progress
  is legible in CloudWatch.
- The script is currently linear (top-to-bottom, one fit). The loop is a real
  refactor: data-load/aggregation can stay shared (or per-plate), but the
  model+contrast+correction+posterior block must become per-plate. NOTE: this
  needs care — `stand-alone_script.r` was written as a single pass; verify no
  global state leaks between plate iterations (factor levels, rho, gam_df, the
  contrast accumulators, threading set-once blocks).

## Required outputs (per plate, in its own folder) — all needed

- `master_results.csv` — all three families' contrasts: estimate, SE, raw_p,
  adjusted_p (FDR within family), plus M_median, M_95CI, Pr_equiv from the
  equivalence test. Tagged by Test_Family (Genotype / Drug_effect / Rescue).
- `summary_statistics_by_group.csv` — with 1 phase group this collapses; repurpose
  to per-family hit counts (n significant by difference-FDR; n equivalent by
  Pr(M<δ) threshold). CONFIRM desired contents.
- `posterior_equivalence_summary.csv` and `posterior_equivalence_draws.csv`
- `aggregated_data.csv`
- `no_movement_animals.csv`

## Open items — RESOLVED

1. **summary_statistics_by_group.csv** → per-family hit counts. Reuse the
   original script's summary logic (`group_by(Test_Family)`, `n_total`,
   `n_sig = sum(adjusted_pvalue < CORRECTION_SPEC$threshold)`, `mean_effect`),
   dropping the per-Group pivot (1 phase group). Computed per plate.
2. **Equivalence threshold** → NONE. The original script only *reports*
   `Pr_equiv` (= mean(M < δ)) and M median/95% CI; it never hard-thresholds
   equivalence. The initial-screen script does the same — reports Pr_equiv per
   pair, no baked-in cutoff. The rescue classification (drug-effect FDR-sig AND
   high Pr_equiv-to-WT) is done downstream by the user at interpretation time.
   Difference-test hits DO use the FDR threshold (adjusted_pvalue < 0.05) per
   family, exactly as the original summary does.
3. **Aggregation** → runs ONCE (writes aggregated_data.csv), then the per-plate
   loop subsets that frame by plate. No redundant .xls re-parsing.

## Verification

1. Parse-check `initial_screen_script.r` with the installed Rscript.
2. Dry tractability: run Plate 1 only first (env/flag) — confirm it fits with
   2-pass rho, produces all output files, and the genotype M is sane.
3. Sanity: genotype contrast (MUT+DMSO vs WT+DMSO) significant by difference;
   positive control NOT equivalent to WT (low Pr_equiv) — else δ too loose.
4. Then run all 3 plates in one looped job.

## Out of scope
- The confirmation screen / `stand-alone_script.r` (untouched).
- Cross-plate hit comparison (batch confound — per-plate only).
