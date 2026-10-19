# ============================================================================
# corrections.R — Multiple-testing correction algorithms for Finomena
# ============================================================================
#
# All functions take a tagged data frame of raw p-values (`raw_pvalue`) plus
# any tree metadata they need, and return the same data frame augmented with:
#
#   adjusted_pvalue    numeric, one per row
#   correction_method  string, identical across all rows
#   tree_level         integer or NA  (NA for flat methods)
#   tree_status        string or NA   ("passed" / "gate_failed" / "not_reached"
#                                       for tree; NA for flat)
#
# Sourced by TweedieAR1_BAM.R after the master table is assembled.
#
# Inputs
# ------
# `master_df`  Required columns: raw_pvalue, Test_Family, Group, Split_By,
#              Passed_Contrasts. For tree methods, also a `tree_level` column
#              assigning each row to a level and a `tree_parent_id` column
#              identifying its parent in the tree (NA at the root).
# `q`         Threshold (q for FDR, alpha for FWER).
# `tree_spec` (tree methods only) List describing the tree structure — see
#              the assign_tree_metadata() helper.
#
# Returns the augmented master_df.
# ============================================================================


# ---------------------------------------------------------------------------
# Flat correction: BH or Holm within Test_Family.
# ---------------------------------------------------------------------------
#
# Always corrects within each Test_Family because each family represents a
# distinct scientific question (Genotype_Effect / Drug_Effect /
# Full_Interaction etc.). The historic "pooled across all" scope was removed
# as it conflated unrelated hypotheses into a single overly conservative
# correction. The `scope` arg is accepted but ignored for back-compat with
# any existing correction.json files; per-family is enforced.
# ---------------------------------------------------------------------------
.flat_correct <- function(master_df, method, q, scope = "per_family") {
  method <- match.arg(method, c("BH", "holm"))
  if (!identical(scope, "per_family")) {
    warning("'pooled' scope is no longer supported — using per_family.")
  }
  master_df <- master_df %>%
    group_by(Test_Family) %>%
    mutate(adjusted_pvalue = p.adjust(raw_pvalue, method = method)) %>%
    ungroup()
  master_df$correction_method <- if (method == "BH") "BH" else "Holm"
  master_df$tree_level        <- NA_integer_
  master_df$tree_status       <- NA_character_
  master_df
}


# ---------------------------------------------------------------------------
# TreeBH — Yekutieli (2008) hierarchical FDR control on a rooted tree.
# ---------------------------------------------------------------------------
#
# Algorithm:
#   Step 1.   Apply BH at level q to the root's children (the L1 family).
#   Step >1.  For each rejected parent at level L, apply BH at level q to its
#             children at level L+1.
#
# A child whose parent was NOT rejected is marked tree_status = "not_reached"
# (its raw p-value isn't tested; adjusted_pvalue stays NA).
#
# `tree_df`: must contain tree_level (integer >= 1) and tree_parent_id (string
# matching another row's tree_id, or NA at level 1) columns. tree_id is a
# unique label per row used to link parent→child.
# ---------------------------------------------------------------------------
.treebh_correct <- function(tree_df, q) {
  tree_df$adjusted_pvalue <- NA_real_
  # Off-tree rows (NA tree_level): the tree spec didn't capture them, so
  # tree-FDR has no opinion. Tag separately from "not_reached" (which means
  # "in the tree but upstream parent didn't reject").
  tree_df$tree_status <- ifelse(is.na(tree_df$tree_level),
                                "off_tree", "not_reached")

  levels_sorted <- sort(unique(tree_df$tree_level))
  if (length(levels_sorted) == 0) return(tree_df)

  # ── L1: BH across all root-level rows ────────────────────────────────────
  l1 <- which(tree_df$tree_level == levels_sorted[1])
  if (length(l1) == 0) return(tree_df)

  l1_p <- tree_df$raw_pvalue[l1]
  l1_adj <- p.adjust(l1_p, method = "BH")
  tree_df$adjusted_pvalue[l1] <- l1_adj
  tree_df$tree_status[l1] <- ifelse(l1_adj < q, "passed", "gate_failed")

  # ── Subsequent levels: per-parent BH on children of rejected parents ────
  for (L in levels_sorted[-1]) {
    rows_at_L <- which(tree_df$tree_level == L)
    if (length(rows_at_L) == 0) next

    for (parent_id in unique(tree_df$tree_parent_id[rows_at_L])) {
      if (is.na(parent_id)) next
      parent_row <- which(tree_df$tree_id == parent_id)
      if (length(parent_row) == 0) next
      parent_status   <- tree_df$tree_status[parent_row]
      parent_rejected <- isTRUE(parent_status == "passed") &&
                         isTRUE(tree_df$adjusted_pvalue[parent_row] < q)
      if (!parent_rejected) {
        # Children of a non-rejected parent stay "not_reached".
        next
      }
      child_rows <- rows_at_L[tree_df$tree_parent_id[rows_at_L] == parent_id]
      child_p    <- tree_df$raw_pvalue[child_rows]
      child_adj  <- p.adjust(child_p, method = "BH")
      tree_df$adjusted_pvalue[child_rows] <- child_adj
      tree_df$tree_status[child_rows] <-
        ifelse(child_adj < q, "passed", "gate_failed")
    }
  }
  tree_df$correction_method <- "TreeBH"
  tree_df
}


# ---------------------------------------------------------------------------
# Holm-gatekeeping — fixed-sequence Holm on a rooted tree (FWER analogue of
# TreeBH).
# ---------------------------------------------------------------------------
#
# Identical traversal pattern to TreeBH but applies p.adjust(method="holm") at
# each within-family step. Rejected = adjusted_pvalue < alpha. Children of a
# non-rejected parent stay "not_reached".
# ---------------------------------------------------------------------------
.holm_gatekeep_correct <- function(tree_df, alpha) {
  tree_df$adjusted_pvalue <- NA_real_
  tree_df$tree_status <- ifelse(is.na(tree_df$tree_level),
                                "off_tree", "not_reached")

  levels_sorted <- sort(unique(tree_df$tree_level))
  if (length(levels_sorted) == 0) return(tree_df)

  l1 <- which(tree_df$tree_level == levels_sorted[1])
  if (length(l1) == 0) return(tree_df)

  l1_adj <- p.adjust(tree_df$raw_pvalue[l1], method = "holm")
  tree_df$adjusted_pvalue[l1] <- l1_adj
  tree_df$tree_status[l1] <- ifelse(l1_adj < alpha, "passed", "gate_failed")

  for (L in levels_sorted[-1]) {
    rows_at_L <- which(tree_df$tree_level == L)
    if (length(rows_at_L) == 0) next

    for (parent_id in unique(tree_df$tree_parent_id[rows_at_L])) {
      if (is.na(parent_id)) next
      parent_row <- which(tree_df$tree_id == parent_id)
      if (length(parent_row) == 0) next
      parent_status   <- tree_df$tree_status[parent_row]
      parent_rejected <- isTRUE(parent_status == "passed") &&
                         isTRUE(tree_df$adjusted_pvalue[parent_row] < alpha)
      if (!parent_rejected) next

      child_rows <- rows_at_L[tree_df$tree_parent_id[rows_at_L] == parent_id]
      child_adj  <- p.adjust(tree_df$raw_pvalue[child_rows], method = "holm")
      tree_df$adjusted_pvalue[child_rows] <- child_adj
      tree_df$tree_status[child_rows] <-
        ifelse(child_adj < alpha, "passed", "gate_failed")
    }
  }
  tree_df$correction_method <- "Holm-gatekeeping"
  tree_df
}


# ---------------------------------------------------------------------------
# graphicalMCP correlation-aware FWER procedure on the tree.
# ---------------------------------------------------------------------------
#
# Built on the `graphicalMCP` CRAN package (pure-R, no Java unlike gMCP).
# Uses `graph_test_closure(test_types = "parametric", test_corr = ...)` to get
# correlation-aware critical values for ref-anchored Dunnett-style contrasts.
#
# Tree → graph:
#   * Nodes:    one per row of tree_df
#   * Initial α-weights: all mass at L1, split uniformly across L1 siblings
#   * Edges:    parent → its children with uniform weight 1/(n_children);
#               within-level siblings cross-recycle (Holm-style) so a
#               rejection in one sibling propagates to peers before falling
#               down to the children
#   * test_groups: one group per tree level so closure-principle testing
#                   respects the hierarchy
#
# `correlation_matrix`: full joint correlation matrix of contrast estimates
# (derived from emmeans's vcov()). Same row order as tree_df. If NULL or
# wrong size, falls back to Bonferroni-style within-group testing.
#
# Falls back to Holm-gatekeeping if graphicalMCP isn't installed or fails.
# ---------------------------------------------------------------------------
.gmcp_correct <- function(tree_df, alpha, correlation_matrix = NULL) {
  if (!requireNamespace("graphicalMCP", quietly = TRUE)) {
    warning("graphicalMCP package not installed — falling back to Holm-gatekeeping.")
    out <- .holm_gatekeep_correct(tree_df, alpha)
    out$correction_method <- "Holm-gatekeeping (graphicalMCP unavailable)"
    return(out)
  }

  # Pre-set defaults for the WHOLE table. Rows that aren't part of any tree
  # level get tree_status = "off_tree" and stay out of the graphicalMCP call —
  # graphicalMCP enforces "each hypothesis in exactly one group", so off-tree
  # rows MUST be excluded from the call rather than passed with no group.
  tree_df$adjusted_pvalue <- NA_real_
  tree_df$tree_status <- ifelse(is.na(tree_df$tree_level),
                                "off_tree", "not_reached")

  on_tree_idx <- which(!is.na(tree_df$tree_level))
  if (length(on_tree_idx) == 0) {
    warning("No rows assigned to any tree level — nothing to test. ",
            "Check that the correction.json tree spec matches the contrast ",
            "Test_Family / Split_By values in master_results.")
    tree_df$correction_method <- "graphicalMCP (no on-tree rows)"
    return(tree_df)
  }

  # ── Work on a sub-table containing ONLY on-tree rows ─────────────────────
  sub_df <- tree_df[on_tree_idx, , drop = FALSE]
  n_sub  <- nrow(sub_df)
  levels_sorted <- sort(unique(sub_df$tree_level))

  # Initial α-weights: all mass at L1, uniform across siblings
  l1_sub <- which(sub_df$tree_level == levels_sorted[1])
  weights <- numeric(n_sub)
  if (length(l1_sub) > 0) weights[l1_sub] <- 1 / length(l1_sub)

  # Transition matrix: parent → children (uniform) + Holm-style sibling
  # recycling within a level so power isn't wasted when one sibling rejects.
  m <- matrix(0, n_sub, n_sub)
  for (i in seq_len(n_sub)) {
    children <- which(sub_df$tree_parent_id == sub_df$tree_id[i])
    if (length(children) > 0) {
      m[i, children] <- 1 / length(children)
    }
  }
  for (L in levels_sorted) {
    parents <- unique(sub_df$tree_parent_id[sub_df$tree_level == L])
    for (parent_id in parents) {
      if (is.na(parent_id)) next
      siblings <- which(sub_df$tree_level == L &
                        sub_df$tree_parent_id == parent_id)
      if (length(siblings) > 1) {
        for (i in siblings) {
          others <- setdiff(siblings, i)
          has_children <- any(sub_df$tree_parent_id == sub_df$tree_id[i],
                              na.rm = TRUE)
          frac_siblings <- if (has_children) 0.5 else 1.0
          frac_children <- 1 - frac_siblings
          m[i, ] <- m[i, ] * frac_children
          m[i, others] <- frac_siblings / length(others)
        }
      }
    }
  }

  g <- tryCatch(
    graphicalMCP::graph_create(hypotheses = weights, transitions = m),
    error = function(e) NULL
  )
  if (is.null(g)) {
    warning("graphicalMCP::graph_create failed — falling back to Holm-gatekeeping.")
    out <- .holm_gatekeep_correct(tree_df, alpha)
    out$correction_method <- "Holm-gatekeeping (graph_create failed)"
    return(out)
  }

  # test_groups: one group per tree level (relative to sub_df). graphicalMCP
  # enforces that every node belongs to exactly one group — every sub_df row
  # has a level, so the partition is complete here.
  test_groups <- lapply(levels_sorted, function(L) which(sub_df$tree_level == L))

  # Build the optional per-level correlation submatrices.
  use_parametric <- !is.null(correlation_matrix) &&
                    is.matrix(correlation_matrix) &&
                    nrow(correlation_matrix) >= max(on_tree_idx) &&
                    ncol(correlation_matrix) >= max(on_tree_idx)
  if (use_parametric) {
    # Translate sub_df group indices into rows of the full correlation matrix
    # via on_tree_idx, then extract the submatrix per group.
    test_types <- rep("parametric", length(test_groups))
    test_corr  <- lapply(test_groups, function(idx_sub) {
      full_idx <- on_tree_idx[idx_sub]
      sub_corr <- correlation_matrix[full_idx, full_idx, drop = FALSE]
      if (anyNA(sub_corr) || any(abs(sub_corr) > 1 + 1e-8)) NA else sub_corr
    })
  } else {
    test_types <- rep("bonferroni", length(test_groups))
    test_corr  <- rep(list(NA), length(test_groups))
  }

  result <- tryCatch(
    graphicalMCP::graph_test_closure(
      graph       = g,
      p           = sub_df$raw_pvalue,
      alpha       = alpha,
      test_groups = test_groups,
      test_types  = test_types,
      test_corr   = test_corr
    ),
    error = function(e) {
      warning("graphicalMCP::graph_test_closure failed (",
              conditionMessage(e),
              ") — falling back to Holm-gatekeeping.")
      NULL
    }
  )

  if (is.null(result)) {
    out <- .holm_gatekeep_correct(tree_df, alpha)
    out$correction_method <- "Holm-gatekeeping (graphicalMCP failed)"
    return(out)
  }

  # graphicalMCP returns a logical vector of rejections, length n_sub.
  rejected_sub <- as.logical(result$rejected)
  if (length(rejected_sub) != n_sub) {
    rejected_sub <- as.logical(result$outputs$rejected %||% rep(FALSE, n_sub))
  }

  # Map rejections back into the full tree_df row positions, walking the
  # tree level-by-level so we can distinguish tested vs not_reached based on
  # parent rejection status.
  for (L in levels_sorted) {
    rows_at_L_sub <- which(sub_df$tree_level == L)
    if (length(rows_at_L_sub) == 0) next
    full_rows <- on_tree_idx[rows_at_L_sub]
    if (L == levels_sorted[1]) {
      tree_df$tree_status[full_rows] <-
        ifelse(rejected_sub[rows_at_L_sub], "passed", "gate_failed")
      tree_df$adjusted_pvalue[full_rows] <-
        ifelse(rejected_sub[rows_at_L_sub],
               tree_df$raw_pvalue[full_rows],
               pmax(tree_df$raw_pvalue[full_rows], alpha))
    } else {
      for (k in seq_along(rows_at_L_sub)) {
        i_sub  <- rows_at_L_sub[k]
        i_full <- full_rows[k]
        parent_id <- sub_df$tree_parent_id[i_sub]
        if (is.na(parent_id)) next
        parent_sub_idx <- which(sub_df$tree_id == parent_id)
        if (length(parent_sub_idx) == 0) next
        if (!isTRUE(rejected_sub[parent_sub_idx])) {
          tree_df$tree_status[i_full] <- "not_reached"
          next
        }
        tree_df$tree_status[i_full] <-
          if (rejected_sub[i_sub]) "passed" else "gate_failed"
        tree_df$adjusted_pvalue[i_full] <-
          if (rejected_sub[i_sub]) tree_df$raw_pvalue[i_full]
          else max(tree_df$raw_pvalue[i_full], alpha)
      }
    }
  }
  tree_df$correction_method <- "graphicalMCP"
  tree_df
}


# ---------------------------------------------------------------------------
# Dispatcher — picks the correct algorithm based on the correction spec.
# ---------------------------------------------------------------------------
#
# `correction_spec` shape (read from correction.json by the caller):
#   list(
#     strategy      = "flat" | "tree",
#     error_rate    = "FDR" | "FWER",
#     contrast_set  = "all_pairs" | "ref_only",   # affects emmeans, not here
#     threshold     = 0.05,
#     scope         = "per_family",                # flat only (always per_family;
#                                                   # pooled was removed)
#     tree          = list(...)                   # tree only (level specs)
#   )
#
# For tree strategies the caller must have ALREADY annotated master_df with
# tree_id, tree_level, tree_parent_id columns. Helper `assign_tree_metadata()`
# (to come in Phase A3) does that from the tree spec.
# ---------------------------------------------------------------------------
apply_correction <- function(master_df, correction_spec,
                             correlation_matrix = NULL) {
  strat  <- correction_spec$strategy   %||% "flat"
  err    <- correction_spec$error_rate %||% "FDR"
  thr    <- as.numeric(correction_spec$threshold %||% 0.05)
  scope  <- correction_spec$scope %||% "per_family"

  if (strat == "flat") {
    method <- if (err == "FDR") "BH" else "holm"
    return(.flat_correct(master_df, method = method, q = thr, scope = scope))
  }

  # tree
  if (err == "FDR") {
    return(.treebh_correct(master_df, q = thr))
  }
  # FWER tree:
  cs <- correction_spec$contrast_set %||% "all_pairs"
  if (cs == "ref_only") {
    return(.gmcp_correct(master_df, alpha = thr,
                         correlation_matrix = correlation_matrix))
  }
  return(.holm_gatekeep_correct(master_df, alpha = thr))
}
