# ============================================================================
# tree_metadata.R — Assigning tree_id / tree_level / tree_parent_id to rows
# ============================================================================
#
# Tree-based correction algorithms need each contrast row in master_df to
# carry three columns:
#
#   tree_id          unique label per row (used to identify parent→child)
#   tree_level       integer level number (1 = root level)
#   tree_parent_id   tree_id of the row's parent at level (level - 1), or NA
#                    if this is a root-level row
#
# This helper takes a `tree_spec` list (one entry per level, in order) and
# the assembled master_df, and emits these columns. The `tree_spec` shape is:
#
#   tree_spec = list(
#     list(
#       name          = "L1 phase gate",
#       Test_Family   = "Genotype_Effect",
#       split_by      = "DMSO",                   # optional Split_By filter
#       linkage_key   = NULL                      # L1 has no parent
#     ),
#     list(
#       name          = "L2 rescue",
#       Test_Family   = "Full_Interaction",
#       split_by      = NULL,
#       linkage_key   = list(parent_var = "Group",  child_var = "Group")
#     ),
#     list(
#       name          = "L3 specificity",
#       Test_Family   = "Drug_Effect",
#       split_by      = NULL,                     # both WT and KO contexts
#       linkage_key   = list(parent_var = "Drug",   child_var = "Drug")
#     )
#   )
#
# linkage_key tells us which variable's identity threads from parent to child.
# For L1→L2, the obvious linkage is `Group` (phase): a phase gate is the
# parent of the rescue claims in THAT phase. For L2→L3, it's `Drug` (the drug
# being tested): a confirmed rescue at L2 is the parent of its specificity
# tests at L3.
#
# Rows in master_df that don't match any level's filter are returned with
# tree_level / tree_id / tree_parent_id = NA and a tree_status = "off_tree".
# Correction algorithms ignore those rows (they're computed but not gated).
# ============================================================================


# ---------------------------------------------------------------------------
# Apply a single level's filter to master_df. Returns row indices.
# ---------------------------------------------------------------------------
.match_level_rows <- function(master_df, level) {
  rows <- which(master_df$Test_Family == level$Test_Family)
  if (length(rows) == 0) return(integer(0))

  if (!is.null(level$split_by) && nzchar(level$split_by)) {
    # Split_By string may contain " + " for multi-variable splits. Match the
    # WHOLE Split_By column against the user's literal request.
    keep <- master_df$Split_By[rows] == level$split_by
    rows <- rows[keep]
  }
  rows
}


# ---------------------------------------------------------------------------
# Public: assign tree metadata columns to master_df given a tree_spec.
# `var_names` is the experimental design's variable ordering (e.g.
# c("Genotype", "Drug")) and is forwarded into .extract_link_value() so
# linkage values can be parsed positionally out of Full_Interaction contrast
# strings.
# ---------------------------------------------------------------------------
assign_tree_metadata <- function(master_df, tree_spec, var_names = character(0),
                                 condition_lookup = NULL) {
  n <- nrow(master_df)
  master_df$tree_id        <- NA_character_
  master_df$tree_level     <- NA_integer_
  master_df$tree_parent_id <- NA_character_

  if (is.null(tree_spec) || length(tree_spec) == 0) return(master_df)

  # First pass: assign tree_level + a unique tree_id per matched row.
  for (Lidx in seq_along(tree_spec)) {
    level <- tree_spec[[Lidx]]
    rows  <- .match_level_rows(master_df, level)
    if (length(rows) == 0) next

    # Build a deterministic id per row: "L{level}__{Passed_Contrasts}__G{Group}__{Split_By}"
    new_ids <- sprintf("L%d__%s__G%s__%s",
                       Lidx,
                       master_df$Passed_Contrasts[rows],
                       as.character(master_df$Group[rows]),
                       master_df$Split_By[rows])
    master_df$tree_id[rows]    <- new_ids
    master_df$tree_level[rows] <- Lidx
  }

  # Second pass: link each row at level L > 1 to its parent at level L-1
  # via the linkage_key. Phase (Group) ALWAYS has to match as well — a child
  # in phase 3 can only be gated by a parent in phase 3, never a parent in
  # phase 1, regardless of what the linkage_key is. This is essential for
  # the rescue paradigm: per-(phase, drug) gating means a Drug_X specificity
  # test in phase 3 must wait for the rescue claim in phase 3, NOT the one
  # in phase 1. The L1→L2 default linkage of "Group" already implies phase
  # matching; this additional check ensures phase matching survives when the
  # linkage_key is a design variable (e.g. Drug) instead.
  for (Lidx in seq_along(tree_spec)) {
    if (Lidx == 1L) next                       # L1 has no parent
    level   <- tree_spec[[Lidx]]
    linkage <- level$linkage_key
    if (is.null(linkage)) next

    child_var  <- linkage$child_var  %||% "Group"
    parent_var <- linkage$parent_var %||% child_var

    child_rows  <- which(master_df$tree_level == Lidx)
    parent_rows <- which(master_df$tree_level == (Lidx - 1L))
    if (length(child_rows) == 0 || length(parent_rows) == 0) next

    for (cr in child_rows) {
      child_value <- .extract_link_value(master_df, cr, child_var,
                                         var_names, condition_lookup)
      # Skip if no linkage info (all elements NA or zero-length)
      if (length(child_value) == 0L || all(is.na(child_value))) next
      child_group <- as.character(master_df$Group[cr])

      # Walk parent rows, finding one with overlapping linkage values AND
      # matching Group (phase). Set-overlap (not equality) lets design-
      # variable linkage handle pairwise vs trt.vs.ctrl directionality:
      # L2 Drug_Effect's "compound - DMSO" returns {compound, DMSO} and L3
      # Full_Interaction's "WT+DMSO - WT+compound" also returns {DMSO,
      # compound} at the Drug position. The intersection is non-empty so
      # they match — without caring which side was LHS vs RHS.
      for (pr in parent_rows) {
        parent_value <- .extract_link_value(master_df, pr, parent_var,
                                            var_names, condition_lookup)
        if (length(parent_value) == 0L || all(is.na(parent_value))) next
        parent_group <- as.character(master_df$Group[pr])
        if (.link_values_overlap(child_value, parent_value) &&
            identical(child_group, parent_group)) {
          master_df$tree_parent_id[cr] <- master_df$tree_id[pr]
          break
        }
      }
    }
  }

  master_df
}


# ---------------------------------------------------------------------------
# Return-type contract: .extract_link_value() returns a character vector
# (length 1 or 2). The matcher then uses set-overlap (intersect of length>0)
# to decide whether two rows link.
#
# Why a vector: pairwise-style contrasts (Full_Interaction) have arbitrary
# LHS/RHS direction, while trt.vs.ctrl-style (per-variable Effect) always
# put treatment on LHS. A contrast varying along a variable carries BOTH
# values; returning the set of varying values + matching by intersection
# lets the two formats link without depending on directionality.
#
# Cases:
#   "Group"               → master_df$Group[row] (length 1)
#   "Passed_Contrasts" /
#   "Tested_Level" (alias) → full contrast string (length 1)
#   any column in master_df → master_df[[col]][row] (length 1)
#   one of var_names      → DESIGN-VARIABLE PARSER. Uses condition_lookup
#                              (Python-supplied) for multi-variable contrasts
#                              and Test_Family / Split_By for single-variable.
#                              Returns:
#                              • length-2 c(LHS, RHS) when contrast varies
#                                along this variable
#                              • length-1 split-by value when held constant
#                                via Split_By context
#                              • NA when contrast doesn't involve the variable
#   anything else         → NA (no fallback string parsing)
# ---------------------------------------------------------------------------
.extract_link_value <- function(master_df, row, var_name, var_names = character(0),
                                condition_lookup = NULL) {
  if (is.null(var_name) || is.na(var_name) || !nzchar(var_name)) {
    return(NA_character_)
  }
  if (var_name == "Group") return(as.character(master_df$Group[row]))
  # Treat "Tested_Level" as the legacy alias for the renamed column
  if (var_name %in% c("Passed_Contrasts", "Tested_Level")) {
    return(as.character(master_df$Passed_Contrasts[row]))
  }
  if (var_name %in% colnames(master_df)) {
    return(as.character(master_df[[var_name]][row]))
  }
  # Design-variable parsing branch — needs condition_lookup for multi-variable
  # contrasts (Full_Interaction). Falls through to NA for an unknown linkage.
  if (length(var_names) > 0 && var_name %in% var_names) {
    return(.parse_design_var_from_contrast(master_df, row, var_name,
                                           var_names, condition_lookup))
  }
  NA_character_
}


# ---------------------------------------------------------------------------
# Design-variable parser. Given a row in master_df and a target design
# variable, return the value of that variable for this contrast — or NA if
# the contrast doesn't differentiate along that variable.
#
# Single-variable contrasts (e.g. Drug_Effect's "Drug_X - DMSO"):
#   • If Test_Family == <var>_Effect: return the LHS (treatment side)
#   • If Test_Family is for a different variable: try the Split_By column
#     to find this variable's value in the split-by context
#
# Multi-variable contrasts (Full_Interaction's "WT+Drug_X - WT+DMSO"):
#   • Split LHS and RHS on "+" using var_names as positional reference
#   • If the variable's LHS value == RHS value, it's held constant → NA
#   • Otherwise return the LHS value (treatment side of the difference)
# ---------------------------------------------------------------------------
.parse_design_var_from_contrast <- function(master_df, row, var_name, var_names,
                                            condition_lookup = NULL) {
  tl <- as.character(master_df$Passed_Contrasts[row])
  parts <- strsplit(tl, " - ", fixed = TRUE)[[1]]
  if (length(parts) != 2L) return(NA_character_)
  lhs <- trimws(parts[1])
  rhs <- trimws(parts[2])

  # CASE A — multi-variable contrast: lhs/rhs are full condition names like
  # "WT+DMSO". Decode via the Python-supplied lookup (no regex parsing).
  if (!is.null(condition_lookup) &&
      !is.null(condition_lookup[[lhs]]) &&
      !is.null(condition_lookup[[rhs]])) {
    lhs_val <- condition_lookup[[lhs]][[var_name]]
    rhs_val <- condition_lookup[[rhs]][[var_name]]
    if (is.null(lhs_val) || is.null(rhs_val)) return(NA_character_)
    if (lhs_val == rhs_val) return(NA_character_)   # variable held constant
    return(c(as.character(lhs_val), as.character(rhs_val)))
  }

  # CASE B — single-variable contrast (e.g. Drug_Effect's "compound - DMSO"):
  # lhs/rhs ARE values of the focal variable. No condition lookup needed.
  test_family <- as.character(master_df$Test_Family[row])
  focal_var <- sub("_Effect$", "", test_family)
  if (focal_var == var_name) {
    return(c(lhs, rhs))
  }

  # CASE C — var_name is in the row's split_by context. Pull the value from
  # the Split_By column positionally.
  if ("Split_By" %in% colnames(master_df)) {
    sb <- as.character(master_df$Split_By[row])
    if (!sb %in% c("", "None", "All", "(none)")) {
      other_vars <- setdiff(var_names, focal_var)
      sb_parts <- strsplit(sb, " + ", fixed = TRUE)[[1]]
      idx <- which(other_vars == var_name)
      if (length(idx) == 1L && idx <= length(sb_parts)) {
        return(trimws(sb_parts[idx]))
      }
    }
  }

  NA_character_
}


# ---------------------------------------------------------------------------
# Set-overlap match: child and parent share at least one linkage value.
# Used by assign_tree_metadata() to decide if a child row links to a parent.
# Empty / all-NA inputs never match.
# ---------------------------------------------------------------------------
.link_values_overlap <- function(child_value, parent_value) {
  if (length(child_value) == 0L || length(parent_value) == 0L) return(FALSE)
  cv <- as.character(child_value)
  pv <- as.character(parent_value)
  cv <- cv[!is.na(cv)]
  pv <- pv[!is.na(pv)]
  if (length(cv) == 0L || length(pv) == 0L) return(FALSE)
  length(intersect(cv, pv)) > 0L
}


# ---------------------------------------------------------------------------
# Build the DEFAULT rescue-paradigm tree spec from the run's variable names.
# Used when correction.json sidecar is missing — the rescue tree pre-populates
# (matches what the UI's "first open" will offer).
# ---------------------------------------------------------------------------
default_rescue_tree_spec <- function(var_names, ref_map) {
  if (length(var_names) < 2) return(NULL)
  # Convention: first variable is the gating variable (e.g. Genotype),
  # second variable is the screened variable (e.g. Drug).
  gate_var   <- var_names[1]
  screen_var <- var_names[2]
  ref_gate   <- ref_map[[gate_var]]    # e.g. "WT"
  ref_screen <- ref_map[[screen_var]]  # e.g. "DMSO"

  # Default tree (rescue paradigm), simpler → more complex hypothesis order:
  #   L1 = gate (Genotype_Effect in DMSO context)
  #     "Does the disease phenotype exist?"
  #   L2 = effect (Drug_Effect, both genotype contexts)
  #     "Does the drug do anything?"
  #   L3 = specificity (Full_Interaction)
  #     "Is the drug's effect specifically rescue (differs by genotype)?"
  # All inter-level linkages use Group → per-phase gating throughout.
  list(
    list(
      name        = paste0("L1 ", gate_var, " gate"),
      Test_Family = paste0(gate_var, "_Effect"),
      split_by    = ref_screen,            # gate in the reference-screen context
      linkage_key = NULL
    ),
    list(
      name        = paste0("L2 ", screen_var, " effect"),
      Test_Family = paste0(screen_var, "_Effect"),
      split_by    = NULL,                  # both genotype contexts
      linkage_key = list(parent_var = "Group", child_var = "Group")
    ),
    list(
      name        = "L3 rescue specificity",
      Test_Family = "Full_Interaction",
      split_by    = NULL,
      linkage_key = list(parent_var = "Group", child_var = "Group")
    )
  )
}
