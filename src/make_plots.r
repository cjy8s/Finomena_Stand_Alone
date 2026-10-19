# =============================================================================
# Initial-screen result plots — publication-ready R / ggplot2 figures.
# -----------------------------------------------------------------------------
# The SageMaker R image has no pip/matplotlib, so plotting lives in R (installed
# the same way the BAM script installs its packages). Reads:
#   ALL_PLATES_full_results.csv   model contrasts + equivalence
#   SSMD_vs_model_ALL.csv         model-free effect sizes + model
#   aggregated_data.csv           per-second data w/ bio_plate (trajectories)
#
# Layout (updated for the 14-plate screen): the scatter, waterfall and Check-A
# figures are ONE ggplot each with facet_wrap(~Plate, 3 columns) and a SINGLE
# shared colourbar — so 14 plates read as one coherent figure instead of 14
# stitched panels with 14 legends. The disease window and (Check A) Spearman r
# are drawn as per-panel corner annotations. The trajectory figures stay
# per-panel (each plate's top movers differ, so a shared legend would carry ~42
# entries) and are tiled 3-wide with gridExtra; every trajectory figure renders
# from the SAME traj_panel() code path so their style is identical by
# construction. All plots are PER PLATE (the batch confound means cross-plate
# magnitudes aren't comparable). ggplot2 is upgraded to >= 3.5 at load so
# ggnewscale works (two colour scales on the trajectory panels).
# Env paths mirror the rest of the pipeline (FINOMENA_ADDON_PROD_DIR…).
# =============================================================================

# ── PACKAGE LIBRARY (mirror initial_screen_script.r) ─────────────────────────
.conda_prefix <- Sys.getenv("CONDA_PREFIX")
.lib <- if (nzchar(.conda_prefix)) {
  file.path(.conda_prefix, "lib", "R", "library")
} else {
  file.path(Sys.getenv("HOME", unset = "/tmp"), "R", "library")
}
dir.create(.lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(.lib, .libPaths()))
cran_mirror <- Sys.getenv("CRAN_MIRROR", unset = "https://cloud.r-project.org")

# CRITICAL — do NOT load any package (library/requireNamespace) until every
# install below is done. Version-check with .have(), which reads the DESCRIPTION
# metadata and does NOT load the namespace. The trap this avoids: requireNamespace
# ("ggplot2") loads the image's OLD ggplot2 just to check its version, and with it
# pulls scales 1.3.0 into the session. The current Posit mirror serves ggplot2
# 4.x, which Imports scales >= 1.4.0 — but a namespace once loaded can't be
# upgraded in-session, so the new ggplot2 then fails to load ("namespace 'scales'
# 1.3.0 is already loaded, but >= 1.4.0 is required"). Checking versions without
# loading, and deferring all library() calls to the very end, sidesteps it.
.have <- function(pkg, min_version = NULL) {
  v <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
  if (is.null(v)) return(FALSE)
  is.null(min_version) || v >= min_version
}
# The base image bakes in a pre-3.5 ggplot2 that blocks modern helpers
# (ggnewscale needs ggplot2 >= 3.5.0). Install the newest ggplot2 into the user
# lib (prepended to .libPaths above), which shadows the image copy — binary from
# the Posit mirror, no compilation. install.packages pulls its up-to-date deps
# (scales, rlang, S7) into the same lib.
if (!.have("ggplot2", "3.5.0")) {
  cat("Upgrading ggplot2 to >= 3.5.0 (image ships an older copy)\n")
  install.packages("ggplot2", repos = cran_mirror)
}
# Belt-and-suspenders: ensure a scales new enough for ggplot2 4.x is present in
# the user lib (normally pulled as a ggplot2 dep above; installed explicitly here
# in case ggplot2 was already >= 3.5 but the image's scales is stale).
if (!.have("scales", "1.4.0")) {
  cat("Installing scales >= 1.4.0 (required by ggplot2 4.x)\n")
  install.packages("scales", repos = cran_mirror)
}
for (pkg in c("data.table", "gridExtra", "ggnewscale")) {
  if (!.have(pkg)) {
    cat("Installing ", pkg, "\n", sep = "")
    install.packages(pkg, repos = cran_mirror)
  }
}
# NOTE: the top-hit labels (figures 1b–1e) are placed DETERMINISTICALLY with
# geom_segment leaders + geom_text (see attach_labels() / rescue_facet()), NOT
# ggrepel. ggrepel de-overlaps label BOXES via a physics sim that knows nothing
# about leader LINES, so with these tightly-clustered diagonal picks it kept
# shoving a label out of the label row and neighbouring leaders slashed across it.
# The straight-leader fan below is crossing-free by construction, so ggrepel (and
# its optional install / failure mode) is gone entirely.
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(gridExtra); library(grid)
  library(ggnewscale)
})
options(width = 320)

# ── PATHS ────────────────────────────────────────────────────────────────────
PROD_DIR <- Sys.getenv("FINOMENA_ADDON_PROD_DIR",
                       unset = "Data/Initial screen/Production Output")
OUT      <- Sys.getenv("FINOMENA_PLOTS_DIR",
                       unset = "Data/Initial screen/Plots")
AGG_FULL <- Sys.getenv("FINOMENA_PLOT_AGG_CSV",
                       unset = file.path(PROD_DIR, "aggregated_data.csv"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
# Clear stale PNGs so the Plots dir reflects EXACTLY the current figure set.
# make_plots.r is the only writer here and regenerates every figure each run, so
# any pre-existing .png is from an older version (e.g. seeded resume output still
# carrying the removed 5_/6_ limma figures). Leaving them would misrepresent the
# run as producing figures it no longer makes.
invisible(file.remove(list.files(OUT, pattern = "\\.png$", full.names = TRUE)))
POS_CONTROL <- Sys.getenv("FINOMENA_POS_CONTROL", unset = "KO+DMSO")
NEG_CONTROL <- Sys.getenv("FINOMENA_NEG_CONTROL", unset = "WT+DMSO")
# Match rescue_addons.r's ashr toggle. When OFF, the *_ashr_* columns hold the
# raw (unshrunken) values, so the Check figures show ONLY the raw series — no
# redundant duplicate "ashr" series/legend entry.
USE_ASHR <- tolower(Sys.getenv("FINOMENA_USE_ASHR", unset = "0")) %in%
            c("1", "true", "yes", "on")

full <- fread(file.path(PROD_DIR, "ALL_PLATES_full_results.csv"))
# Augmented per-plate table (movement effect + Glass's delta + landing/equiv per
# drug). Used by Check A for the robust Glass's delta column.
aug  <- fread(file.path(PROD_DIR, "augmented_per_plate.csv"))

# ── PLATE LIST + NATURAL ORDER ───────────────────────────────────────────────
# The plate list is derived from the data actually present (single source of
# truth), NOT a hardcoded count — run_pipeline.sh derives FINOMENA_PLATES from
# the BAM output folders, but we also fall back to the plates present in the
# results so a hand run covers everything. `mixedsort`-style natural order keeps
# "Plate 10" after "Plate 9" (plain sort would put it after "Plate 1").
.natural_order <- function(x) {
  x <- unique(as.character(x))
  key <- suppressWarnings(as.integer(gsub("\\D+", "", x)))
  x[order(key, x)]
}
.env_plates <- trimws(strsplit(Sys.getenv("FINOMENA_PLATES", unset = ""), ",")[[1]])
.env_plates <- .env_plates[nzchar(.env_plates)]
PLATES <- if (length(.env_plates)) .env_plates else .natural_order(full$Plate)
cat("Plates (", length(PLATES), "): ", paste(PLATES, collapse = ", "), "\n", sep = "")

# Facet order is driven by a factor whose levels are the natural-ordered plates.
FP <- function(x) factor(as.character(x), levels = PLATES)
full[, Plate := FP(Plate)]
aug[,  Plate := FP(Plate)]

# ── SHARED PUBLICATION THEME ─────────────────────────────────────────────────
STRIP_FILL <- "#37474F"
theme_pub <- function(base = 12) {
  theme_minimal(base_size = base, base_family = "") +
    theme(
      plot.title    = element_text(face = "bold", size = rel(1.0),
                                   hjust = 0.5, margin = margin(b = 6)),
      strip.background = element_rect(fill = STRIP_FILL, colour = NA),
      # Facet strip ("Plate X") enlarged for publication legibility.
      strip.text    = element_text(colour = "white", face = "bold", size = rel(1.35),
                                   margin = margin(5, 5, 5, 5)),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      panel.border      = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),
      # Axis titles + tick labels sized up so a shrunk-to-column-width figure stays
      # readable in print.
      axis.title    = element_text(colour = "grey20", size = rel(1.2)),
      axis.text     = element_text(size = rel(1.0)),
      legend.position = "right",
      # Colourbar title sits to the RIGHT of the bar, rotated parallel to it and
      # centred along its height — one line, no wrapping. Enlarged to match the
      # bigger bar.
      legend.title  = element_text(face = "bold", size = rel(1.2),
                                   angle = -90, hjust = 0.5),
      legend.text   = element_text(size = rel(1.05)),
      legend.key.height = unit(1.4, "lines"),
      legend.key.width  = unit(1.1, "lines"),
      plot.margin   = margin(6, 8, 6, 8)
    )
}
# Wrap a plot (or a list of per-plate panels) under a shared CENTERED bold title
# + grey subtitle on a white background. `panels` may be a single faceted ggplot
# (pass list(p), ncol = 1) or a list of per-plate panels tiled ncol-wide. Returns
# a grob suitable for save_fig().
assemble <- function(panels, title, subtitle = NULL, ncol = 1) {
  nrow <- ceiling(length(panels) / ncol)
  # Publication-scale header: larger bold title + larger grey subtitle. The band
  # heights are widened to match so the enlarged text never crowds the top panels.
  title_grob <- grid::textGrob(title, gp = grid::gpar(fontface = "bold", cex = 1.9),
                               hjust = 0.5, x = 0.5)
  grobs <- list(title_grob)
  heights <- unit(2.8, "lines")
  if (!is.null(subtitle)) {
    grobs <- c(grobs, list(grid::textGrob(subtitle,
                 gp = grid::gpar(col = "grey35", cex = 1.15), hjust = 0.5, x = 0.5)))
    heights <- grid::unit.c(heights, unit(2.2, "lines"))
  }
  body <- gridExtra::arrangeGrob(grobs = panels, nrow = nrow, ncol = ncol)
  gridExtra::arrangeGrob(
    grobs   = c(grobs, list(body)),
    heights = grid::unit.c(heights, unit(1, "null")),
    ncol    = 1)
}
# Save a grob to PNG with an explicit white background. Uses the png device
# directly + grid.draw so it works for gridExtra arrangeGrob output across
# ggplot2 vintages (ggsave on a gtable can be finicky on older versions).
save_fig <- function(name, grob, width, height) {
  png(file.path(OUT, name), width = width, height = height, units = "in",
      res = 200, bg = "white")
  grid::grid.newpage()
  grid::grid.draw(grob)
  dev.off()
}

PR_LIMITS <- c(0, 0.5)
# Large, legible colourbar; title on the RIGHT, rotated parallel to the bar and
# centred along its height (single line). Plain-language title (no code syntax).
# One shared instance across every faceted figure -> a single legend per figure.
# Bar enlarged (taller + noticeably wider) so it reads clearly next to the 3-wide,
# 5-row facet grid; tick labels enlarged to match.
viridis_fill <- scale_fill_viridis_c(
  name = "Probability of equivalence to wild-type", limits = PR_LIMITS, option = "D",
  oob = scales::squish,
  guide = guide_colourbar(barheight = unit(9, "cm"), barwidth = unit(1.1, "cm"),
                          title.position = "right", title.hjust = 0.5,
                          ticks.colour = "grey30", frame.colour = "grey50"))

# ── DISEASE GAP per plate (= WT - KO). DrugEffect_vsKO - Rescue_vsWT = WT - KO.
geno <- sapply(PLATES, function(p) {
  rs <- full[Plate == p]
  median(rs$DrugEffect_vsKO_est_log2 - rs$Rescue_vsWT_est_log2, na.rm = TRUE)
}); names(geno) <- PLATES
# Trajectory panel title as a plotmath expression so the log2 subscript renders
# as a real subscript (the Unicode "₂" has no glyph in the container font and
# showed up as a box; plotmath log[2] is a true typeset subscript that works).
panel_title <- function(p)
  bquote(.(p) ~ "·  disease window =" ~ .(sprintf("%.2f", geno[[p]])) ~ log[2]*"FC")

# Per-plate data frames the faceted figures need (one row per plate). The band /
# KO-line depend on that plate's disease gap, so they are supplied as data (an
# annotate() would draw the SAME gap in every facet — wrong). The corner label
# carries the disease window as a parseable plotmath string (real log2 subscript).
band <- data.table(Plate = FP(PLATES), gap = as.numeric(unlist(geno[PLATES])))
gap_lab <- data.table(
  Plate = FP(PLATES),
  lab = sprintf('"disease window %s"~log[2]*"FC"',
                sprintf("%.2f", as.numeric(unlist(geno[PLATES])))))

# ── CANDIDATE SELECTORS ──────────────────────────────────────────────────────
# top_drugs() returns the top-N drugs on a plate under an ordering function
# (controls excluded). Defined once and reused for BOTH the rescue-scatter labels
# (figures 1c/1d/1e) and the trajectory-curve picks (figures 3/3b/3c), so the two
# views always agree on "the top 3".
#   ORDER_CLOSENESS  : smallest |Rescue_vsWT| (landed nearest WT), tie-broken by
#                      highest posterior equivalence to WT. The rescue criterion.
#   ORDER_MAX_EFFECT : largest positive DrugEffect_vsKO (biggest movement off the
#                      KO baseline — the waterfall's ranking), same tie-break.
#   ORDER_EQUIV      : highest posterior equivalence to WT (equiv_Pr_vsWT — the
#                      point-colour metric), tie-broken by closeness to WT. The
#                      model's own similarity-to-WT score, vs the point-estimate
#                      distance ORDER_CLOSENESS uses.
top_drugs <- function(p, order_fun, n = 3) {
  d <- full[Plate == p & !drug %in% c(POS_CONTROL, NEG_CONTROL)]
  head(d[order_fun(d)]$drug, n)
}
ORDER_CLOSENESS  <- function(d) order(abs(d$Rescue_vsWT_est_log2), -d$equiv_Pr_vsWT)
ORDER_MAX_EFFECT <- function(d) order(-d$DrugEffect_vsKO_est_log2, -d$equiv_Pr_vsWT)
ORDER_EQUIV      <- function(d) order(-d$equiv_Pr_vsWT, abs(d$Rescue_vsWT_est_log2))

# Per-plate axis extents, used to place labels a PROPORTIONATE distance off their
# anchor on every panel (plates span very different ranges under scales="free", so
# a fixed nudge would be huge on one plate and invisible on another). Carries min/
# max as well as span so attach_labels() can build an in-panel horizontal fan.
plate_ranges <- full[, .(
  xmin = min(DrugEffect_vsKO_est_log2, na.rm = TRUE),
  xmax = max(DrugEffect_vsKO_est_log2, na.rm = TRUE),
  ymin = min(Rescue_vsWT_est_log2,     na.rm = TRUE),
  ymax = max(Rescue_vsWT_est_log2,     na.rm = TRUE)), by = Plate]
plate_ranges[, xr := xmax - xmin][, yr := ymax - ymin]

# Compute an explicit, GUARANTEED-NON-CROSSING label position (target_x/target_y)
# for every hit, per plate. The top picks cluster tightly on the diagonal near the
# WT line, so the labels must be lifted off the cloud into the open upper region.
# The crossing-free layout is purely geometric — no physics engine (this is why
# ggrepel was dropped: its box de-overlap moved labels off the row and let leaders
# cross). The rule:
#   1. All labels sit in ONE horizontal ROW at a common height ABOVE every anchor,
#      so every leader points strictly upward.
#   2. Targets are fanned evenly across the panel width and assigned so their
#      left→right order MATCHES the anchors' left→right order.
# Two upward leaders whose TOP endpoints preserve their BOTTOM endpoints' x-order
# cannot cross — and because rescue_facet() draws the leaders as straight
# geom_segments to exactly these targets (nothing perturbs them afterward), the
# guarantee holds in the rendered figure, not just in intent.
attach_labels <- function(t) {
  t <- copy(t)
  t[, well := sub("KO\\+", "", drug)]
  t <- merge(t, plate_ranges, by = "Plate", all.x = TRUE, sort = FALSE)
  t[, `:=`(target_x = NA_real_, target_y = NA_real_)]
  for (p in unique(t$Plate)) {
    idx <- which(t$Plate == p)
    k <- length(idx)
    xr <- t$xr[idx][1]; yr <- t$yr[idx][1]
    xmin <- t$xmin[idx][1]; xmax <- t$xmax[idx][1]
    if (!is.finite(xr) || xr == 0) xr <- 1
    if (!is.finite(yr) || yr == 0) yr <- 1
    ax <- t$DrugEffect_vsKO_est_log2[idx]   # anchor x
    ay <- t$Rescue_vsWT_est_log2[idx]        # anchor y
    # Rank anchors left→right (ties by y) and lay the label targets in that SAME
    # order, so target x-order == anchor x-order — the non-crossing guarantee.
    ord <- order(ax, ay)
    # Fan the targets horizontally, spaced a FIXED ~0.30*x-span apart so a 5-char
    # well ID never touches its neighbour (labels render at fixed size, so the gap
    # must be a real fraction of panel width). The fan is a RIGID block: we place it
    # at even spacing then SHIFT the whole block (never squeeze it) to sit inside the
    # panel. Per-label clamping was wrong — when the top picks are the rightmost
    # points (the usual case: best rescues move most), clamping crushed the fan
    # against the right edge and labels overlapped. A rigid shift keeps every gap at
    # 0.30*xr no matter where the anchors sit; the empty upper-LEFT region absorbs it.
    step <- 0.30 * xr
    tx   <- mean(ax) + (seq_len(k) - (k + 1) / 2) * step   # evenly spaced, centred on anchors
    lo <- xmin - 0.05 * xr; hi <- xmax + 0.05 * xr          # panel's expanded x-range
    if (tx[k] > hi) tx <- tx - (tx[k] - hi)                 # block overflows right → shift left
    if (tx[1] < lo) tx <- tx + (lo - tx[1])                 # block overflows left  → shift right
    # (If the fan is wider than the panel it can't fit either way; the shifts above
    # still keep it roughly centred and the spacing intact — headroom expansion below
    # widens the panel enough that k=3 always fits.)
    # Common target height: a clear margin ABOVE the highest anchor (not the whole
    # plate), so labels sit in a tidy row above the picks and leaders point up. The
    # top headroom expansion (0.34*yr above the plate max) guarantees room for it.
    ty <- max(ay) + 0.18 * yr
    t$target_x[idx[ord]] <- tx
    t$target_y[idx[ord]] <- ty
  }
  t[!is.finite(target_x), target_x := DrugEffect_vsKO_est_log2]
  t[!is.finite(target_y), target_y := Rescue_vsWT_est_log2]
  t[, Plate := FP(Plate)]
  t[]
}
# Build a ring/label table for the top-N drugs per plate under an ordering — the
# shape rescue_facet()'s `hits` expects.
top_hits_table <- function(order_fun, n = 3) {
  tbl <- rbindlist(lapply(PLATES, function(p) full[Plate == p & drug %in% top_drugs(p, order_fun, n)]))
  attach_labels(tbl)
}

# =============================================================================
# 1. RESCUE SCATTER — movement (x) vs landing (y), coloured by P(equiv to WT)
# =============================================================================
# ONE faceted figure. `hits` is an optional data.table of wells to call out with
# a ring + well-ID label (needs `well`, the x/y anchor cols and target_x/target_y
# from attach_labels(), all carrying a Plate factor so each facet is laid out
# independently); NULL leaves the plain scatter. Factored out so the plain and
# every labeled figure render from ONE code path — a labeled figure is figure 1
# plus annotations.
#   ring_colour : colour of the ring, leader segment, and label text. Crimson
#                 (#B71C1C) = confirmed hit (fig 1b); red (#D32F2F) = top-ranked
#                 call-out (figs 1c/1d/1e).
#   headroom    : pad the axes (esp. the top) so the label row has empty room
#                 INSIDE the panel. Labels are clipped to the panel, so the
#                 headroom is what keeps them from crowding the topmost points or
#                 spilling across the panel border.
rescue_facet <- function(hits = NULL, ring_colour = "#B71C1C", headroom = FALSE) {
  # Both axes are RAW point estimates. Movement and landing differ by a fixed
  # per-plate constant (the disease gap: movement - landing = WT - KO), so they
  # are algebraically locked on a 1:1 diagonal. ashr-shrinking only the movement
  # axis (landing must never be shrunk — that biases the WT-equivalence call)
  # would break that exact relationship and slide points off the true diagonal
  # for non-biological reasons. So the scatter shows raw/raw; ashr is represented
  # where it is paired and meaningful (Check A) and via the confidence s-values.
  pl <- ggplot(full, aes(DrugEffect_vsKO_est_log2, Rescue_vsWT_est_log2)) +
    # Disease-window band + KO line are per-plate (data-driven, not annotate()).
    geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = -gap, ymax = 0),
              fill = "#FFB300", alpha = 0.10, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, colour = "#2E7D32", linetype = "dashed", linewidth = 0.5) +
    geom_hline(data = band, aes(yintercept = -gap), colour = "#C62828",
               linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(aes(fill = equiv_Pr_vsWT), shape = 21, colour = "grey20",
               size = 1.8, stroke = 0.2, alpha = 0.9) +
    viridis_fill +
    # Disease window in each panel's BOTTOM-RIGHT corner (empty region: points lie
    # on the landing = movement - gap diagonal, so the lower-right triangle is
    # always clear, and the labeled top picks sit upper-right — so this corner
    # never collides). parse=TRUE renders the real log2 subscript; x/y = Inf/-Inf
    # pins it to each facet's own bottom-right corner under free scales; hjust≈1 /
    # vjust<0 inset it just inside the panel.
    geom_text(data = gap_lab, aes(x = Inf, y = -Inf, label = lab), parse = TRUE,
              inherit.aes = FALSE, hjust = 1.03, vjust = -1.0, size = 4.6,
              fontface = "bold", colour = "grey30") +
    facet_wrap(~ Plate, ncol = 3, scales = "free") +
    labs(x = expression(Log[2]*"FC effect estimate from KO with DMSO"),
         y = expression(Log[2]*"FC effect estimate from wild-type with DMSO")) +
    theme_pub()

  # Extra breathing room so the label row seats cleanly INSIDE the panel. More at
  # the TOP, where the strongest candidates cluster.
  if (headroom)
    pl <- pl +
      scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
      # Generous TOP room: labels target ~0.18*y-span above the highest anchor
      # (often the plate max), and the text needs a little more above that — so
      # expand the top by 0.34*y-span to keep the label row off the panel border.
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.34)))

  if (!is.null(hits) && nrow(hits) > 0) {
    # Unfilled ring drawn OVER the point (viridis fill still reads through) plus a
    # straight leader + well-ID label. Drawn last so ring + label sit above the cloud.
    pl <- pl +
      geom_point(data = hits,
                 aes(DrugEffect_vsKO_est_log2, Rescue_vsWT_est_log2),
                 shape = 21, fill = NA, colour = ring_colour,
                 size = 4.2, stroke = 1.1, inherit.aes = FALSE)
    # DETERMINISTIC leaders + labels (no ggrepel). attach_labels() has placed every
    # label at target_x/target_y — a single horizontal row above the anchors, in
    # anchor x-order — so the leaders form a clean upward fan that CANNOT cross (two
    # upward segments preserving their base x-order can't intersect). We draw exactly
    # that: a thin straight leader from each anchor to its label, then the well ID at
    # the target. Nothing repositions them afterward, so what's computed is what's
    # drawn. vjust = -0.3 seats the text just ABOVE target_y, so the leader tip meets
    # the bottom of the label instead of poking through it.
    pl <- pl +
      geom_segment(data = hits,
                   aes(x = DrugEffect_vsKO_est_log2, y = Rescue_vsWT_est_log2,
                       xend = target_x, yend = target_y),
                   inherit.aes = FALSE, colour = ring_colour,
                   linewidth = 0.35, alpha = 0.9) +
      geom_text(data = hits,
                aes(x = target_x, y = target_y, label = well),
                inherit.aes = FALSE, colour = ring_colour, fontface = "bold",
                size = 3.7, vjust = -0.3)
  }
  pl
}

# Figure sizes for the 3-wide faceted grids (14 plates → 5 rows). Generous so
# each panel is large enough to seat the label row without crowding.
FACET_W <- 15; FACET_H <- 20

p1 <- assemble(list(rescue_facet()),
  "Rescue triage: movement off the disease baseline vs landing at wild-type",
  "Each point is a compound; rescues sit upper-right. Green = wild-type, red = knockout.")
save_fig("1_rescue_scatter.png", p1, FACET_W, FACET_H)

# =============================================================================
# 1b. RESCUE SCATTER + CONFIRMED HITS — identical to figure 1, with a crimson
#     ring and well-ID label on wells later validated as true rescues.
# -----------------------------------------------------------------------------
# `drug` is stored as "KO+<well>"; the well ID is what the bench refers to. Add
# more confirmed hits here as they are validated (plates with no confirmed hit
# are simply left unlabeled).
CONFIRMED_HITS <- data.table(
  Plate = c("Plate 1", "Plate 3"),
  well  = c("01F05",   "03E06"))
CONFIRMED_HITS[, drug := paste0("KO+", well)]
# Key on plate+drug as plain strings to sidestep any factor/character merge
# coercion, then report any confirmed hit absent from the results (typo or a
# plate that wasn't fit) rather than silently dropping it.
.hit_key <- CONFIRMED_HITS[, paste(Plate, drug)]
hits_full <- full[paste(as.character(Plate), drug) %in% .hit_key]
found_key <- hits_full[, paste(as.character(Plate), drug)]
missing <- CONFIRMED_HITS[!.hit_key %in% found_key]
if (nrow(missing) > 0)
  cat("WARNING: confirmed hit(s) not found in results, cannot label: ",
      paste(missing$Plate, missing$well, sep = "/", collapse = ", "), "\n", sep = "")
hits_full <- attach_labels(hits_full)

p1b <- assemble(list(rescue_facet(hits = hits_full, headroom = TRUE)),
  "Rescue triage with confirmed hits highlighted",
  "Ringed = validated rescues. Green = wild-type, red = knockout.")
save_fig("1b_rescue_scatter_hits.png", p1b, FACET_W, FACET_H)

# =============================================================================
# 1c/1d/1e. RESCUE SCATTER + TOP-3 LABELS — identical to figure 1, with each
#     plate's top 3 compounds ringed and labeled by well ID. Three rankings, one
#     figure each: labels sit in a straight-leader fan above the picks (see
#     attach_labels()), so no leader crosses a label. Rankings: 1c closeness to WT,
#     1d movement effect, 1e equivalence.
# =============================================================================
p1c <- assemble(list(rescue_facet(hits = top_hits_table(ORDER_CLOSENESS),
                                   ring_colour = "#D32F2F", headroom = TRUE)),
  "Rescue triage: top 3 per plate by closeness to wild-type",
  "Ringed = 3 compounds landing nearest wild-type. Green = wild-type, red = knockout.")
save_fig("1c_rescue_scatter_top3_closeness.png", p1c, FACET_W, FACET_H)

p1d <- assemble(list(rescue_facet(hits = top_hits_table(ORDER_MAX_EFFECT),
                                   ring_colour = "#D32F2F", headroom = TRUE)),
  "Rescue triage: top 3 per plate by movement effect",
  "Ringed = 3 strongest movers off the knockout baseline. Green = wild-type, red = knockout.")
save_fig("1d_rescue_scatter_top3_max_effect.png", p1d, FACET_W, FACET_H)

# 1e. Top 3 by equivalence to WT — highest equiv_Pr_vsWT, the model's posterior
# similarity-to-WT score (the point-colour metric). This ranks by the model's
# confidence that a compound is indistinguishable from wild-type, as opposed to
# 1c's raw point-estimate distance; the two usually agree but can diverge when a
# point sits near WT with a wide posterior (small distance, modest probability).
p1e <- assemble(list(rescue_facet(hits = top_hits_table(ORDER_EQUIV),
                                   ring_colour = "#D32F2F", headroom = TRUE)),
  "Rescue triage: top 3 per plate by equivalence to wild-type",
  "Ringed = 3 compounds most likely equivalent to wild-type. Green = wild-type, red = knockout.")
save_fig("1e_rescue_scatter_top3_equivalence.png", p1e, FACET_W, FACET_H)

# =============================================================================
# 2. WATERFALL — drugs ranked by movement, coloured by P(equiv to WT)
# =============================================================================
# Bars = RAW movement effect, matching the scatter (see note there on why the
# displayed point estimates are raw, not ashr-shrunken). Rank is computed WITHIN
# each plate. Y-axis = movement FROM KO (drug - KO): 0 is the KO baseline (red)
# and a full disease gap UP (+gap) reaches WT (green). One faceted figure.
wf <- copy(full)
wf[, rank := frank(DrugEffect_vsKO_est_log2, ties.method = "first"), by = Plate]
p2_plot <- ggplot(wf, aes(rank, DrugEffect_vsKO_est_log2, fill = equiv_Pr_vsWT)) +
  geom_col(width = 1) +
  geom_hline(yintercept = 0, colour = "#C62828", linetype = "dashed", linewidth = 0.5) +
  geom_hline(data = band, aes(yintercept = gap), colour = "#2E7D32",
             linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
  viridis_fill +
  geom_text(data = gap_lab, aes(x = -Inf, y = Inf, label = lab), parse = TRUE,
            inherit.aes = FALSE, hjust = -0.04, vjust = 1.5, size = 2.7, colour = "grey35") +
  facet_wrap(~ Plate, ncol = 3, scales = "free") +
  labs(x = "Compounds, ranked by effect on mutant movement",
       y = expression(log[2]*"FC vs KO  (0 = KO,  +gap = WT)")) +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())
p2 <- assemble(list(p2_plot),
  "Per-compound effect landscape",
  "Each bar is a compound; above zero moves toward wild-type. Colour = P(equivalent to wild-type).")
save_fig("2_waterfall.png", p2, FACET_W, FACET_H)

# NOTE: the old "method agreement" plot (BAM landing vs Glass, vs WT) was removed
# — it is the same concept as Check A (BAM vs Glass on the movement contrast),
# which is the preferred validation. See figure 4 below.

# =============================================================================
# 3 / 3b / 3c. TRAJECTORIES — top movers vs KO & WT mean curves (per-plate).
# -----------------------------------------------------------------------------
# These stay PER-PANEL (not facet_wrap): each plate's top movers are different
# drugs, so a single shared legend would carry ~42 entries. Each panel keeps its
# own two-part legend (Controls / Top movers). All THREE trajectory figures
# render from the SAME traj_panel() code path and the SAME assemble() tiling, so
# their style is identical by construction — they differ only in which drugs are
# selected (closeness / max-effect / equivalence), mirroring the 1c/1d/1e set.
# =============================================================================
agg <- fread(AGG_FULL)
if (!"bio_plate" %in% names(agg)) agg[, bio_plate := sub("_Replicate.*$", "", plate)]
curve <- agg[, .(mv = mean(pixel_diff, na.rm = TRUE)),
             by = .(bio_plate, Condition, time_sec)]
# Reuse top_drugs() so a curve's picks can never drift from its scatter's labels.
top_closeness  <- setNames(lapply(PLATES, top_drugs, order_fun = ORDER_CLOSENESS),  PLATES)
top_max_effect <- setNames(lapply(PLATES, top_drugs, order_fun = ORDER_MAX_EFFECT), PLATES)
top_equiv      <- setNames(lapply(PLATES, top_drugs, order_fun = ORDER_EQUIV),      PLATES)

CTRL_LEVELS <- c("Knockout", "Wild-type")
CTRL_PAL    <- setNames(c("#C62828", "#2E7D32"), CTRL_LEVELS)   # red KO, green WT
# Candidate palette: three hues chosen to NOT collide with the red/green controls
# (blue, amber, violet — all far from red and green in hue).
CAND_PAL_BASE <- c("#1565C0", "#E69F00", "#6A1B9A")
# Per-plate trajectory panel. `top_cand` is the named-by-plate list of chosen
# drugs; factored out so all three trajectory figures render from ONE code path.
traj_panel <- function(p, top_cand) {
  cand <- top_cand[[p]]
  d <- curve[bio_plate == p & Condition %in% c(POS_CONTROL, NEG_CONTROL, cand)]
  if (nrow(d) == 0) return(grid::nullGrob())
  d <- copy(d)
  d[, kind := fifelse(Condition %in% c(POS_CONTROL, NEG_CONTROL), "control", "candidate")]
  d[, ctrl := fifelse(Condition == POS_CONTROL, CTRL_LEVELS[1],
              fifelse(Condition == NEG_CONTROL, CTRL_LEVELS[2], NA_character_))]
  d[, drug := fifelse(kind == "candidate", sub("KO\\+", "", Condition), NA_character_)]
  d[, ctrl := factor(ctrl, levels = CTRL_LEVELS)]
  setorder(d, Condition, time_sec)
  d[, mv_s := frollmean(mv, 5, align = "center", fill = NA), by = Condition]
  d[is.na(mv_s), mv_s := mv]
  cand_lev <- sub("KO\\+", "", cand)
  d[, drug := factor(drug, levels = cand_lev)]
  cand_pal <- setNames(CAND_PAL_BASE[seq_along(cand)], cand_lev)
  # ggnewscale gives TWO real colour scales -> two separately-titled legends:
  # "Controls" (red/green) and "Top movers" (blue/amber/violet). Stacked
  # vertically with each group's title on the LEFT, keys laid out horizontally.
  ggplot(mapping = aes(time_sec, mv_s, group = Condition)) +
    geom_line(data = d[kind == "control"], aes(colour = ctrl), linewidth = 1.6) +
    scale_colour_manual(name = "Controls", values = CTRL_PAL, breaks = CTRL_LEVELS,
                        guide = guide_legend(order = 1, direction = "horizontal",
                          title.position = "left", title.hjust = 1, nrow = 1,
                          override.aes = list(linewidth = 1.6))) +
    ggnewscale::new_scale_colour() +
    geom_line(data = d[kind == "candidate"], aes(colour = drug),
              linewidth = 0.9, alpha = 0.9) +
    scale_colour_manual(name = "Top movers", values = cand_pal, breaks = cand_lev,
                        guide = guide_legend(order = 2, direction = "horizontal",
                          title.position = "left", title.hjust = 1, nrow = 1,
                          override.aes = list(linewidth = 1.3))) +
    labs(title = panel_title(p),
         x = "Time (s)", y = "Mean movement (pixels)") +
    theme_pub() +
    # Two legends stacked vertically (legend.box = "vertical"), each group's row
    # of keys CENTERED horizontally within the panel; titles on the left of their
    # row. Override the rotated colourbar title from theme_pub.
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.box.just = "center", legend.justification = "center",
          legend.spacing.y = unit(2, "pt"),
          legend.title = element_text(face = "bold", size = rel(0.8),
                                      angle = 0, hjust = 1, vjust = 0.5),
          legend.key.width = unit(1.6, "lines"),
          legend.text = element_text(size = rel(0.78)))
}

# Trajectory grids are 3-wide (14 plates → 5 rows). Taller than the scatters
# because each panel carries a two-row bottom legend.
TRAJ_W <- 16; TRAJ_H <- 24

# 3. closeness-to-WT selection — the original rescue view (landed nearest WT).
p3 <- assemble(lapply(PLATES, function(p) traj_panel(p, top_closeness)),
  "Movement over time: top rescue candidates vs the knockout and wild-type baselines",
  "Top movers = closest to wild-type. Red = knockout, green = wild-type; 5 s rolling mean.",
  ncol = 3)
save_fig("3_trajectories.png", p3, TRAJ_W, TRAJ_H)

# 3b. max-effect selection — identical layout/style, candidates are the strongest
# movers off the KO baseline (largest movement effect vs KO), not closeness.
p3b <- assemble(lapply(PLATES, function(p) traj_panel(p, top_max_effect)),
  "Movement over time: strongest movers off the disease baseline vs KO and wild-type",
  "Top movers = strongest off the knockout baseline. Red = knockout, green = wild-type; 5 s rolling mean.",
  ncol = 3)
save_fig("3b_trajectories_max_effect.png", p3b, TRAJ_W, TRAJ_H)

# 3c. equivalence selection — identical layout/style, candidates are the compounds
# of highest posterior equivalence to WT (the model's similarity-to-WT score),
# mirroring the 1e scatter.
p3c <- assemble(lapply(PLATES, function(p) traj_panel(p, top_equiv)),
  "Movement over time: compounds most equivalent to wild-type vs KO and wild-type",
  "Top movers = most likely equivalent to wild-type. Red = knockout, green = wild-type; 5 s rolling mean.",
  ncol = 3)
save_fig("3c_trajectories_equivalence.png", p3c, TRAJ_W, TRAJ_H)

# =============================================================================
# 4. CHECK A — model-free robust Glass's delta vs the BAM movement effect.
#    Points = the augmented per-plate table; Spearman r per panel from the
#    check_A CSV (fall back to computing it if the CSV is absent).
#    (Check B and the limma volcano were removed with the limma engine — on this
#    underpowered screen limma only re-confirmed the BAM, adding no signal.)
# =============================================================================
chkA_path <- file.path(PROD_DIR, "check_A_glass_vs_bam.csv")
chkA <- if (file.exists(chkA_path)) fread(chkA_path) else NULL
.spear <- function(x, y) suppressWarnings(cor(x, y, method = "spearman",
                                              use = "complete.obs"))
# Per-panel Spearman r as a corner annotation (Glass is in SD units, BAM in
# log2FC, so a rank check is the right comparison — points won't lie on a 1:1
# line and shouldn't be expected to).
spear_lab <- rbindlist(lapply(PLATES, function(p) {
  d <- aug[Plate == p]
  r <- if (!is.null(chkA) && nrow(chkA[Plate == p])) chkA[Plate == p, spearman[1]]
       else .spear(d$DrugEffect_vsKO_est_log2, d$Glass_vsKO_log)
  data.table(Plate = p, lab = sprintf("Spearman %.2f", r))
}))
spear_lab[, Plate := FP(Plate)]

pA_plot <- ggplot(aug, aes(DrugEffect_vsKO_est_log2, Glass_vsKO_log)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#1565C0", fill = "#1565C0",
              linewidth = 0.7, alpha = 0.18, formula = y ~ x) +
  geom_point(size = 1.7, alpha = 0.65, colour = "#1565C0", stroke = 0) +
  geom_text(data = spear_lab, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.6, size = 3, fontface = "bold", colour = "grey25") +
  facet_wrap(~ Plate, ncol = 3, scales = "free") +
  labs(x = expression("BAM movement effect  ("*log[2]*"FC vs KO)"),
       y = "Robust Glass's delta vs KO") +
  theme_pub()
pA <- assemble(list(pA_plot),
  "Check A: model-free Glass's delta vs the BAM movement effect",
  "Each point is a compound. High Spearman = model-free and BAM rankings agree.")
save_fig("4_checkA_glass_vs_bam.png", pA, FACET_W, FACET_H)

cat("Saved 10 figures to ", OUT, "\n", sep = "")

# =============================================================================
# POWER DIAGNOSTIC — why difference-test significance (a volcano) is unusable
# =============================================================================
pd_rows <- rbindlist(lapply(PLATES, function(p) {
  rs <- full[Plate == p]; n <- nrow(rs)
  sigKO <- sum(rs$DrugEffect_vsKO_adj_p_FDR < 0.05, na.rm = TRUE)
  sigWT <- sum(rs$Rescue_vsWT_adj_p_FDR    < 0.05, na.rm = TRUE)
  sigKO_nsWT <- sum(rs$DrugEffect_vsKO_adj_p_FDR < 0.05 &
                    !(rs$Rescue_vsWT_adj_p_FDR < 0.05), na.rm = TRUE)
  rescue_dir <- sum(rs$DrugEffect_vsKO_adj_p_FDR < 0.05 &
                    !(rs$Rescue_vsWT_adj_p_FDR < 0.05) &
                    rs$DrugEffect_vsKO_est_log2 > 0, na.rm = TRUE)
  sigKO_suppr <- sum(rs$DrugEffect_vsKO_adj_p_FDR < 0.05 &
                     rs$DrugEffect_vsKO_est_log2 < 0, na.rm = TRUE)
  data.table(Plate = p, n_drugs = n, sig_vs_KO = sigKO, ns_vs_KO = n - sigKO,
             sig_vs_WT = sigWT, ns_vs_WT = n - sigWT,
             sig_KO_and_ns_WT = sigKO_nsWT, rescue_direction = rescue_dir,
             sig_KO_suppressive = sigKO_suppr)
}))
pd_tot <- pd_rows[, lapply(.SD, sum), .SDcols = setdiff(names(pd_rows), "Plate")]
pd_tot[, Plate := "ALL"]
pd_all <- rbind(pd_rows, pd_tot, fill = TRUE)
setcolorder(pd_all, c("Plate", setdiff(names(pd_all), "Plate")))
fwrite(pd_all, file.path(PROD_DIR, "power_diagnostic_volcano.csv"))

cat("\n", strrep("=", 78), "\n",
    "POWER DIAGNOSTIC: why a volcano (difference-test significance) is not usable\n",
    strrep("=", 78), "\n", sep = "")
print(as.data.frame(pd_all), row.names = FALSE)
cat("-> If rescue_direction == 0, every sig-vs-KO hit is suppressive;",
    "significance tracks effect MAGNITUDE, not rescue.\n")
cat("Saved: ", file.path(PROD_DIR, "power_diagnostic_volcano.csv"), "\n", sep = "")
