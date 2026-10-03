#!/usr/bin/env Rscript
# =============================================================================
# mvp_ms_sim_figures.R -- manuscript figures for the simulation section (SS-Clines, 600
# replicates). Rebuilt 2026-10-03 (Phase 5 of the SS-Clines re-analysis) on the offset13 panels:
# top 0.25 % of SNPs per method and per predictor (RDA, one multivariate p for both predictors,
# 0.5 %), RDA without structure correction as the third method, 5 kb agreement window. The
# 2026-09-25 set (journal-16 panels, offset12) is archived read-only in FIG_OUT/_offset12/; git
# history holds the script that drew it.
#
# PLOTTING ONLY. Every number comes from tables already on disk (compute once, reuse):
#   GALLERY_DIR  mvp_j16_gallery.R run on offset13 (work/p5_figures.sh): detection_per_seed,
#                offset_delta_per_seed, detection_stats, C_delta_stats, D5_offset_stats
#   RANK_DIR     mvp_panel_rank.R: rank_long, rank_summary (within-replicate rank, per engine);
#                mvp_rank_reliability.R: reliability_{split_half,engines,variance} (numbers only)
#   SIZE_DIR     mvp_oracle_stats.R with MVP_SIZE_PANELS=1 on block b1: delta_per_seed (numbers only)
#   OFFSET_DIR   offset13 (name under benchmarks/mvp_eval): panel_pr_recomputed, snp_sets_summary
#   RDA_FIG_DIR  mvp_ms_rda_correction.R: S_rda_correction.{svg,png} (copied as is), numbers.tsv
#                (caption facts) and detection_table.tsv -- the remeasure harness's medians, an
#                independent tie for the panel builder's counts plotted here
# Medians and counts re-derived here are regression ties or caption facts only. NO INPUT HAS A
# DEFAULT: older arms (offset12_ssclines_pooled, figures_ssclines_j16_gallery) still exist on disk
# and a forgotten variable would give a complete, successful, wrong figure set.
#
#   main  A1_detection_plane       TP-vs-FP plane: 50 % bags + medians, own x-axis per panel
#         R1_rank_distribution     rank of the six panels inside each replicate, per offset engine
#                                  (main offset figure since the Phase 4b gate, 2026-10-03)
#   supp  C2_offset_distribution   absolute offset accuracy per panel + a separate causal-loci track
#         R2_rank_by_architecture  R1 split by genic level
#         D5_detection_cells       A1 view, genic level x (pleiotropy x selection regime)
#         D5_offset_cells          C2 view, same 12 cells
#         D6_detection_demography  A1 view, the five demography blocks
#         S_rda_correction         why RDA enters the panels WITHOUT structure correction (copied)
# Dropped at the 2026-10-03 review (user): C5_offset_reach (reach tile) and S_size_curve (its
# numbers stay in numbers.tsv). RDA WITH structure correction as an OFFSET ENGINE is in no figure,
# caption or table here (Phase 4b gate); as a GEA method it appears only in S_rda_correction.
#
# REVIEW FIXES, 2026-10-03 (user):
#   * Detection planes draw a BAG per rule -- the convex hull of the rule's median and the half
#     of the replicates closest to it -- instead of 2-D KDE contours. The KDE contours split into
#     disjoint pieces (A1: 15 of 36 contours, 20 pieces without the median; D5 cells: 124
#     orphan pieces): arcs clipped at the density grid's edge plus secondary modes such as the
#     FP = 0 heap. A wider bandwidth made it worse (work/p5_diag/contour_pieces.R). A bag is one
#     closed region that contains its median by construction; the script stops if one does not.
#   * Count axes: ticks EVENLY SPACED on the square-root scale (c * k^2), one step c per panel.
#   * Panel order of every per-panel figure is DATA-DRIVEN, not hand-set: fewest replicates in
#     which the panel ranks last (median of the three engines), ties by the last two ranks.
#   * Offset distributions show ABSOLUTE accuracy, with the causal loci as their own track below
#     a separator, instead of delta from the causal loci.
#
# LOOK: benchmarks/ms_style.R -- the shared manuscript style (ltc "minou" palette, 7 pt, print
# size, editable SVG via ms_save()). Rule colours: combined rules green/blue (2/3 teal, 1/3 sage,
# 3/3 navy), single methods red/orange (LFMM red, RDA amber, EMMAX orange), references grey (the
# causal-loci track). Every detection plane labels its medians by rule name; every offset
# distribution is the same raincloud; every rank distribution is the V4 view of mvp_panel_rank.R.
#
# Also writes captions.md, numbers.tsv (the numbers a caption or the Results text needs, long
# format, with source), README.md (text-only ranking-reliability result + how to reproduce every
# input of this report), index.html (served via work/journal/sim_figures on port 8099) and
# sim_figures_ms.zip (every SVG + PNG + captions.md + numbers.tsv).
#
#   FIG_OUT      default benchmarks/mvp_eval/figures_ssclines_ms
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT via benchmarks/mvp_arm.R (mandatory, 600)
#   Driver: work/p5_figures.sh ms
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggdist); library(ggrepel)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
need <- function(v) {
    x <- Sys.getenv(v, "")
    if (!nzchar(x)) stop("set ", v, " explicitly (no default: older arms exist on disk)")
    x
}
GAL  <- need("GALLERY_DIR")
RNK  <- need("RANK_DIR")
SIZ  <- need("SIZE_DIR")
ODIR <- file.path(EVAL, need("OFFSET_DIR"))
RDA  <- need("RDA_FIG_DIR")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_ms"))
JITTER_SEED <- 16L
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ vocabulary
# Labels and rule order as mvp_j16_gallery.R:57-70.
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
BLOCK_LAB <- c("N variable, m variable", "N cline north-south", "N equal, m constant",
               "N cline centre-edge", "N equal, m breaks")
RULES <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
# Hue = group (combined rules green/blue, single methods red/orange); ms_style.R palette.
RULE_COL <- c("2/3 methods" = MINOU[["teal"]], "1/3 methods" = MINOU[["sage"]],
              "3/3 methods" = MINOU[["navy"]], "LFMM" = MINOU[["red"]],
              "RDA" = MINOU[["amber"]], "EMMAX" = MS_ORANGE)
# Offset engines as mvp_panel_rank.R writes them -> manuscript labels.
ENG_LAB <- c(gradientForest = "Gradient Forest", LFMM2 = "LFMM2", RDA = "RDA",
             `median of 3` = "median of 3")
DELTA_LAB <- "Δ accuracy (panel − causal loci)"
N_PANEL <- length(RULES)

TEXT_TONES <- c("#FFFFFF" = "#FFFFFF", setNames(MS_INK, MS_INK))
ZERO <- geom_vline(xintercept = 0, linetype = "dashed", colour = MS_REF, linewidth = 0.35)

FIGS <- list()                                  # stem -> list(role, w, h, source)
save_ms <- function(stem, p, w, h, role, source) {
    ms_save(file.path(OUT, stem), p, w, h)
    FIGS[[stem]] <<- list(role = role, w = w, h = h, source = source)
    message(sprintf("  %-26s %3.0f x %3.0f mm", stem, w, h))
}
NUM <- list()                                   # numbers.tsv, appended per figure
num <- function(dt, figure, source) {
    NUM[[length(NUM) + 1L]] <<- cbind(data.table(figure = figure), dt, source = source)
}

# ------------------------------------------------------------------ inputs
man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
rd <- function(dir, stem, ...) {
    f <- file.path(dir, paste0(stem, ".tsv"))
    if (!file.exists(f)) stop("MISSING: ", f)
    fread(f, ...)
}
fac <- function(D) {
    D[, arch := factor(arch, levels = ARCH_LABELS)]
    if ("regime" %in% names(D)) D[, regime := factor(regime, levels = c("equal-S", "unequal-S"))]
    if ("pleio"  %in% names(D)) D[, pleio  := factor(pleio,  levels = c("no pleiotropy", "pleiotropy"))]
    if ("block"  %in% names(D)) D[, block  := factor(block,  levels = BLOCK_LAB)]
    D[, label := factor(label, levels = RULES)]
    stopifnot(!anyNA(D$arch), !anyNA(D$label))
    D[]
}
DET <- fac(rd(GAL, "detection_per_seed", colClasses = c(seed = "character"))[label %in% RULES])
OFF <- fac(rd(GAL, "offset_delta_per_seed", colClasses = c(seed = "character"))[label %in% RULES])
DS  <- rd(GAL, "detection_stats")[label %in% RULES]
CS  <- rd(GAL, "C_delta_stats")[label %in% RULES]
D5S <- fac(rd(GAL, "D5_offset_stats")[label %in% RULES])
stopifnot(!anyNA(DET$block), !anyNA(OFF$regime), !anyNA(D5S$pleio))
# One wrapped panel per design cell (gene level first, so ncol = 4 gives rows = gene level and
# columns = pleiotropy x regime). facet_wrap rather than facet_grid because only facet_wrap can
# give every panel its own x-axis without ggh4x, and both D5 figures use it so they pair up.
CELL_LEVELS <- CJ(arch = ARCH_LABELS, pleio = c("no pleiotropy", "pleiotropy"),
                  regime = c("equal-S", "unequal-S"), sorted = FALSE)[
    , paste0(arch, "\n", pleio, ", ", regime)]
for (D in list(DET, OFF, D5S))
    D[, panel := factor(paste0(arch, "\n", pleio, ", ", regime), levels = CELL_LEVELS)]
stopifnot(!anyNA(DET$panel), !anyNA(OFF$panel), !anyNA(D5S$panel))

# Corpus guards: every seed is a primary SS-Clines seed; detection keeps empty panels (all
# 600 per rule); design cells hold 200 / 50 / 120 replicates. Offset drops empty panels only.
NREP <- nrow(PRIM)
stopifnot(all(DET$seed %in% PRIM$seed), all(OFF$seed %in% PRIM$seed),
          uniqueN(DET$seed) == NREP, uniqueN(OFF$seed) == NREP,
          DET[, .N, by = label][, all(N == nrow(PRIM))],
          DET[, .N, by = .(label, arch)][, all(N == nrow(PRIM) / 3)],
          DET[, .N, by = .(label, arch, pleio, regime)][, all(N == nrow(PRIM) / 12)],
          DET[, .N, by = .(label, block)][, all(N == nrow(PRIM) / 5)])
EMPTY <- DET[n == 0, .N, by = .(label = as.character(label))]
off_missing <- merge(CJ(seed = PRIM$seed, label = RULES), OFF[, .(seed, label = as.character(label), x = 1L)],
                     by = c("seed", "label"), all.x = TRUE)[is.na(x)]
stopifnot(nrow(off_missing) == sum(EMPTY$N),
          nrow(merge(off_missing, DET[n == 0, .(seed, label = as.character(label))],
                     by = c("seed", "label"))) == nrow(off_missing))
message(sprintf("empty panels: %s (no offset, ranked last)",
                if (nrow(EMPTY)) paste(sprintf("%s x %d", EMPTY$label, EMPTY$N), collapse = ", ") else "none"))

# Regression tie: medians re-derived from the per-seed tables equal the gallery's stats tables.
tie <- function(a, b, keys, what) {
    x <- merge(a, b, by = keys, all = TRUE)
    stopifnot(nrow(x) == nrow(a), !anyNA(x$m), !anyNA(x$median),
              max(abs(x$m - x$median)) < 1e-9, all(x$n.x == x$n.y))
    message(sprintf("tie %-24s %3d medians, max |diff| %.1e", what, nrow(x), max(abs(x$m - x$median))))
}
tie(OFF[, .(m = median(delta), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
    CS[arch != "pooled", .(label, arch, median, n)], c("label", "arch"), "C offset by arch")
tie(OFF[, .(m = median(delta), n = .N), by = .(label = as.character(label))],
    CS[arch == "pooled", .(label, median, n)], "label", "C offset pooled")
tie(OFF[, .(m = median(delta), n = .N), by = .(label, arch, pleio, regime)],
    D5S[, .(label, arch, pleio, regime, median, n)], c("label", "arch", "pleio", "regime"), "D5 offset")
tie(DET[, .(m = median(TP), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
    DS[arch != "pooled", .(label, arch, median = TP_median, n = n_seeds)], c("label", "arch"), "detection TP")
tie(DET[, .(m = median(FP), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
    DS[arch != "pooled", .(label, arch, median = FP_median, n = n_seeds)], c("label", "arch"), "detection FP")

# Independent tie: the panel builder's counts (offset13 -> gallery) against the remeasure
# harness's medians in the RDA brief (remeasure600/rdaunc_rda2x, top 0.25 %, RDA x 2, 5 kb).
# Two code paths, one operating point: this fails if either arm is not the one the captions name.
DTB <- fread(file.path(RDA, "detection_table.tsv"))
dmed <- function(D) D[, .(n_med = median(n), causal = median(n_causal), linked = median(n_linked),
                          background = median(FP), TP = median(TP), empty = sum(n == 0), n_rep = .N)]
DM <- rbind(DET[, dmed(.SD), by = .(rule = as.character(label), arch = as.character(arch))],
            DET[, dmed(.SD), by = .(rule = as.character(label))][, arch := "pooled"], use.names = TRUE)
x <- merge(DM, DTB[, .(rule, arch, b_n = n_median, b_causal = causal_median, b_linked = linked_median,
                       b_background = background_median, b_TP = TP_A1_median, b_empty = empty,
                       b_rep = n_rep)], by = c("rule", "arch"), all = TRUE)
stopifnot(nrow(x) == N_PANEL * 4L, !anyNA(x),
          max(abs(as.matrix(x[, .(n_med - b_n, causal - b_causal, linked - b_linked,
                                  background - b_background, TP - b_TP)]))) < 1e-9,
          all(x$empty == x$b_empty), all(x$n_rep == x$b_rep))
message(sprintf("tie detection vs RDA brief (harness)  %d rule x stratum rows, 5 medians each: exact", nrow(x)))

# =============================================================================
# Panel order -- one data-driven rule for every per-panel figure (user, 2026-10-03: a hand-set
# order with 2/3 on top reads as cherry-picking). Rule: share of replicates in which the panel
# ranks LAST among the six, under the median of the three offset engines, fewest first; ties
# by the share in the last two ranks. It is the stability criterion fixed at the Phase 4b gate
# ("rarely the panel you would regret"). The order is computed, never asserted.
# =============================================================================
RL <- rd(RNK, "rank_long", colClasses = c(seed = "character"))
RS <- rd(RNK, "rank_summary")
stopifnot(setequal(unique(RL$engine), names(ENG_LAB)), setequal(unique(RL$panel), RULES),
          all(RL$seed %in% PRIM$seed), uniqueN(RL$seed) == NREP,
          nrow(RL) == NREP * length(ENG_LAB) * N_PANEL, !anyDuplicated(RL[, .(seed, engine, panel)]))
ORD <- dcast(RS[stratum == "pooled" & engine == "median of 3" & metric %in% c("bottom1", "bottom2")],
             panel ~ metric, value.var = "value")
stopifnot(nrow(ORD) == N_PANEL, setequal(ORD$panel, RULES), !anyNA(ORD))
setorder(ORD, bottom1, bottom2)
PANEL_ORDER <- ORD$panel
message("panel order (fewest replicates ranked last, median of 3 engines): ",
        paste(sprintf("%s %.1f%%", ORD$panel, ORD$bottom1), collapse = " < "))

# =============================================================================
# Detection plane (A1, D5-det, D6-det)
# =============================================================================
# Ticks evenly spaced on the square-root scale: breaks c * k^2, k = 0, 1, 2, ..., with the
# smallest step c from a round set that gives at most 6 ticks over the panel's range. Called per
# panel (free x), so every panel gets its own even ticks.
SQRT_STEPS <- c(1, 2, 5, 10, 20, 25, 50, 100)
sqrt_even_breaks <- function(limits) {
    m <- max(limits, na.rm = TRUE)
    for (c in SQRT_STEPS) {
        n <- floor(sqrt(m / c)) + 1L
        if (n <= 6L) return(c * (seq_len(n) - 1L)^2)
    }
    c(0, m)
}
# Point inside a polygon, or on one of its edges (ray casting + edge test).
in_polygon <- function(px, py, x, y, tol = 1e-9) {
    n <- length(px); j <- n; inside <- FALSE
    for (i in seq_len(n)) {
        ex <- px[j] - px[i]; ey <- py[j] - py[i]
        if (abs((x - px[i]) * ey - (y - py[i]) * ex) < tol &&
            x >= min(px[i], px[j]) - tol && x <= max(px[i], px[j]) + tol &&
            y >= min(py[i], py[j]) - tol && y <= max(py[i], py[j]) + tol) return(TRUE)
        if (((py[i] > y) != (py[j] > y)) && (x < px[i] + (y - py[i]) * ex / ey)) inside <- !inside
        j <- i
    }
    inside
}
# BAG: the convex hull of the median and the half of the replicates closest to it, distance
# measured on the plotted (square-root) axes, each scaled by its interquartile range. The median
# is a hull point because a single-method panel has a near-fixed size (TP + FP ~ constant): its
# replicates lie along an anti-diagonal and the component-wise median falls just off that line
# (3 of 72 D5 cells, n = 50), outside a hull of the replicates alone. NULL when there is
# no spread to draw (fewer than 10 replicates, a zero IQR, or a collinear half): that rule is
# drawn by its median only, and the caption names it.
BAG_SHARE <- 0.5
bag_hull <- function(fp, tp) {
    x <- sqrt(fp); y <- sqrt(tp)
    sx <- IQR(x); sy <- IQR(y)
    if (length(x) < 10L || sx == 0 || sy == 0) return(NULL)
    mx <- sqrt(median(fp)); my <- sqrt(median(tp))
    k <- order(((x - mx) / sx)^2 + ((y - my) / sy)^2)[seq_len(ceiling(BAG_SHARE * length(x)))]
    hx <- c(x[k], mx); hy <- c(y[k], my)
    h <- chull(hx, hy)
    if (length(h) < 3L) return(NULL)
    list(FP = hx[h]^2, TP = hy[h]^2, has_median = in_polygon(hx[h], hy[h], mx, my))
}
bags_of <- function(D, keys) {
    B <- D[, { b <- bag_hull(FP, TP); if (is.null(b)) NULL else b }, by = keys]
    bad <- unique(B[has_median == FALSE, keys, with = FALSE])
    if (nrow(bad)) stop("bag without its median (would be an orphan region): ",
                        paste(do.call(paste, c(bad, sep = " / ")), collapse = "; "))
    B[, grp := do.call(paste, c(.SD, sep = "|")), .SDcols = keys][]
}

plane_ms <- function(D, facet) {
    keys <- c("label", facet$vars)
    M  <- D[, .(TP = median(TP), FP = median(FP)), by = keys]
    B  <- bags_of(D, keys)
    p <- ggplot(D, aes(FP, TP)) +
        geom_polygon(data = B, aes(group = grp, colour = label), fill = NA, linewidth = 0.35) +
        geom_point(data = M, aes(fill = label), shape = 21, size = 1.7, colour = "black", stroke = 0.3) +
        # Median labels are the rule names on solid tags in the rule colour; the text is white or
        # dark, whichever contrasts more with that colour (ms_text_on). The text tone rides on the
        # same colour scale as the bags, under its own two keys.
        geom_label_repel(data = M[, tone := ms_text_on(RULE_COL[as.character(label)])],
                         aes(label = label, fill = label, colour = tone), size = MS_LAB,
                         fontface = "bold", label.size = 0, label.padding = 0.12, box.padding = 0.3,
                         min.segment.length = 0, segment.size = 0.25, segment.colour = MS_INK,
                         seed = JITTER_SEED, show.legend = FALSE, max.overlaps = Inf)
    # Both count axes start at 0 (the scales train on bags + medians only).
    p + facet$layer + expand_limits(x = 0, y = 0) +
        scale_x_sqrt(breaks = sqrt_even_breaks) + scale_y_sqrt(breaks = sqrt_even_breaks) +
        scale_colour_manual(values = c(RULE_COL, TEXT_TONES), guide = "none") +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        labs(x = "False positives", y = "True positives") +
        theme_ms()
}
# Caption clause naming the rules drawn by their median only; empty when every rule has a bag.
median_only <- function(D, vars) {
    keys <- c("label", vars)
    has <- unique(bags_of(D, keys)[, keys, with = FALSE])[, bag := TRUE]
    x <- merge(unique(D[, keys, with = FALSE]), has, by = keys, all.x = TRUE)[is.na(bag), .N, by = label]
    x <- x[order(match(label, RULES))]
    if (!nrow(x)) return("")
    paste0(" Shown by the median only (no spread for a bag): ",
           paste(sprintf("%s (%d of %d panels)", x$label, x$N,
                         nrow(unique(D[, vars, with = FALSE]))), collapse = ", "), ".")
}
# Every panel of a detection plane: tick labels present and evenly spaced on the sqrt scale.
check_ticks <- function(p, what) {
    for (pp in ggplot_build(p)$layout$panel_params) for (ax in c("x", "y")) {
        b <- suppressWarnings(as.numeric(pp[[ax]]$get_labels()))
        b <- b[!is.na(b)]
        g <- diff(sqrt(b))
        if (length(b) < 2L || max(abs(g - g[1])) > 1e-9)
            stop(what, ": ", ax, " ticks not evenly spaced on the sqrt scale: ", paste(b, collapse = " "))
    }
    message(sprintf("%-24s every panel: x and y ticks evenly spaced on the sqrt scale", what))
}

message("=== detection")
# Every detection plane gives each panel its own false-positive axis (scales = "free_x").
pA1 <- plane_ms(DET, list(vars = "arch", layer = facet_grid(cols = vars(arch), scales = "free_x")))
for (pp in ggplot_build(pA1)$layout$panel_params)
    message("A1 panel axis labels  x: ", paste(pp$x$get_labels(), collapse = " "),
            "   y: ", paste(pp$y$get_labels(), collapse = " "))
check_ticks(pA1, "A1_detection_plane")
save_ms("A1_detection_plane", pA1, 180, 62, "main", "detection_per_seed.tsv")

pD5d <- plane_ms(DET, list(vars = "panel", layer = facet_wrap(vars(panel), ncol = 4, scales = "free_x")))
check_ticks(pD5d, "D5_detection_cells")
save_ms("D5_detection_cells", pD5d, 180, 165, "supp", "detection_per_seed.tsv")

pD6 <- plane_ms(DET, list(vars = "block", layer = facet_wrap(vars(block), nrow = 2, scales = "free_x")))
check_ticks(pD6, "D6_detection_demography")
save_ms("D6_detection_demography", pD6, 180, 110, "supp", "detection_per_seed.tsv")

num(melt(DS[, .(rule = label, stratum = arch, n = n_seeds, n_empty = as.numeric(n_empty), TP_median, TP_q25, TP_q75,
                FP_median, FP_q25, FP_q75, prec_cl_median, F1_median)],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "A1 / D5 / D6 detection", "figures_ssclines_gallery13/detection_stats.tsv")

# =============================================================================
# Offset accuracy (C2, D5-offset): absolute accuracy per panel + the causal loci as a track
# =============================================================================
message("=== offset accuracy")
ORACLE_LAB <- "causal loci"
ARCH_WRAP <- c("oligogenic" = "oligogenic", "moderately polygenic" = "moderately\npolygenic",
               "highly polygenic" = "highly\npolygenic")
# accuracy - oracle == delta on every row: the two columns are the gallery's own decomposition.
stopifnot(max(abs(OFF$accuracy - OFF$oracle - OFF$delta)) < 1e-12)
ORC <- unique(OFF[, .(seed, arch, pleio, regime, panel, accuracy = oracle)])
stopifnot(nrow(ORC) == NREP, uniqueN(ORC$seed) == NREP)      # one causal-loci value per replicate
ACC <- rbind(OFF[, .(seed, arch, pleio, regime, panel, label = as.character(label), accuracy)],
             ORC[, label := ORACLE_LAB], use.names = TRUE)
# Rows top to bottom: the six panels in PANEL_ORDER, then -- below a separator -- the causal loci.
Y_AT  <- setNames(c(rev(seq_len(N_PANEL)) + 1, 0.5), c(PANEL_ORDER, ORACLE_LAB))
SEP_Y <- 1.3
ACC_COL <- c(RULE_COL, setNames(MS_REF, ORACLE_LAB))
acc_cloud <- function(D, facet_layer) {
    RC <- copy(D)[, ypos := Y_AT[label]]
    stopifnot(!anyNA(RC$ypos))
    ggplot(RC, aes(accuracy, ypos)) +
        geom_hline(yintercept = SEP_Y, colour = MS_REF, linewidth = 0.3) +
        stat_halfeye(aes(fill = label, group = label), orientation = "horizontal", adjust = 0.8,
                     height = 0.55, justification = -0.2, .width = 0, point_colour = NA,
                     slab_alpha = 0.75) +
        geom_boxplot(aes(group = label), orientation = "y", width = 0.14, outlier.shape = NA,
                     fill = "white", colour = MS_INK, linewidth = 0.25) +
        geom_point(aes(y = ypos - 0.25, colour = label),
                   position = position_jitter(height = 0.08, width = 0, seed = JITTER_SEED),
                   size = 0.6, alpha = 0.6, stroke = 0) +
        facet_layer +
        scale_y_continuous(breaks = unname(Y_AT), labels = names(Y_AT)) +
        scale_fill_manual(values = ACC_COL, guide = "none") +
        scale_colour_manual(values = ACC_COL, guide = "none") +
        labs(x = "Offset accuracy (−Kendall τ)", y = NULL) + theme_ms()
}
save_ms("C2_offset_distribution",
        acc_cloud(ACC, facet_grid(cols = vars(arch), labeller = labeller(arch = ARCH_WRAP))),
        114, 66, "supp", "offset_delta_per_seed.tsv (accuracy, oracle)")
# D5 offset: the same view per design cell; the x-axis stays SHARED across cells.
save_ms("D5_offset_cells", acc_cloud(ACC, facet_wrap(vars(panel), ncol = 4)),
        180, 180, "supp", "offset_delta_per_seed.tsv (accuracy, oracle)")

acc_stats <- function(D) D[, .(n = .N, median = median(accuracy), q25 = quantile(accuracy, .25),
                               q75 = quantile(accuracy, .75))]
ACS <- rbind(ACC[, acc_stats(.SD), by = .(rule = label, stratum = as.character(arch))],
             ACC[, acc_stats(.SD), by = .(rule = label)][, stratum := "pooled"], use.names = TRUE)
num(melt(ACS, id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "C2 / D5 offset accuracy (median of 3 engines)", "figures_ssclines_gallery13/offset_delta_per_seed.tsv")
# Delta from the causal loci: no figure any more, the Results text still quotes it.
num(melt(CS[, .(rule = label, stratum = arch, n, median, ci_lo, ci_hi, q25, q75, win_pct,
                reach_pct, p_wilcox_two_sided)],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "offset delta vs causal loci (text only, no figure)", "figures_ssclines_gallery13/C_delta_stats.tsv")

# =============================================================================
# Within-replicate rank (R1 main, R2 supp) -- the V4 view of mvp_panel_rank.R
# =============================================================================
message("=== rank")
# Empty panels in the rank table are exactly the empty detection panels, one rank per engine.
stopifnot(setequal(RL[empty == TRUE, unique(paste(seed, panel))],
                   DET[n == 0, paste(seed, label)]))
# Tie: mean rank re-derived here equals mvp_panel_rank.R's pooled summary, every engine x panel.
mr <- merge(RL[, .(m = mean(rank)), by = .(engine, panel)],
            RS[stratum == "pooled" & metric == "mean_rank", .(engine, panel, value)],
            by = c("engine", "panel"), all = TRUE)
stopifnot(nrow(mr) == length(ENG_LAB) * N_PANEL, !anyNA(mr), max(abs(mr$m - mr$value)) < 1e-9)
message(sprintf("tie %-24s %3d mean ranks, max |diff| %.1e", "rank summary", nrow(mr), max(abs(mr$m - mr$value))))
RL[, key := fifelse(empty, -Inf, accuracy)]
N_TIED <- RL[, .(tied = any(duplicated(key))), by = .(seed, engine)][, sum(tied)]
N_UNITS <- uniqueN(RL[, .(seed, engine)])
RL[, `:=`(arch = factor(arch_lab, levels = ARCH_LABELS),
          engine_lab = factor(ENG_LAB[engine], levels = ENG_LAB),
          # a tied pair shares the average rank; drawn at the nearest whole rank, as V4
          r = pmin(N_PANEL, pmax(1, round(rank))))]
stopifnot(!anyNA(RL$arch), !anyNA(RL$engine_lab))

RANK_COL <- setNames(colorRampPalette(c(MINOU[["red"]], "#F2F2F2", MINOU[["teal"]]))(N_PANEL),
                     N_PANEL:1)
rank_dist <- function(D, extra = NULL, facet_layer) {
    R <- D[, .(n = .N), by = c("engine_lab", "panel", "r", extra)]
    R[, pct := 100 * n / sum(n), by = c("engine_lab", "panel", extra)]
    R[, r := factor(r, levels = N_PANEL:1)]
    R[, panel := factor(panel, levels = rev(PANEL_ORDER))]
    ggplot(R, aes(pct, panel, fill = r)) +
        geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
        facet_layer +
        scale_fill_manual(values = RANK_COL, name = "rank (1 = most accurate)",
                          breaks = as.character(1:N_PANEL)) +
        scale_x_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
        labs(x = "% replicates at each rank", y = NULL) + theme_ms()
}
save_ms("R1_rank_distribution", rank_dist(RL, NULL, facet_wrap(~ engine_lab, nrow = 1)),
        180, 66, "main", "rank_long.tsv")
save_ms("R2_rank_by_architecture",
        rank_dist(RL, "arch", facet_grid(rows = vars(arch), cols = vars(engine_lab))),
        180, 120, "supp", "rank_long.tsv")

num(RS[stratum %in% c("pooled", ARCH_LABELS) &
       metric %in% c("top1", "top2", "bottom1", "bottom2", "mean_rank"),
       .(rule = panel, stratum, n, engine = ENG_LAB[engine], metric, value, lo, hi, expected)],
    "R1 / R2 within-replicate rank", "figures_ssclines_offset13/rank/rank_summary.tsv")
num(data.table(metric = c("units_with_tie", "units"), value = c(N_TIED, N_UNITS)),
    "R1 / R2 within-replicate rank", "figures_ssclines_offset13/rank/rank_long.tsv")
# Ranking reliability (mvp_rank_reliability.R): text only, no figure (user, 2026-10-03).
RSH <- rd(RNK, "reliability_split_half"); REN <- rd(RNK, "reliability_engines"); RVP <- rd(RNK, "reliability_variance")
stopifnot(nrow(RSH) == 4L, nrow(REN) == 4L, nrow(RVP) == 3L)
num(melt(RSH[, .(engine, n = n_rep, r_half, r_half_lo, r_half_hi, r_full, worst_same_pct, worst_same_lo,
                 worst_same_hi, best_same_pct, chance_pct, n_splits = as.numeric(n_splits))],
         id.vars = c("engine", "n"), variable.name = "metric", value.name = "value")[, stratum := "pooled"],
    "text: ranking reliability, split-half (Methods/Results)",
    "figures_ssclines_offset13/rank/reliability_split_half.tsv")
num(rbind(melt(REN[test == "kendall_w", .(engine = engine_a, n = n_rep, kendall_w_mean = value, kendall_w_median = median,
                                          null_mean, null_lo, null_hi, p_perm, n_perm = as.numeric(n_perm))],
               id.vars = c("engine", "n"), variable.name = "metric", value.name = "value"),
          melt(REN[test == "pair", .(engine = paste(engine_a, "vs", engine_b), n = n_rep, r_centred = value,
                                     worst_same_pct, best_same_pct, chance_pct)],
               id.vars = c("engine", "n"), variable.name = "metric", value.name = "value"))[, stratum := "pooled"],
    "text: ranking reliability, engine agreement (Methods/Results)",
    "figures_ssclines_offset13/rank/reliability_engines.tsv")
num(rbind(RVP[, .(rule = source, metric = "variance_share_pct", value = share_pct, n = n_rep)],
          RVP[1, .(rule = "all panels", metric = c("panel_median_range", "within_replicate_range_median"),
                   value = c(panel_median_range, within_rep_range_median), n = n_rep)])[
              , `:=`(engine = "median of 3", stratum = "pooled")],
    "text: ranking reliability, variance partition (Methods/Results)",
    "figures_ssclines_offset13/rank/reliability_variance.tsv")
num(data.table(rule = PANEL_ORDER, metric = "panel_order", value = seq_along(PANEL_ORDER),
               stratum = "pooled", engine = "median of 3"),
    "panel order (all per-panel figures)", "figures_ssclines_offset13/rank/rank_summary.tsv (bottom1, bottom2)")

# =============================================================================
# Size curve (block b1): numbers only (figure dropped at the 2026-10-03 review)
# =============================================================================
message("=== size curve (numbers only)")
RUNG <- c("1/3 methods, top 0.1%" = 0.1, "1/3 methods" = 0.25, "1/3 methods, top 0.5%" = 0.5,
          "1/3 methods, top 1%" = 1, "1/3 methods, top 2%" = 2)
SET_RUNG <- c(union_top0.1pct = 0.1, union = 0.25, union_top0.5pct = 0.5,
              union_top1pct = 1, union_top2pct = 2)
B1_TAG <- "ssclines_nvar_mvar"
SB1 <- PRIM[added == B1_TAG, seed]
SZ <- rd(SIZ, "delta_per_seed", colClasses = c(seed = "character"))[label %in% names(RUNG)]
SZ[, share := RUNG[label]]
stopifnot(length(SB1) == NREP / 5, setequal(unique(SZ$seed), SB1),
          SZ[, .N, by = share][, .N == length(RUNG) && all(N == length(SB1))])
# Tie: the 0.25 % rung IS the main 1/3 panel; its deltas on b1 equal the gallery's.
tz <- merge(SZ[share == 0.25, .(seed, d = delta)],
            OFF[label == "1/3 methods" & seed %in% SB1, .(seed, delta)], by = "seed", all = TRUE)
stopifnot(nrow(tz) == length(SB1), !anyNA(tz), max(abs(tz$d - tz$delta)) < 1e-9)
message(sprintf("tie %-24s %3d deltas, max |diff| %.1e", "size 0.25 % vs gallery", nrow(tz), max(abs(tz$d - tz$delta))))
PPR <- fread(file.path(ODIR, "panel_pr_recomputed.tsv"), colClasses = c(seed = "character"))
SZN <- PPR[seed %in% SB1 & set %in% names(SET_RUNG),
           .(n_snps = median(n), n_rep = .N), by = set][, share := SET_RUNG[set]]
stopifnot(nrow(SZN) == length(SET_RUNG), all(SZN$n_rep == length(SB1)))
fw <- dcast(SZ, seed ~ share, value.var = "delta")
stopifnot(!anyNA(fw))
FR <- friedman.test(as.matrix(fw[, -1]))
SZS <- merge(SZ[, .(n = .N, median = median(delta), q25 = quantile(delta, .25),
                    q75 = quantile(delta, .75), below_pct = 100 * mean(delta < 0)), by = share],
             SZN[, .(share, n_snps)], by = "share")
setorder(SZS, share)
message(sprintf("size curve: Friedman chi2 = %.2f, df = %d, p = %.3f; medians %s", FR$statistic,
                FR$parameter, FR$p.value, paste(sprintf("%.3f", SZS$median), collapse = " / ")))
num(rbind(melt(SZS[, .(stratum = sprintf("top %s%%", share), n, n_snps, median, q25, q75, below_pct)],
               id.vars = c("stratum", "n"), variable.name = "metric", value.name = "value"),
          data.table(stratum = "all five shares", n = length(SB1),
                     metric = c("friedman_chi2", "friedman_df", "friedman_p"),
                     value = c(unname(FR$statistic), unname(FR$parameter), FR$p.value)))[
              , rule := "1/3 methods"],
    "size curve, block N variable, m variable (text only, no figure)",
    "figures_ssclines_offset13/size_b1/delta_per_seed.tsv + offset13/panel_pr_recomputed.tsv")

# =============================================================================
# S_rda_correction: copied as built by mvp_ms_rda_correction.R (its own ties live there);
# its caption argues for the uncorrected RDA, from that script's numbers.tsv.
# =============================================================================
for (ext in c("svg", "png")) {
    src <- file.path(RDA, paste0("S_rda_correction.", ext))
    if (!file.exists(src)) stop("MISSING: ", src)
    stopifnot(file.copy(src, file.path(OUT, basename(src)), overwrite = TRUE))
}
FIGS[["S_rda_correction"]] <- list(role = "supp", w = 180, h = 60,
                                   source = "figures_ssclines_rda/ (mvp_ms_rda_correction.R)")
message(sprintf("  %-26s %3.0f x %3.0f mm (copied)", "S_rda_correction", 180, 60))
RNUM <- fread(file.path(RDA, "numbers.tsv"), colClasses = "character")
rn <- function(k) {
    v <- RNUM[key == k, shown]
    if (length(v) != 1L || !nzchar(v)) stop("RDA numbers.tsv: key ", k, " missing or not unique")
    v
}
RDA_FACTS <- list(wins = rn("wins_unc_pooled"), n = rn("n_pooled"),
                  unc = rn("aucpr_RDA_uncorrected_pooled_median"),
                  cor = rn("aucpr_RDA_corrected_pooled_median"), ratio = rn("ratio_med_pooled"),
                  kept = rn("kept_median"), vif10 = rn("vif_ge10"))

# =============================================================================
# Operating-point justification: share of causal loci among tested SNPs (offset13 'all' panel =
# every tested SNP). 0.25 % per method and predictor was set to the corpus median (Phase 4b gate).
# =============================================================================
ARCH_OF <- c("oliogenic" = "oligogenic", "mod-polygenic" = "moderately polygenic",
             "highly-polygenic" = "highly polygenic")              # corpus ships the typo
SSA <- fread(file.path(ODIR, "snp_sets_summary.tsv"), colClasses = c(seed = "character"))[set == "all"]
SSA <- merge(SSA[, .(seed, n_snps, n_causal)], PRIM[, .(seed, arch = ARCH_OF[arch_level])], by = "seed")
stopifnot(nrow(SSA) == NREP, !anyNA(SSA$arch), all(SSA$n_snps > 0))
SSA[, cshare := 100 * n_causal / n_snps]
cs_stats <- function(D) D[, .(n = .N, median = median(cshare), q25 = quantile(cshare, .25),
                              q75 = quantile(cshare, .75), geomean = exp(mean(log(cshare))),
                              tested_median = median(n_snps), causal_median = median(n_causal))]
CSH <- rbind(SSA[, cs_stats(.SD), by = .(stratum = arch)], cs_stats(SSA)[, stratum := "pooled"], use.names = TRUE)
CS_MED <- CSH[stratum == "pooled", median]
message(sprintf("causal share of tested SNPs: median %.3f %% (geomean %.2f %%)", CS_MED,
                CSH[stratum == "pooled", geomean]))
num(melt(CSH, id.vars = c("stratum", "n"), variable.name = "metric", value.name = "value")[
        , rule := "causal loci, % of tested SNPs"],
    "operating point (Methods)", "offset13/snp_sets_summary.tsv (set all) + mvp_seeds.tsv")

# =============================================================================
# Captions (all explanation lives here, not on the figures)
# =============================================================================
nP <- function(rule) as.integer(CS[label == rule & arch == "pooled", n])
EMPTY_CLAUSE <- if (nrow(EMPTY)) paste0(
    paste(sprintf("The %s panel is empty in %d replicate%s", EMPTY$label, EMPTY$N,
                  ifelse(EMPTY$N == 1, "", "s")), collapse = "; "), ".") else ""
OPERATING_DEF <- paste0(
    "Operating point, fixed for all replicates: each method keeps the 0.25% of tested SNPs with the ",
    "lowest p-value for each of the two environmental predictors (LFMM and EMMAX test the predictors ",
    "separately and their two sets are pooled; RDA, without structure correction, tests both in one ",
    "multivariate p and keeps 0.5%), so each method's panel holds about 0.5% of the tested SNPs; ",
    sprintf("0.25%% is the median share of causal loci among tested SNPs over the %d replicates (%.2f%%). ", NREP, CS_MED),
    "1/3 methods: every SNP called by at least one method; 2/3 and 3/3 methods: SNPs within 5 kb of ",
    "calls by at least two / all three methods (each method keeps its own SNPs, so a 3/3 panel can be ",
    "larger than any single-method panel).")
DETECTION_DEF <- paste0(
    "True positives are selected loci that are causal or linked-neutral (linkage groups 1–10); ",
    "false positives are selected background-neutral loci (linkage groups 11–20). ",
    "Both axes use a square-root scale with evenly spaced ticks; each panel has its own false-positive axis. ",
    "Outlines (bags) are the convex hull of each rule's median and the half of its replicates closest to it (distance on the ",
    "square-root axes, each scaled by its interquartile range); labelled points are the median (true and false positives medianed separately). ",
    OPERATING_DEF, " ", EMPTY_CLAUSE, if (nzchar(EMPTY_CLAUSE)) " " else "",
    "The causal loci are not shown, because by construction they carry no linked hits.")
ACC_DEF <- paste0(
    "Accuracy is −Kendall τ between the predicted genetic offset and the fitness of the 100 ",
    "source populations in a common garden, median over the 100 landscape gardens, for each of three ",
    "offset engines (Gradient Forest, LFMM2 geometric offset, RDA without structure correction).")
ORDER_DEF <- paste0(
    "Panels are ordered by the share of replicates in which they rank last under the median of the ",
    "three engines, fewest first (ties broken by the share in the last two ranks).")
RANK_DEF <- paste0(
    "Within each replicate the six marker panels are ranked by offset accuracy (1 = most accurate), ",
    "separately for each offset engine and for the median of the three. ", ACC_DEF,
    " The causal loci are a constant reference within a replicate, so ranks on accuracy equal ranks on accuracy relative to the causal loci. ",
    sprintf("Under random ordering every panel would sit at each rank in %.1f%% of replicates. ", 100 / N_PANEL),
    sprintf("Tied panels share the average rank, drawn at the nearest whole rank (%d of %s replicate × engine units carry a tie). ",
            N_TIED, format(N_UNITS, big.mark = ",")),
    if (nrow(EMPTY)) "An empty panel has no offset and is ranked last. " else "",
    ORDER_DEF, " Panels as in Main a.")
ACC3_DEF <- paste0(
    "Accuracy is −Kendall τ between the predicted genetic offset and the fitness of the 100 source ",
    "populations in a common garden, median over the 100 landscape gardens and then over three offset ",
    "engines (Gradient Forest, LFMM2 geometric offset, RDA without structure correction). ",
    "Half-violins: density; boxes: median and interquartile range; dots: replicates. ",
    "The grey track below the line is the same accuracy computed from the causal loci of each replicate (reference). ",
    if (nrow(EMPTY)) sprintf("Empty panels have no offset and are omitted (%s). ",
                             paste(sprintf("%s: n = %d of %d", EMPTY$label,
                                           vapply(EMPTY$label, nP, 1L), NREP), collapse = "; ")) else "",
    "Panels ordered as in Main b.")
CAP <- list(
    A1_detection_plane = paste0(
        sprintf("Detection of adaptive loci by six marker-selection rules, by genetic architecture (n = %d replicates per level). ", NREP / 3),
        DETECTION_DEF, median_only(DET, "arch")),
    R1_rank_distribution = paste0(
        sprintf("Rank of each marker panel among the six within a replicate, per offset engine (n = %d replicates). ", NREP),
        "Bars: share of replicates at each rank. ", RANK_DEF),
    C2_offset_distribution = paste0(
        sprintf("Genetic-offset accuracy per replicate for each marker panel and for the causal loci, by genetic architecture (n = %d per level). ", NREP / 3),
        ACC3_DEF),
    R2_rank_by_architecture = paste0(
        sprintf("Rank distribution of Main b split by genetic architecture (rows; n = %d replicates per level) and offset engine (columns). ", NREP / 3),
        RANK_DEF),
    D5_detection_cells = paste0(
        sprintf("Detection plane (as in Main a) for every design cell: genetic architecture (rows) × pleiotropy × selection regime (columns), n = %d replicates per cell. ", NREP / 12),
        "equal-S: both environmental axes under equal selection; unequal-S: selection on one axis 8-fold weaker. ",
        DETECTION_DEF, median_only(DET, "panel")),
    D5_offset_cells = paste0(
        sprintf("Genetic-offset accuracy for every design cell, drawn as Supplementary 1 (layout as Supplementary 3, n = %d replicates per cell). ", NREP / 12),
        ACC3_DEF),
    D6_detection_demography = paste0(
        sprintf("Detection plane (as in Main a) for the five demographic scenarios of the simulation deposit, named by their deme-size (N) and migration (m) patterns (n = %d replicates each). ", NREP / 5),
        DETECTION_DEF, median_only(DET, "block")),
    S_rda_correction = paste0(
        "Why RDA enters the marker panels without structure correction. ",
        sprintf("Ranking of the causal loci by RDA with and without structure correction, by genetic architecture (n = %d replicates per level), LFMM for reference. ", NREP / 3),
        sprintf("Without correction RDA ranks the causal loci better in %s of %s replicates (median AUC-PR %s vs %s, %s-fold). ",
                RDA_FACTS$wins, RDA_FACTS$n, RDA_FACTS$unc, RDA_FACTS$cor, RDA_FACTS$ratio),
        sprintf("On this landscape population structure follows the climate gradients: the conditioning principal components are collinear with the predictors (largest variance inflation factor ≥ 10 in %s of %s replicates), and the corrected model keeps a median %s of the uncorrected model's adjusted R². ",
                RDA_FACTS$vif10, RDA_FACTS$n, RDA_FACTS$kept),
        "AUC-PR: mean precision at each causal locus in the list of SNPs ranked by p-value (causal loci not retrieved count 0); ",
        "background-neutral loci are the only false positives and linked-neutral loci are excluded. ",
        "With correction: partial RDA conditioned on the first K principal components (K = number of ancestral populations, at least 2); ",
        "without correction: the same model with no conditioning. ",
        "Half-violins: density; boxes: median and interquartile range; dots: replicates. ",
        "Each architecture has its own x-axis range; compare within a panel."))
stopifnot(setequal(names(CAP), names(FIGS)))

ORDER <- c("A1_detection_plane", "R1_rank_distribution", "C2_offset_distribution",
           "R2_rank_by_architecture", "D5_detection_cells", "D5_offset_cells",
           "D6_detection_demography", "S_rda_correction")
stopifnot(setequal(ORDER, names(FIGS)))
TAG <- c(A1_detection_plane = "Main a", R1_rank_distribution = "Main b",
         C2_offset_distribution = "Supplementary 1", R2_rank_by_architecture = "Supplementary 2",
         D5_detection_cells = "Supplementary 3", D5_offset_cells = "Supplementary 4",
         D6_detection_demography = "Supplementary 5", S_rda_correction = "Supplementary 6")
TITLE <- c(A1_detection_plane = "Detection: true vs false positives",
           R1_rank_distribution = "Offset: rank of each panel within a replicate",
           C2_offset_distribution = "Offset accuracy, panels and causal loci",
           R2_rank_by_architecture = "Offset rank by genetic architecture",
           D5_detection_cells = "Detection by design cell",
           D5_offset_cells = "Offset accuracy by design cell",
           D6_detection_demography = "Detection by demography",
           S_rda_correction = "Why uncorrected RDA: GEA ranking with vs without structure correction")
md <- unlist(lapply(ORDER, function(s) c(
    sprintf("## %s — %s", TAG[[s]], TITLE[[s]]), "",
    sprintf("`%s.svg` / `.png` · %.0f × %.0f mm · source: `%s`", s, FIGS[[s]]$w, FIGS[[s]]$h,
            FIGS[[s]]$source), "", CAP[[s]], "")))
writeLines(c("# Simulation figures — draft captions",
             "", sprintf("SS-Clines corpus, %d replicates. Arial 7 pt, drawn at print size. Numbers: `numbers.tsv`; for Supplementary 6, `figures_ssclines_rda/numbers.tsv` and `RDA_CORRECTION_BRIEF.md`.", NREP), "", md),
           file.path(OUT, "captions.md"))

NUMS <- rbindlist(NUM, use.names = TRUE, fill = TRUE)
setcolorder(NUMS, intersect(c("figure", "rule", "stratum", "engine", "metric", "value", "lo", "hi",
                              "expected", "n", "source"), names(NUMS)))
fwrite(NUMS, file.path(OUT, "numbers.tsv"), sep = "\t")
message(sprintf("  numbers.tsv %d rows", nrow(NUMS)))

# Figures dropped at the review must not linger in FIG_OUT (index and zip list ORDER only, but a
# stale file in the folder reads as current).
for (s in c("C5_offset_reach", "S_size_curve")) unlink(file.path(OUT, paste0(s, c(".svg", ".png"))))

# =============================================================================
# Text-only results (ranking reliability) + how to reproduce this report -> README.md + index.html
# Wording is a draft for the Methods / Results; every number comes from the reliability tables.
# =============================================================================
f_ <- function(x, d) formatC(x, format = "f", digits = d)
KW  <- REN[test == "kendall_w"]; PRS <- REN[test == "pair"]; M3R <- RSH[engine == "median of 3"]
stopifnot(nrow(KW) == 1L, nrow(PRS) == 3L, nrow(M3R) == 1L)
TXT_METHODS <- sprintf(paste0(
    "To test whether within-replicate ranks reflect reproducible differences rather than evaluation noise, ",
    "we split the 100 landscape gardens of each replicate at random into two halves (%d splits), recomputed ",
    "each panel's accuracy on each half, centred it on the replicate's mean over its six panels, and correlated ",
    "the two halves (Spearman–Brown-corrected to 100 gardens). We also measured agreement of the panel order ",
    "across the three offset engines with Kendall's W, against a permutation null in which each engine's panel ",
    "values were shuffled within the replicate (%s permutations). Replicates with an empty panel were excluded ",
    "(%d of %d analysed)."),
    M3R$n_splits, format(KW$n_perm, big.mark = ","), M3R$n_rep, NREP)
TXT_RESULTS <- sprintf(paste0(
    "Mean accuracy barely differed among panels (medians within %s; the panel explained %s%% of the variance, ",
    "the replicate %s%%), but within a replicate the panels differed reproducibly: the split-half correlation ",
    "was %s–%s per engine (%s for the median of the three engines at 100 gardens), the least accurate panel was ",
    "the same in both halves in %s–%s%% of replicates (chance %s%%), and the three offset engines agreed on the ",
    "panel order (Kendall's W = %s; %s under the null; P ≤ %s). Which panel was best changed between replicates ",
    "(median range of the six panels within a replicate %s), which is why the pooled distributions overlap ",
    "while the ranks do not."),
    f_(RVP[1, panel_median_range], 3), f_(RVP[source == "panel", share_pct], 2),
    f_(RVP[source == "replicate", share_pct], 1),
    f_(min(RSH[engine != "median of 3", r_half]), 2), f_(max(RSH[engine != "median of 3", r_half]), 2),
    f_(M3R$r_full, 2), f_(min(RSH$worst_same_pct), 0), f_(max(RSH$worst_same_pct), 0), f_(100 / N_PANEL, 0),
    f_(KW$value, 2), f_(KW$null_mean, 2), formatC(KW$p_perm, format = "g", digits = 1),
    f_(RVP[1, within_rep_range_median], 3))
TXT_CAVEAT <- paste0(
    "Not covered by these tests: whether the 2/3 panel ranks mid-field because it is built from the other ",
    "panels' SNPs. Present it as a hedge against the worst case, not as a rule that selects better loci.")

# Reproduce: every script between the offset13 panels and this page. Paths relative to the repo
# root /mnt/data/eugene/ADAPTOGENE; journal = work/journal/18_ssclines_onefit_rule.html.
REPRO <- data.table(
    stage = c(rep("this report (driver: bash work/p5_figures.sh all)", 4), rep("upstream (commands in the journal)", 6)),
    script = c("benchmarks/mvp_j16_gallery.R", "benchmarks/mvp_panel_rank.R",
               "benchmarks/mvp_rank_reliability.R", "benchmarks/mvp_ms_sim_figures.R",
               "benchmarks/mvp_build_snp_sets.R",
               "work/p4b_configs.sh + work/p4b_sweep.sh, then benchmarks/eval_offset_lind.R + benchmarks/mvp_panel_tables.R",
               "benchmarks/mvp_remeasure_600.R (SHARE_MULT=RDA=2)",
               "benchmarks/mvp_oracle_stats.R (+ MVP_SIZE_PANELS=1 on block b1)",
               "benchmarks/mvp_offset_engines.R",
               "benchmarks/mvp_ms_rda_correction.R"),
    output = c("benchmarks/mvp_eval/figures_ssclines_gallery13/ — per-replicate detection and offset tables",
               "benchmarks/mvp_eval/figures_ssclines_offset13/rank/ — rank_long, rank_summary",
               "benchmarks/mvp_eval/figures_ssclines_offset13/rank/reliability_*.tsv",
               "benchmarks/mvp_eval/figures_ssclines_ms/ — figures, captions.md, numbers.tsv, README.md, this page",
               "benchmarks/mvp_eval/offset13/snp_sets_summary.tsv + the six panels per replicate",
               "benchmarks/mvp_eval/offset13/ — garden_performance.tsv, phase1_seed_medians_solo.tsv, panel_pr_recomputed.tsv",
               "benchmarks/mvp_eval/remeasure600/rdaunc_rda2x/ — detection at every share and window",
               "benchmarks/mvp_eval/figures_ssclines_offset13/ (+ size_b1/) — delta vs causal loci, size curve",
               "benchmarks/mvp_eval/figures_ssclines_offset13/ — backward gardens, RDA engine corrected vs uncorrected",
               "benchmarks/mvp_eval/figures_ssclines_rda/ — S_rda_correction, numbers.tsv, RDA_CORRECTION_BRIEF.md"),
    journal = c("Step 19", "Step 18", "Step 21", "Steps 19–21", "Step 13", "Steps 14–15", "Steps 14, 16",
                "Step 17", "Step 17", "Steps 12, 16"))
REPRO_NOTE <- c(
    "Repo root on lab (host fidel): /mnt/data/eugene/ADAPTOGENE. Image cline-go:latest (93bd6025); docker via `nix shell nixpkgs#docker-client -c docker`.",
    "Corpus: MVP_ARM=primary_ssclines, MVP_ADDED = the five SS-Clines block tags, MVP_N_EXPECT=600 (benchmarks/mvp_arm.R). Every arm input is passed explicitly — the scripts' defaults name older arms.",
    "Journal: work/journal/18_ssclines_onefit_rule.html (Steps 12–21). Logs: work/p5_logs/. Seeds: jitter 16, bootstrap 16 / 20261003, splits and permutations 20261003.",
    "Reproduction check 2026-10-03: re-running the rank step reproduced all 8 rank / reliability tables byte for byte.",
    "Superseded set (journal-16 panels): figures_ssclines_ms/_offset12/ (read-only).")
writeLines(c("# Simulation results — report README",
             "", sprintf("SS-Clines, %d replicates. Contains every figure, caption and number needed to write the simulation part of the manuscript's Results (and the matching Methods). Open `index.html` for the figures; numbers in `numbers.tsv`; captions in `captions.md`.", NREP),
             "", "## Text-only result: the within-replicate ranking is signal, not noise", "",
             "**Methods (draft).** ", TXT_METHODS, "", "**Results (draft).** ", TXT_RESULTS, "", paste0("*", TXT_CAVEAT, "*"),
             "", "## Reproduce", "", paste0("- ", REPRO_NOTE), "",
             "| stage | script | output | journal 18 |", "|---|---|---|---|",
             REPRO[, sprintf("| %s | `%s` | %s | %s |", stage, script, output, journal)]),
           file.path(OUT, "README.md"))

# =============================================================================
# Download-all zip + index.html
# =============================================================================
files <- c(paste0(ORDER, ".svg"), paste0(ORDER, ".png"), "captions.md", "numbers.tsv", "README.md")
ZIP <- "sim_figures_ms.zip"
old <- setwd(OUT)
unlink(ZIP)
stopifnot(utils::zip(ZIP, files, flags = "-9Xq") == 0L)
setwd(old)

esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
                     gsub(">", "&gt;", x, fixed = TRUE) }
kb <- function(f) { b <- file.size(file.path(OUT, f)); if (b > 1e6) sprintf("%.1f MB", b / 1e6) else sprintf("%.0f kB", b / 1e3) }
card <- function(s, cls = "") sprintf(paste0(
    '<figure class="card %s" id="%s">',
    '<figcaption class="head"><span class="tag">%s</span><span class="ttl">%s</span>',
    '<span class="dim">%.0f × %.0f mm</span></figcaption>',
    '<div class="img"><img src="%s.svg" alt="%s" style="--w:%.0f" loading="lazy"></div>',
    '<p class="cap">%s</p>',
    '<div class="dl"><a class="btn" href="%s.svg" download>SVG · %s</a>',
    '<a class="btn" href="%s.png" download>PNG 600 dpi · %s</a></div></figure>'),
    cls, s, esc(TAG[[s]]), esc(TITLE[[s]]), FIGS[[s]]$w, FIGS[[s]]$h, s, esc(TITLE[[s]]),
    FIGS[[s]]$w, esc(CAP[[s]]), s, kb(paste0(s, ".svg")), s, kb(paste0(s, ".png")))

html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
'<meta name="viewport" content="width=device-width, initial-scale=1">',
'<title>Simulation figures</title><style>',
':root{--bg:#f5f6f8;--card:#fff;--fg:#1f2933;--muted:#5f6b7a;--line:#dde2e8;--accent:#00798c;--accent-fg:#fff}',
'@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#15181d;--card:#1f242b;--fg:#e6e9ee;--muted:#9aa5b3;--line:#323a44;--accent:#3fa58a;--accent-fg:#0d1512}}',
'*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",Arial,sans-serif}',
'header{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--line);padding:14px 24px;display:flex;flex-wrap:wrap;gap:12px;align-items:center;justify-content:space-between}',
'header h1{font-size:18px;margin:0}header p{margin:2px 0 0;color:var(--muted);font-size:13px}',
'main{max-width:1500px;margin:0 auto;padding:20px 24px 60px}',
'h2{font-size:15px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:28px 0 12px}',
'.note{color:var(--muted);font-size:13px;margin:-6px 0 14px}',
'.stack{display:grid;gap:16px}',
'.card{margin:0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px;display:flex;flex-direction:column;min-width:0}',
'.head{display:flex;flex-wrap:wrap;gap:8px;align-items:baseline;margin-bottom:10px}',
'.tag{background:var(--accent);color:var(--accent-fg);font-size:12px;font-weight:600;padding:2px 8px;border-radius:999px}',
'.ttl{font-weight:600}.dim{color:var(--muted);font-size:12px;margin-left:auto}',
'.img{background:#fff;border-radius:6px;padding:6px;display:flex;justify-content:center;overflow-x:auto}',
'.img img{width:calc(100% * var(--w) / 180);height:auto;display:block}',
'body.print .img img{width:calc(var(--w) * 1mm);max-width:none}',
'.cap{color:var(--muted);font-size:13px;margin:10px 0 12px;flex:1}',
'.dl{display:flex;flex-wrap:wrap;gap:8px}',
'section.card p{margin:0 0 10px}code{font-size:12px;overflow-wrap:anywhere}table.repro{border-collapse:collapse;font-size:13px;width:100%}table.repro th,table.repro td{border-top:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}table.repro th{color:var(--muted);font-weight:600}',
'.btn{display:inline-block;cursor:pointer;font:inherit;text-decoration:none;font-size:13px;padding:6px 12px;border-radius:7px;border:1px solid var(--line);color:var(--fg);background:var(--bg)}',
'.btn:hover{border-color:var(--accent)}.btn.primary{background:var(--accent);color:var(--accent-fg);border-color:var(--accent);font-weight:600}',
'@media (max-width:900px){main,header{padding-left:16px;padding-right:16px}}',
'</style></head><body>',
sprintf('<header><div><h1>Simulation figures — SS-Clines, %d replicates</h1><p>Editable SVG (Arial 7 pt, print size) and PNG 600 dpi; captions are drafts. Panels: top 0.25%% per method and predictor (RDA 0.5%%), RDA without structure correction, 5 kb. Built %s.</p></div>', NREP, format(Sys.time(), "%Y-%m-%d %H:%M")),
sprintf('<div class="dl"><button class="btn" id="ps" type="button" aria-pressed="false">Show print size</button><a class="btn" href="captions.md" download>captions.md</a><a class="btn" href="numbers.tsv" download>numbers.tsv</a><a class="btn" href="README.md" download>README.md</a><a class="btn primary" href="%s" download>Download all (.zip, %s)</a></div></header>', ZIP, kb(ZIP)),
'<main><h2>Main figure — proposed layout</h2>',
'<p class="note">Row 1: detection. Row 2: offset, rank of each panel within a replicate. Both 180 mm wide.</p>',
'<section class="stack">', card("A1_detection_plane"), card("R1_rank_distribution"), '</section>',
'<h2>Supplementary</h2>',
'<section class="stack">', card("C2_offset_distribution"), card("R2_rank_by_architecture"),
card("D5_detection_cells"), card("D5_offset_cells"), card("D6_detection_demography"),
card("S_rda_correction"), '</section>',
'<h2>Text-only result — the within-replicate ranking is signal, not noise</h2>',
'<section class="card"><p class="note" style="margin:0 0 8px">No figure (user decision 2026-10-03): these tests go into the manuscript text. Draft wording; numbers from <code>figures_ssclines_offset13/rank/reliability_*.tsv</code>.</p>',
sprintf('<p><strong>Methods.</strong> %s</p>', esc(TXT_METHODS)),
sprintf('<p><strong>Results.</strong> %s</p>', esc(TXT_RESULTS)),
sprintf('<p class="note"><em>%s</em></p>', esc(TXT_CAVEAT)), '</section>',
'<h2>Reproduce</h2><section class="card">',
paste0('<ul class="note" style="margin:0 0 10px;padding-left:18px">', paste(sprintf('<li>%s</li>', esc(REPRO_NOTE)), collapse = ""), '</ul>'),
'<div style="overflow-x:auto"><table class="repro"><thead><tr><th>stage</th><th>script</th><th>output</th><th>journal 18</th></tr></thead><tbody>',
paste(REPRO[, sprintf('<tr><td>%s</td><td><code>%s</code></td><td>%s</td><td>%s</td></tr>', esc(stage), esc(script), esc(output), esc(journal))], collapse = ""),
'</tbody></table></div></section></main>',
# Every image is shown at one common scale (width in mm / 180 of the row), so relative panel
# sizes stay true; "print size" shows each at its physical width (CSS mm) to judge 7 pt text.
'<script>document.getElementById("ps").addEventListener("click",function(){var on=document.body.classList.toggle("print");this.setAttribute("aria-pressed",on);this.textContent=on?"Fit to page":"Show print size";});</script>',
'</body></html>')
writeLines(html, file.path(OUT, "index.html"))

message(sprintf("OK: %d figures (SVG + PNG), captions.md, numbers.tsv, index.html, %s -> %s",
                length(FIGS), ZIP, OUT))
