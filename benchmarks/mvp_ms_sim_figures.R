#!/usr/bin/env Rscript
# =============================================================================
# mvp_ms_sim_figures.R -- manuscript display items for the simulation section (SS-Clines, 600
# replicates): figures, captions.md, Table S1, Data S2, numbers.tsv, README.md, index.html.
# Rebuilt 2026-10-03 (Phase 5 of the SS-Clines re-analysis) on the offset13 panels: top 0.25 % of
# SNPs per method and per predictor (RDA, one multivariate p for both predictors, 0.5 %), RDA
# without structure correction as the third method, 5 kb agreement window. The 2026-09-25 set
# (journal-16 panels, offset12) is archived read-only in FIG_OUT/_offset12/ and the 2026-10-03 set
# in FIG_OUT/_2026-10-03/; git history holds the scripts that drew them.
#
# [changed 2026-10-14, user] Manuscript numbering and vocabulary, one tie rule, four new items:
#   * numbering: Figure 2 = a (A1) over b (R1) in one file, S1-S8 as listed below, Table S1,
#     Data S2. D5_offset_cells is dropped (its files moved by hand to FIG_OUT/_removed/).
#   * vocabulary: the six panels are CANDIDATE SETS -- three single methods and three AGREEMENT
#     SETS ("1/3 set", "2/3 set", "3/3 set"); true positives are "causal or linked SNPs", false
#     positives "neutral SNPs" (unlinked to any causal locus); 0.25 % is "the cutoff". The data
#     keys ("1/3 methods", ...) stay as the input tables write them; SET_LAB maps them for display,
#     and vocab_check() stops the script on banned wording in any caption, page or figure text.
#   * ties: within-replicate ranks split ties evenly (benchmarks/mvp_rank_ties.R) -- Figure 2b,
#     S7, S8 and Table S1 here, rank_summary.tsv upstream (mvp_panel_rank.R).
#   * new: S2 AUC-PR per method, S3 cutoff sensitivity, S8 least-accurate set by design cell and
#     demographic scenario (+ Kruskal-Wallis), Table S1 (replaces an ad hoc python script that was
#     never kept in the repo), Data S2 (every source TSV + a manifest).
#
# PLOTTING + TABLES. Every number comes from tables already on disk (compute once, reuse); the one
# computation is Table S1's bootstrap CIs and permutation P (seed 20261005, 10,000 each):
#   GALLERY_DIR  mvp_j16_gallery.R run on offset13 (work/p5_figures.sh): detection_per_seed,
#                offset_delta_per_seed, detection_stats, C_delta_stats, D5_offset_stats,
#                B4_aucpr_per_seed + B4_aucpr_stats (S2)
#   RANK_DIR     mvp_panel_rank.R: rank_long, rank_summary (within-replicate rank, per engine);
#                mvp_rank_reliability.R: reliability_{split_half,engines,variance} (numbers only)
#   SIZE_DIR     mvp_oracle_stats.R with MVP_SIZE_PANELS=1 on block b1: delta_per_seed (numbers only)
#   OFFSET_DIR   offset13 (name under benchmarks/mvp_eval): panel_pr_recomputed, snp_sets_summary
#   RDA_FIG_DIR  mvp_ms_rda_correction.R: S_rda_correction.{svg,png} (copied as is), numbers.tsv
#                (caption facts, S2's AUC-PR tie) and detection_table.tsv -- the remeasure
#                harness's medians, an independent tie for the panel builder's counts plotted here
#   REMEASURE_DIR  mvp_remeasure_600.R arm rdaunc_rda2x (RDA x 2): calls_per_method,
#                calls_combine at the five cutoffs (S3) + provenance.tsv
# Medians and counts re-derived here are regression ties or caption facts only. NO INPUT HAS A
# DEFAULT: older arms (offset12_ssclines_pooled, figures_ssclines_j16_gallery) still exist on disk
# and a forgotten variable would give a complete, successful, wrong figure set.
#
#   Figure 2   Fig2_simulation           a = A1_detection_plane (TP vs FP plane: 50 % bags + medians,
#                                        own x-axis per plot) over b = R1_rank_distribution (rank of
#                                        the six sets inside each replicate, per offset engine);
#                                        A1 and R1 are also written as separate files
#   Figure S1  S_rda_correction          why RDA enters WITHOUT structure correction (copied)
#   Figure S2  B4_aucpr_methods          AUC-PR of LFMM, RDA, EMMAX by genetic architecture
#   Figure S3  A5_cutoff_sensitivity     median TP and FP of the six sets at cutoffs 0.1-2 %
#   Figure S4  D5_detection_cells        Figure 2a view, genic level x (pleiotropy x selection regime)
#   Figure S5  D6_detection_demography   Figure 2a view, the five demography blocks
#   Figure S6  C2_offset_distribution    absolute offset accuracy per set + a separate causal-loci track
#   Figure S7  R2_rank_by_architecture   Figure 2b split by genic level
#   Figure S8  R3_least_accurate_strata  % replicates least accurate, by design cell (a) and demography (b)
#   Table S1   Table_S1.{tsv,md}         rank (a), delta vs causal loci (b), stratum ranges (c),
#                                        variance partition (d), Kruskal-Wallis tests (e)
#   Data S2    Data_S2.zip + Data_S2_manifest.tsv
# Dropped at the 2026-10-03 review (user): C5_offset_reach (reach tile) and S_size_curve (its
# numbers stay in numbers.tsv); at the 2026-10-14 renumbering (user): D5_offset_cells. RDA WITH
# structure correction as an OFFSET ENGINE is in no figure, caption or table here (Phase 4b gate);
# as a GEA method it appears only in S_rda_correction.
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
# Also writes captions.md (Figure 2, S1-S8, then Table S1), numbers.tsv (the numbers a caption or
# the Results text needs, long format, with source), Table_S1.{tsv,md}, Data_S2.zip +
# Data_S2_manifest.tsv, README.md (text-only ranking-reliability result + how to reproduce every
# input of this report), index.html (served via work/journal/sim_figures on port 8099) and
# sim_figures_ms.zip (every SVG + PNG + captions.md + numbers.tsv + README.md + Table S1).
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
REM  <- need("REMEASURE_DIR")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_ms"))
JITTER_SEED <- 16L
S1_SEED <- 20261005L                     # Table S1: bootstrap and permutation, set before each
S1_B    <- 10000L                        # Table S1: resamples and permutations
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
source(file.path(ROOT, "benchmarks/mvp_rank_ties.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
# Figures dropped at a review must not linger in FIG_OUT (a stale file reads as current). They are
# MOVED to FIG_OUT/_removed/ by hand, never deleted here (user, 2026-10-14): stop if one is left.
DROPPED <- c("C5_offset_reach", "S_size_curve", "D5_offset_cells")
left_over <- DROPPED[file.exists(file.path(OUT, paste0(DROPPED, ".svg"))) |
                     file.exists(file.path(OUT, paste0(DROPPED, ".png")))]
if (length(left_over)) stop("dropped figure(s) still in FIG_OUT, move them to _removed/: ",
                            paste(left_over, collapse = ", "))

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
# [2026-10-18] display labels follow the manuscript text; the keys stay the data keys.
ENG_LAB <- c(gradientForest = "Gradient Forest", LFMM2 = "Geometric offset", RDA = "RDA-based offset",
             `median of 3` = "Median of the three")
N_PANEL <- length(RULES)
# Manuscript vocabulary (user, 2026-10-14). The six things compared are CANDIDATE SETS: three single
# methods and three AGREEMENT SETS. RULES stay the data keys of every input table; whatever a reader
# sees -- scales, legends, tags, captions, Table S1 -- goes through SET_LAB.
SET_LAB <- c("2/3 methods" = "2/3 set", "1/3 methods" = "1/3 set", "3/3 methods" = "3/3 set",
             "LFMM" = "LFMM", "RDA" = "RDA", "EMMAX" = "EMMAX")
stopifnot(setequal(names(SET_LAB), RULES))
set_lab <- function(x) unname(SET_LAB[as.character(x)])
TP_LAB <- "True positives (causal or linked SNPs)"
FP_LAB <- "False positives (neutral SNPs)"
DELTA_LAB <- "Δ accuracy (candidate set − causal loci)"
# Wording a reader must never see (user, 2026-10-14): "panel" for a candidate set, "rule" for an
# agreement set, "methods" in a set name, "background", "linked-neutral" / "background-neutral",
# linkage-group numbers, "share" for the cutoff or for composition. Figure panels are "plots".
BANNED <- c(set_name = "/3 methods", panel = "\\bpanels?\\b", rule = "(?<!tie )\\brules?\\b",
            background = "background", linked_neutral = "linked-neutral",
            linkage_group = "linkage group", share = "\\bshares?\\b")
# Reader prose minus code spans and file paths (a script name may say "panel").
strip_code <- function(x) {
    x <- gsub("`[^`]*`", " ", x)
    x <- gsub("<code>.*?</code>", " ", x, perl = TRUE)
    gsub("[^ ()]*[_/][^ ()]*\\.(R|tsv|md|html|sh|zip|svg|png)\\b", " ", x, perl = TRUE)
}
vocab_check <- function(x, where) {
    x <- strip_code(as.character(unlist(x)))
    hit <- Filter(length, lapply(BANNED, function(b) grep(b, x, ignore.case = TRUE, perl = TRUE, value = TRUE)))
    if (length(hit)) stop(where, ": banned wording for a reader: ",
                          paste(sprintf("%s in \"%s\"", names(hit), vapply(hit, function(h) substr(h[1], 1, 160), "")),
                                collapse = "; "))
    invisible(length(x))
}
# Every string drawn in a figure: walk the gtable (titles, ticks, strips, legends, text and label
# layers, ggrepel's label trees). Used on each plot before it is saved.
grob_text <- function(g) {
    out <- character()
    walk <- function(x) {
        if (!is.null(x$label) && (is.character(x$label) || is.expression(x$label) || is.factor(x$label)))
            out <<- c(out, as.character(x$label))
        if (is.data.frame(x$data) && "label" %in% names(x$data)) out <<- c(out, as.character(x$data$label))
        if (inherits(x, "gtable")) for (y in x$grobs) walk(y)
        if (length(x$children)) for (y in x$children) walk(y)
    }
    walk(g)
    unique(out[!is.na(out) & nzchar(out)])
}

TEXT_TONES <- c("#FFFFFF" = "#FFFFFF", setNames(MS_INK, MS_INK))
ZERO <- geom_vline(xintercept = 0, linetype = "dashed", colour = MS_REF, linewidth = 0.35)

FIGS <- list()                                  # stem -> list(role, w, h, source)
# Figure text through the gtable, on a null device (no Rplots.pdf from a text measurement).
fig_text <- function(p) { grDevices::pdf(NULL); on.exit(grDevices::dev.off()); grob_text(ggplotGrob(p)) }
# role: "main" / "supp" = a display item with its own caption; "part" = a panel of Figure 2 also
# written as its own file (in the zip and on the Figure 2 card, captioned under Figure 2).
save_ms <- function(stem, p, w, h, role, source) {
    n_txt <- vocab_check(fig_text(p), paste0(stem, " (figure text)"))
    ms_save(file.path(OUT, stem), p, w, h)
    FIGS[[stem]] <<- list(role = role, w = w, h = h, source = source)
    message(sprintf("  %-26s %3.0f x %3.0f mm  (%d text strings checked)", stem, w, h, n_txt))
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

# Numbers the manuscript already quotes, fixed before the 2026-10-14 display work; that work must
# not move them (user). AUC-PR is asserted in the S2 section, 557 / 600 at the S1 copy.
FIX_DET <- data.table(label = c("1/3 methods", "2/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX"),
                      TP = c(213, 178, 111, 109, 108, 71), FP = c(89, 36, 13.5, 16.5, 16, 62))
fx <- merge(DET[, .(TPm = median(TP), FPm = median(FP)), by = .(label = as.character(label))],
            FIX_DET, by = "label")
stopifnot(nrow(fx) == N_PANEL, all(fx$TPm == fx$TP), all(fx$FPm == fx$FP),
          nrow(EMPTY) == 1L, EMPTY$label == "3/3 methods", EMPTY$N == 2L)
message("fixed numbers: Figure 2a medians (TP / FP) of all six sets and 2 empty 3/3 sets -- unchanged")

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
# [2026-10-14] The tie rule must not reorder the sets. The order under the previous count (last =
# rank >= 6, last two = rank >= 5, re-derived here from rank_long) must equal the order above
# (rank_summary.tsv, ties split evenly); the script stops if it does not.
ORD_OLD <- RL[engine == "median of 3", .(bottom1 = 100 * mean(rank >= N_PANEL),
                                         bottom2 = 100 * mean(rank >= N_PANEL - 1)), by = panel]
setorder(ORD_OLD, bottom1, bottom2)
message("set order, previous count (rank >= 6, >= 5): ",
        paste(sprintf("%s %.3f / %.3f%%", ORD_OLD$panel, ORD_OLD$bottom1, ORD_OLD$bottom2), collapse = " < "))
message("set order, ties split evenly:                ",
        paste(sprintf("%s %.3f / %.3f%%", ORD$panel, ORD$bottom1, ORD$bottom2), collapse = " < "))
if (!identical(ORD_OLD$panel, PANEL_ORDER)) stop("the tie rule changed the order of the candidate sets")
# [2026-10-18] the order of Figure 2b, as in the 2026-10-14 build (display work must not move it)
stopifnot(identical(PANEL_ORDER, c("2/3 methods", "1/3 methods", "LFMM", "3/3 methods", "RDA", "EMMAX")))
message("set order unchanged since 2026-10-14: ", paste(PANEL_ORDER, collapse = " < "))

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

# y0 = FALSE (A1 only, 2026-10-18): the true-positive axis starts just below the lowest bag instead of at 0;
# the false-positive axis keeps its 0.
plane_ms <- function(D, facet, y_lab = TP_LAB, y0 = TRUE) {
    keys <- c("label", facet$vars)
    M  <- D[, .(TP = median(TP), FP = median(FP)), by = keys]
    B  <- bags_of(D, keys)
    p <- ggplot(D, aes(FP, TP)) +
        geom_polygon(data = B, aes(group = grp, colour = label), fill = NA, linewidth = 0.35) +
        geom_point(data = M, aes(fill = label), shape = 21, size = 1.7, colour = "black", stroke = 0.3) +
        # Median labels are the candidate-set names (SET_LAB) on solid tags in the set colour; the
        # text is white or dark, whichever contrasts more with that colour (ms_text_on). The text
        # tone rides on the same colour scale as the bags, under its own two keys.
        geom_label_repel(data = M[, tone := ms_text_on(RULE_COL[as.character(label)])],
                         aes(label = set_lab(label), fill = label, colour = tone), size = MS_LAB,
                         fontface = "bold", label.size = 0, label.padding = 0.12, box.padding = 0.3,
                         min.segment.length = 0, segment.size = 0.25, segment.colour = MS_INK,
                         seed = JITTER_SEED, show.legend = FALSE, max.overlaps = Inf)
    # Both count axes start at 0 (the scales train on bags + medians only); with y0 = FALSE only x does.
    p + facet$layer + (if (y0) expand_limits(x = 0, y = 0) else expand_limits(x = 0)) +
        scale_x_sqrt(breaks = sqrt_even_breaks) + scale_y_sqrt(breaks = sqrt_even_breaks) +
        scale_colour_manual(values = c(RULE_COL, TEXT_TONES), guide = "none") +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        labs(x = FP_LAB, y = y_lab) +
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
           paste(sprintf("%s (%d of %d plots)", set_lab(x$label), x$N,
                         nrow(unique(D[, vars, with = FALSE]))), collapse = ", "), ".")
}
# Every panel of a detection plane: tick labels present and evenly spaced on the sqrt scale
# (`axes`: the square-root axes to check; S3's x axis is logarithmic).
check_ticks <- function(p, what, axes = c("x", "y")) {
    for (pp in ggplot_build(p)$layout$panel_params) for (ax in axes) {
        b <- suppressWarnings(as.numeric(pp[[ax]]$get_labels()))
        b <- b[!is.na(b)]
        g <- diff(sqrt(b))
        if (length(b) < 2L || max(abs(g - g[1])) > 1e-9)
            stop(what, ": ", ax, " ticks not evenly spaced on the sqrt scale: ", paste(b, collapse = " "))
    }
    message(sprintf("%-24s every panel: %s ticks evenly spaced on the sqrt scale", what,
                    paste(axes, collapse = " and ")))
}

message("=== detection")
# The long true-positive title has to fit the 62 mm plot: its width at 7 pt (null pdf device:
# Helvetica, metric-compatible with the Arial / Nimbus Sans of the SVG) against the height left
# for the panels in the built gtable. Wrapped onto two lines only if it does not fit.
y_title_fits <- function(p, h_mm, label) {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    g <- ggplotGrob(p)
    fixed <- sum(grid::convertHeight(g$heights[grid::unitType(g$heights) != "null"], "mm", valueOnly = TRUE))
    need  <- grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = MS_BASE))),
                                "mm", valueOnly = TRUE)
    message(sprintf("A1 y title: %.1f mm of text, %.1f mm of panel height", need, h_mm - fixed))
    need <= h_mm - fixed
}
# Every detection plane gives each panel its own false-positive axis (scales = "free_x").
pA1 <- plane_ms(DET, list(vars = "arch", layer = facet_grid(cols = vars(arch), scales = "free_x")), y0 = FALSE)
if (!y_title_fits(pA1, 62, TP_LAB))
    pA1 <- plane_ms(DET, list(vars = "arch", layer = facet_grid(cols = vars(arch), scales = "free_x")),
                    y_lab = sub(" (", "\n(", TP_LAB, fixed = TRUE), y0 = FALSE)
for (pp in ggplot_build(pA1)$layout$panel_params)
    message("A1 panel axis labels  x: ", paste(pp$x$get_labels(), collapse = " "),
            "   y: ", paste(pp$y$get_labels(), collapse = " "))
# [2026-10-18] y no longer pinned at 0: it must start above 0 and below the lowest bag vertex / median
A1_LOW <- min(bags_of(DET, c("label", "arch"))$TP, DET[, .(TP = median(TP)), by = .(label, arch)]$TP)
for (pp in ggplot_build(pA1)$layout$panel_params) {
    y_lo <- pp$y$continuous_range[1]^2                    # range is on the sqrt scale
    message(sprintf("A1 panel y range  %.1f to %.1f (lowest bag vertex or median %.1f)",
                    y_lo, pp$y$continuous_range[2]^2, A1_LOW))
    stopifnot(pp$y$continuous_range[1] > 0, y_lo < A1_LOW)
}
check_ticks(pA1, "A1_detection_plane")
save_ms("A1_detection_plane", pA1, 180, 62, "part", "detection_per_seed.tsv")

pD5d <- plane_ms(DET, list(vars = "panel", layer = facet_wrap(vars(panel), ncol = 4, scales = "free_x")))
check_ticks(pD5d, "D5_detection_cells")
save_ms("D5_detection_cells", pD5d, 180, 165, "supp", "detection_per_seed.tsv")

pD6 <- plane_ms(DET, list(vars = "block", layer = facet_wrap(vars(block), nrow = 2, scales = "free_x")))
check_ticks(pD6, "D6_detection_demography")
save_ms("D6_detection_demography", pD6, 180, 110, "supp", "detection_per_seed.tsv")

num(melt(DS[, .(rule = label, stratum = arch, n = n_seeds, n_empty = as.numeric(n_empty), TP_median, TP_q25, TP_q75,
                FP_median, FP_q25, FP_q75, prec_cl_median, F1_median)],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "Figure 2a / S4 / S5 detection", "figures_ssclines_gallery13/detection_stats.tsv")

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
        scale_y_continuous(breaks = unname(Y_AT), labels = c(set_lab(PANEL_ORDER), ORACLE_LAB)) +
        scale_fill_manual(values = ACC_COL, guide = "none") +
        scale_colour_manual(values = ACC_COL, guide = "none") +
        labs(x = "Offset accuracy (−Kendall τ)", y = NULL) + theme_ms()
}
save_ms("C2_offset_distribution",
        acc_cloud(ACC, facet_grid(cols = vars(arch), labeller = labeller(arch = ARCH_WRAP))),
        114, 66, "supp", "offset_delta_per_seed.tsv (accuracy, oracle)")
# D5_offset_cells (the same view per design cell) was dropped at the 2026-10-14 renumbering (user).

acc_stats <- function(D) D[, .(n = .N, median = median(accuracy), q25 = quantile(accuracy, .25),
                               q75 = quantile(accuracy, .75))]
ACS <- rbind(ACC[, acc_stats(.SD), by = .(rule = label, stratum = as.character(arch))],
             ACC[, acc_stats(.SD), by = .(rule = label)][, stratum := "pooled"], use.names = TRUE)
num(melt(ACS, id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "Figure S6 offset accuracy (median of 3 engines)", "figures_ssclines_gallery13/offset_delta_per_seed.tsv")
# Delta from the causal loci: no figure any more; the Results text and Table S1 (b) quote it.
num(melt(CS[, .(rule = label, stratum = arch, n, median, ci_lo, ci_hi, q25, q75, win_pct,
                reach_pct, p_wilcox_two_sided)],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "offset delta vs causal loci (Table S1 b + text; reach_pct = delta >= 0, win_pct = delta > 0)",
    "figures_ssclines_gallery13/C_delta_stats.tsv")

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
          engine_lab = factor(ENG_LAB[engine], levels = ENG_LAB))]
stopifnot(!anyNA(RL$arch), !anyNA(RL$engine_lab))
# [changed 2026-10-14] Ranks as (rank, weight) rows under the one tie rule (mvp_rank_ties.R): a set
# tied for ranks k and k + 1 counts half at each. Was: drawn at round(rank), half to even (4.5 -> 4,
# 5.5 -> 6). The tie groups found from the ranks must be the ties of the raw accuracies (key), and
# an untied set is one row at its own rank with weight 1 -- exactly what the bars drew before.
RL[, g_key := .GRP, by = .(seed, engine, key)][, g_rank := .GRP, by = .(seed, engine, rank)]
stopifnot(uniqueN(RL[, .(g_key, g_rank)]) == uniqueN(RL$g_key), uniqueN(RL$g_key) == uniqueN(RL$g_rank))
RL[, c("g_key", "g_rank") := NULL]
RX <- rank_split(RL, c("seed", "engine"))
stopifnot(RX[tie_n == 1L, all(r == rank & w == 1)])
TG <- unique(RX[tie_n > 1L, .(seed, engine, rank, tie_n)])          # one row per tie group
N_TIE2 <- TG[tie_n == 2L, .N]; N_TIE3 <- TG[tie_n > 2L, .N]; N_HALF <- RL[, sum(rank != round(rank))]
stopifnot(uniqueN(TG[, .(seed, engine)]) == N_TIED)
message(sprintf("ties: %d of %d replicate x engine units; tie groups of two sets %d, of three or more %d; %d x.5 ranks",
                N_TIED, N_UNITS, N_TIE2, N_TIE3, N_HALF))
RX <- merge(RX, RL[, .(seed, engine, panel, arch, engine_lab)], by = c("seed", "engine", "panel"))

RANK_COL <- setNames(colorRampPalette(c(MINOU[["red"]], "#F2F2F2", MINOU[["teal"]]))(N_PANEL),
                     N_PANEL:1)
rank_dist <- function(D, extra = NULL, facet_layer) {
    R <- D[, .(n = sum(w)), by = c("engine_lab", "panel", "r", extra)]
    R[, pct := 100 * n / sum(n), by = c("engine_lab", "panel", extra)]
    R[, r := factor(r, levels = N_PANEL:1)]
    R[, panel := factor(panel, levels = rev(PANEL_ORDER))]
    ggplot(R, aes(pct, panel, fill = r)) +
        geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
        facet_layer +
        scale_y_discrete(labels = SET_LAB) +
        scale_fill_manual(values = RANK_COL, name = "rank (1 = most accurate)",
                          breaks = as.character(1:N_PANEL)) +
        # [2026-10-18] ticks 0 / 50 / 100 without "%" (the axis title carries the unit) and room between the
        # facets, so neighbouring tick labels never touch ("100%0%")
        scale_x_continuous(breaks = c(0, 50, 100), labels = c("0", "50", "100"), expand = c(0, 0)) +
        labs(x = "% replicates at each rank", y = NULL) + theme_ms() +
        theme(panel.spacing.x = unit(4, "mm"))
}
pR1 <- rank_dist(RX, NULL, facet_wrap(~ engine_lab, nrow = 1))
save_ms("R1_rank_distribution", pR1, 180, 66, "part", "rank_long.tsv")
save_ms("R2_rank_by_architecture",
        rank_dist(RX, "arch", facet_grid(rows = vars(arch), cols = vars(engine_lab))),
        180, 120, "supp", "rank_long.tsv")
# Figure 2 (user, 2026-10-14): a = detection over b = rank, one file; a and b stay as their own
# files too. The panel letters take the family of theme_ms() (the device default the text of every
# plot resolves to), so ms_save()'s font check sees one family.
FIG2 <- cowplot::plot_grid(pA1, pR1, ncol = 1, rel_heights = c(62, 66), labels = c("a", "b"),
                           label_size = 8, label_fontface = "bold",
                           label_fontfamily = theme_ms()$text$family, align = "v", axis = "lr")
save_ms("Fig2_simulation", FIG2, 180, 128, "main", "detection_per_seed.tsv (a), rank_long.tsv (b)")

num(RS[stratum %in% c("pooled", ARCH_LABELS) &
       metric %in% c("top1", "top2", "bottom1", "bottom2", "mean_rank"),
       .(rule = panel, stratum, n, engine = ENG_LAB[engine], metric, value, lo, hi, expected)],
    "Figure 2b / S7 within-replicate rank (CIs from 2,000 resamples; the text quotes Table S1's 10,000)",
    "figures_ssclines_offset13/rank/rank_summary.tsv")
num(data.table(metric = c("units_with_tie", "units", "tie_groups_two_sets", "tie_groups_three_or_more",
                          "half_ranks"),
               value = c(N_TIED, N_UNITS, N_TIE2, N_TIE3, N_HALF)),
    "tie rule: tied sets split the tied ranks evenly (Figure 2b, S7, S8, Table S1)",
    "figures_ssclines_offset13/rank/rank_long.tsv")
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
    "candidate-set order (every per-set figure)", "figures_ssclines_offset13/rank/rank_summary.tsv (bottom1, bottom2)")

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
# Fixed number the manuscript quotes (user, 2026-10-14): uncorrected RDA better in 557 of 600.
stopifnot(RDA_FACTS$wins == "557", RDA_FACTS$n == as.character(NREP))

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
# Figure S2: AUC-PR of the three single methods (gallery code B4), drawn as S1 (user, 2026-10-14)
# =============================================================================
message("=== S2 AUC-PR per method")
AUC <- rd(GAL, "B4_aucpr_per_seed", colClasses = c(seed = "character"))
AUS <- rd(GAL, "B4_aucpr_stats")
SINGLE <- intersect(PANEL_ORDER, c("LFMM", "RDA", "EMMAX"))      # the data-driven order of every figure
AUC[, `:=`(arch = factor(arch, levels = ARCH_LABELS), label = factor(label, levels = SINGLE))]
stopifnot(nrow(AUC) == 3L * NREP, !anyNA(AUC$arch), !anyNA(AUC$label), !anyNA(AUC$aucpr),
          all(AUC$seed %in% PRIM$seed), AUC[, .N, by = .(label, arch)][, all(N == NREP / 3)])
AUM <- rbind(AUC[, .(m = median(aucpr), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
             AUC[, .(m = median(aucpr), n = .N), by = .(label = as.character(label))][, arch := "pooled"])
x <- merge(AUM, AUS[, .(label, arch, median, n)], by = c("label", "arch"), all = TRUE)
stopifnot(nrow(x) == 12L, !anyNA(x), max(abs(x$m - x$median)) < 1e-12, all(x$n.x == x$n.y))
# Tie to the RDA brief (figures_ssclines_rda/numbers.tsv: the same rank_metrics read by another
# script) and to the three values the manuscript quotes (fixed, user 2026-10-14).
AU_KEY <- c(LFMM = "aucpr_LFMM_pooled_median", RDA = "aucpr_RDA_uncorrected_pooled_median",
            EMMAX = "aucpr_EMMAX_pooled_median")
AU_FIX <- c(LFMM = "0.249", RDA = "0.217", EMMAX = "0.047")
for (meth in names(AU_KEY)) {
    v <- AUM[label == meth & arch == "pooled", m]
    vb <- as.numeric(RNUM[key == AU_KEY[[meth]], value])
    if (length(vb) != 1L || abs(v - vb) > 1e-12 || sprintf("%.3f", v) != AU_FIX[[meth]])
        stop("S2: pooled AUC-PR of ", meth, " = ", v, " differs from figures_ssclines_rda/numbers.tsv or ", AU_FIX[[meth]])
}
message(sprintf("tie S2 AUC-PR: 12 medians = B4_aucpr_stats; pooled = figures_ssclines_rda/numbers.tsv (%s)",
                paste(sprintf("%s %.3f", names(AU_KEY), AUM[arch == "pooled"][match(names(AU_KEY), label), m]),
                      collapse = ", ")))
AUC[, ypos := as.numeric(factor(label, levels = rev(SINGLE)))]
pS2 <- ggplot(AUC, aes(aucpr, ypos)) +
    stat_halfeye(aes(fill = label, group = label), orientation = "horizontal", adjust = 0.8,
                 height = 0.55, justification = -0.2, .width = 0, point_colour = NA, slab_alpha = 0.75) +
    geom_boxplot(aes(group = label), orientation = "y", width = 0.14, outlier.shape = NA,
                 fill = "white", colour = MS_INK, linewidth = 0.25) +
    geom_point(aes(y = ypos - 0.25, colour = label),
               position = position_jitter(height = 0.08, width = 0, seed = JITTER_SEED),
               size = 0.6, alpha = 0.6, stroke = 0) +
    # free x per architecture, as S1: AUC-PR spans ~0-1 (oligogenic) but ~0-0.15 (highly polygenic)
    facet_grid(cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP)) +
    scale_y_continuous(breaks = seq_along(SINGLE), labels = rev(set_lab(SINGLE))) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    scale_fill_manual(values = RULE_COL, guide = "none") +
    scale_colour_manual(values = RULE_COL, guide = "none") +
    labs(x = "AUC-PR (causal loci vs neutral SNPs)", y = NULL) + theme_ms()
save_ms("B4_aucpr_methods", pS2, 180, 60, "supp", "B4_aucpr_per_seed.tsv")
num(melt(AUS[, .(rule = label, stratum = arch, n, median, q25, q75)], id.vars = c("rule", "stratum", "n"),
         variable.name = "metric", value.name = "value"),
    "Figure S2 AUC-PR per method (RDA without structure correction)",
    "figures_ssclines_gallery13/B4_aucpr_stats.tsv (= B4_aucpr_per_seed.tsv)")

# =============================================================================
# Figure S3: detection at five cutoffs -- median true and false positives of the six sets
# (remeasure harness, RDA at twice each cutoff, 5 kb window; user, 2026-10-14)
# =============================================================================
message("=== S3 cutoff sensitivity")
pv <- readLines(file.path(REM, "provenance.tsv"))
PV <- setNames(sub("^[^\t]*\t", "", pv), sub("\t.*$", "", pv))
CUTS <- c(top_0.001 = 0.1, top_0.0025 = 0.25, top_0.005 = 0.5, top_0.01 = 1, top_0.02 = 2)  # % per method and predictor
WIN_KB <- 5
stopifnot(PV[["share_mult"]] == "LFMM=1,EMMAX=1,RDA=2", PV[["n_replicates"]] == as.character(NREP),
          PV[["corpus"]] == mvp_arm_label(), setequal(strsplit(PV[["grid"]], ",")[[1]], names(CUTS)),
          as.character(WIN_KB) %in% strsplit(PV[["windows_kb"]], ",")[[1]])
cpm <- fread(file.path(REM, "calls_per_method.tsv"), colClasses = c(seed = "character"))
ccb <- fread(file.path(REM, "calls_combine.tsv"), colClasses = c(seed = "character"))
for (D in list(cpm, ccb))
    stopifnot(setequal(D$seed, PRIM$seed), D[, all(n_called == tp_any_causal + expected_linked + fp_background)])
stopifnot(ccb[combine == "union", uniqueN(n_called), by = .(seed, rung)][, all(V1 == 1L)])   # 1/3: no window
cnt <- function(D, lab) D[, .(seed, rung, label = lab, n = n_called,
                              TP = tp_any_causal + expected_linked, FP = fp_background)]
CUT <- rbind(cnt(ccb[combine == "union" & window_kb == WIN_KB], "1/3 methods"),
             cnt(ccb[combine == "ge2" & window_kb == WIN_KB], "2/3 methods"),
             cnt(ccb[combine == "all3" & window_kb == WIN_KB], "3/3 methods"),
             cpm[, .(seed, rung, label = method, n = n_called,
                     TP = tp_any_causal + expected_linked, FP = fp_background)])
CUT[, cutoff := unname(CUTS[rung])]
stopifnot(!anyNA(CUT$cutoff), setequal(CUT$label, RULES),
          CUT[, .N, by = .(cutoff, label)][, .N == length(CUTS) * N_PANEL && all(N == NREP)])
CM <- CUT[, .(TP = as.numeric(median(TP)), FP = as.numeric(median(FP)), empty = sum(n == 0L), n_rep = .N),
          by = .(cutoff, label)]
# Tie: at the cutoff used throughout (0.25 %) the medians are Figure 2a's (detection_per_seed, the
# panel builder) -- proves the columns and the window are the right ones.
t25 <- merge(CM[cutoff == 0.25, .(label, TP, FP)],
             DET[, .(TPd = median(TP), FPd = median(FP)), by = .(label = as.character(label))], by = "label")
stopifnot(nrow(t25) == N_PANEL, all(t25$TP == t25$TPd), all(t25$FP == t25$FPd))
# The manuscript: the order of the sets by median true positives is the same from 0.25 to 2 % (strict).
ORD_TP <- CM[order(cutoff, -TP), .(order = paste(set_lab(label), collapse = " > "), strict = !anyDuplicated(TP)),
             by = cutoff]
for (i in seq_len(nrow(ORD_TP))) message(sprintf("  S3 cutoff %4.2f %%: %s", ORD_TP$cutoff[i], ORD_TP$order[i]))
k <- ORD_TP[cutoff %in% c(0.25, 0.5, 1, 2)]
if (!all(k$strict) || uniqueN(k$order) != 1L)
    stop("S3: the order of the sets by median true positives is not the same at 0.25, 0.5, 1 and 2 %")
S3_ORDER <- k$order[1]
S3_ORDER_LOW <- ORD_TP[cutoff == min(CUTS), order]
message("tie S3 at 0.25 %: medians = Figure 2a (6 sets, TP and FP); order by median TP identical at 0.25-2 %")
CL <- melt(CM, id.vars = c("cutoff", "label"), measure.vars = c("TP", "FP"), variable.name = "what",
           value.name = "count")
CL[, `:=`(what = factor(c(TP = TP_LAB, FP = FP_LAB)[as.character(what)], levels = c(TP_LAB, FP_LAB)),
          label = factor(label, levels = PANEL_ORDER))]
pS3 <- ggplot(CL, aes(cutoff, count, colour = label)) +
    geom_vline(xintercept = 0.25, linetype = "dashed", colour = MS_REF, linewidth = 0.35) +
    geom_line(linewidth = 0.45) + geom_point(size = 1.2) +
    facet_wrap(vars(what), nrow = 1, scales = "free_y") +
    expand_limits(y = 0) +
    # counts on the square-root scale with even ticks, as Figure 2a (a log axis would drop the 3/3
    # set's zero median false positives at 0.1 %); the cutoff axis is logarithmic
    scale_x_log10(breaks = unname(CUTS), labels = function(x) paste0(x, "%")) +
    scale_y_sqrt(breaks = sqrt_even_breaks) +
    scale_colour_manual(values = RULE_COL, breaks = PANEL_ORDER, labels = SET_LAB, name = NULL) +
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "Cutoff (% of tested SNPs per method and predictor)", y = "Median number of SNPs") + theme_ms()
check_ticks(pS3, "A5_cutoff_sensitivity", axes = "y")
save_ms("A5_cutoff_sensitivity", pS3, 180, 70, "supp",
        "remeasure600/rdaunc_rda2x/{calls_per_method,calls_combine}.tsv")
num(melt(CM[, .(rule = label, stratum = sprintf("cutoff %s%%", cutoff), n = n_rep, TP_median = as.numeric(TP),
                FP_median = as.numeric(FP), n_empty = as.numeric(empty))],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value"),
    "Figure S3 detection at five cutoffs (5 kb window; RDA at twice each cutoff)",
    "remeasure600/rdaunc_rda2x/{calls_per_method,calls_combine}.tsv")

# =============================================================================
# Figure S8: % of replicates in which each set is the least accurate (median of 3 engines, ties
# split evenly) by design cell and by demographic scenario, + Kruskal-Wallis tests on the
# within-replicate rank (user, 2026-10-14)
# =============================================================================
message("=== S8 least accurate set by design cell and demographic scenario")
COV <- unique(DET[, .(seed, arch, pleio, regime, block)])
stopifnot(nrow(COV) == NREP, !anyNA(COV))
# rank_long names the demography by its block tag, the gallery by its label: one to one
TAGB <- unique(merge(RL[, .(seed, tag = block)], COV[, .(seed, block)], by = "seed")[, .(tag, block)])
stopifnot(nrow(TAGB) == 5L, uniqueN(TAGB$tag) == 5L, uniqueN(TAGB$block) == 5L)
M3 <- merge(RL[engine == "median of 3", .(seed, panel, rank)], COV, by = "seed")
stopifnot(nrow(M3) == NREP * N_PANEL)
I3 <- merge(rank_indicators(rank_split(M3, "seed"), "seed", n_set = N_PANEL), COV, by = "seed")
cell_of <- function(a, p, r) paste(a, p, r, sep = ", ")
strat <- function(cols, kind) I3[, .(least_pct = 100 * mean(bottom1), last2_pct = 100 * mean(bottom2), n = .N),
                                 by = c("panel", cols)][, kind := kind][]
S8A <- strat("arch", "genetic architecture")
S8B <- strat("block", "demographic scenario")
S8C <- strat(c("arch", "pleio", "regime"), "design cell")
stopifnot(nrow(S8A) == 3L * N_PANEL, nrow(S8B) == 5L * N_PANEL, nrow(S8C) == 12L * N_PANEL,
          all(S8A$n == NREP / 3), all(S8B$n == NREP / 5), all(S8C$n == NREP / 12))
# Tie: the architecture and scenario values are rank_summary.tsv's (mvp_panel_rank.R, the same tie rule).
t8 <- rbind(S8A[, .(panel, stratum = as.character(arch), bottom1 = least_pct, bottom2 = last2_pct)],
            merge(S8B, TAGB, by = "block")[, .(panel, stratum = tag, bottom1 = least_pct, bottom2 = last2_pct)])
t8 <- melt(t8, id.vars = c("panel", "stratum"), variable.name = "metric", value.name = "v")[, metric := as.character(metric)]
t8 <- merge(t8, RS[engine == "median of 3" & metric %in% c("bottom1", "bottom2"), .(panel, stratum, metric, value)],
            by = c("panel", "stratum", "metric"))
stopifnot(nrow(t8) == 8L * N_PANEL * 2L, max(abs(t8$v - t8$value)) < 1e-9)
message(sprintf("tie S8: %d architecture / scenario values = rank_summary.tsv, max |diff| %.1e",
                nrow(t8), max(abs(t8$v - t8$value))))
KWT <- rbindlist(lapply(c("design cell", "demographic scenario"), function(f) rbindlist(lapply(PANEL_ORDER, function(s) {
    D <- M3[panel == s]
    g <- if (f == "design cell") factor(cell_of(D$arch, D$pleio, D$regime)) else factor(D$block)
    kt <- kruskal.test(D$rank, g)
    data.table(factor = f, panel = s, groups = nlevels(g), n = nrow(D), H = unname(kt$statistic),
               df = unname(kt$parameter), p = kt$p.value)
}))))
KWT[, p_holm := p.adjust(p, method = "holm")]
stopifnot(nrow(KWT) == 12L, all(KWT[factor == "design cell", df] == 11), all(KWT[factor == "demographic scenario", df] == 4))
print(KWT[, .(factor, set = set_lab(panel), H = round(H, 2), df, p = signif(p, 3), p_holm = signif(p_holm, 3))])
ALLS <- rbind(S8A[, .(kind, stratum = as.character(arch), panel, least_pct, last2_pct, n)],
              S8B[, .(kind, stratum = as.character(block), panel, least_pct, last2_pct, n)],
              S8C[, .(kind, stratum = cell_of(arch, pleio, regime), panel, least_pct, last2_pct, n)])
MX23 <- ALLS[panel == "2/3 methods"][order(-least_pct)][1]
message(sprintf("2/3 set: least accurate in at most %.1f %% of the replicates of any of %d strata (%s: %s, n = %d)",
                MX23$least_pct, ALLS[panel == "2/3 methods", .N], MX23$kind, MX23$stratum, MX23$n))
if (MX23$least_pct > 6 + 1e-9) stop("S8: the 2/3 set is least accurate in more than 6 % of a stratum (the manuscript says at most 6 %)")
# Heat maps: white at chance (1 in 6), teal below, red above; value printed in each tile
CHANCE <- 100 / N_PANEL
HEAT_MAX <- max(20, 10 * ceiling(max(S8B$least_pct, S8C$least_pct) / 10))
HEAT_COL <- c(MINOU[["teal"]], "#FFFFFF", MINOU[["red"]])
HEAT_AT  <- scales::rescale(c(0, CHANCE, HEAT_MAX))
heat_fill <- function(v) scales::gradient_n_pal(HEAT_COL, HEAT_AT)(scales::rescale(v, from = c(0, HEAT_MAX)))
heat <- function(D, x, facet_layer = NULL, legend = TRUE) {
    D <- copy(D)[, `:=`(panel = factor(panel, levels = rev(PANEL_ORDER)), tone = ms_text_on(heat_fill(least_pct)))]
    ggplot(D, aes(.data[[x]], panel, fill = least_pct)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = sprintf("%.1f", least_pct), colour = tone), size = MS_LAB) +
        facet_layer +
        scale_fill_gradientn(colours = HEAT_COL, values = HEAT_AT, limits = c(0, HEAT_MAX),
                             breaks = seq(0, HEAT_MAX, 10), labels = function(b) paste0(b, "%"),
                             name = "Replicates in which the set is least accurate (white: chance, 16.7%)",
                             guide = if (legend) guide_colourbar(theme = theme(
                                 legend.key.width = grid::unit(45, "mm"), legend.key.height = grid::unit(2.5, "mm"),
                                 legend.title.position = "top")) else "none") +
        scale_colour_identity() +
        scale_y_discrete(labels = SET_LAB) +
        labs(x = NULL, y = NULL) + theme_ms() +
        theme(axis.line = element_blank(), axis.ticks = element_blank())
}
pS8a <- heat(S8C, "regime", facet_grid(cols = vars(arch, pleio), labeller = labeller(arch = ARCH_WRAP)), legend = FALSE)
pS8b <- heat(S8B, "block")
for (pp in list(a = pS8a, b = pS8b)) vocab_check(fig_text(pp), "R3_least_accurate_strata part (figure text)")
pS8 <- cowplot::plot_grid(pS8a, pS8b, ncol = 1, rel_heights = c(56, 54), labels = c("a", "b"),
                          label_size = 8, label_fontface = "bold", label_fontfamily = theme_ms()$text$family,
                          align = "v", axis = "lr")
save_ms("R3_least_accurate_strata", pS8, 180, 110, "supp",
        "rank_long.tsv (median of 3) + detection_per_seed.tsv (design cell, demographic scenario)")
num(rbind(ALLS[, .(rule = panel, stratum = paste0(kind, ": ", stratum), n, metric = "least_accurate_pct", value = least_pct)],
          ALLS[, .(rule = panel, stratum = paste0(kind, ": ", stratum), n, metric = "last_two_pct", value = last2_pct)])[
              , engine := "median of 3"],
    "Figure S8 / Table S1 (c) % replicates least accurate by stratum (ties split evenly)",
    "figures_ssclines_offset13/rank/rank_long.tsv + figures_ssclines_gallery13/detection_per_seed.tsv")
num(melt(KWT[, .(rule = panel, stratum = factor, n, groups = as.numeric(groups), H, df = as.numeric(df), p, p_holm)],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value")[, engine := "median of 3"],
    "Figure S8 / Table S1 (e) Kruskal-Wallis test of the within-replicate rank across strata (Holm over 12)",
    "figures_ssclines_offset13/rank/rank_long.tsv + figures_ssclines_gallery13/detection_per_seed.tsv")

# =============================================================================
# Table S1 (user, 2026-10-14; replaces the ad hoc python behind the first version, never kept)
#   (a) % least / two least / most / two most accurate per set x engine: ties split evenly, an
#       empty set last; 95 % percentile bootstrap CI (replicates resampled jointly, S1_B) and a
#       two-sided permutation P against chance (set labels shuffled within each replicate, S1_B,
#       P = (k + 1) / (S1_B + 1)); set.seed(S1_SEED) before each of the two
#   (b) delta accuracy (set - causal loci): median, min and max over replicates, % delta >= 0
#   (c) range of the 2/3 set's % least / two least accurate over the strata of S8 (+ architecture)
#   (d) variance of accuracy explained by replicate and by set (reliability_variance.tsv)
#   (e) the S8 Kruskal-Wallis tests
# =============================================================================
message("=== Table S1")
A_MET <- c(bottom1 = "least accurate", bottom2 = "among the two least accurate",
           top1 = "most accurate", top2 = "among the two most accurate")
A_CHANCE <- 100 * c(bottom1 = 1, bottom2 = 2, top1 = 1, top2 = 2) / N_PANEL
ENGS <- c("median of 3", "gradientForest", "LFMM2", "RDA")
stopifnot(setequal(ENGS, names(ENG_LAB)))
TI <- rank_indicators(RX, c("seed", "engine"), n_set = N_PANEL)
SEEDS <- sort(unique(TI$seed), method = "radix")
NS <- length(SEEDS); NE <- length(ENGS); NM <- length(A_MET)
stopifnot(NS == NREP)
MAT <- lapply(setNames(ENGS, ENGS), function(e) lapply(setNames(names(A_MET), names(A_MET)), function(m) {
    w <- dcast(TI[engine == e], seed ~ panel, value.var = m)
    stopifnot(identical(w$seed, SEEDS), !anyNA(w))
    as.matrix(w[, RULES, with = FALSE])
}))
OBS <- rbindlist(lapply(ENGS, function(e) rbindlist(lapply(names(A_MET), function(m)
    data.table(engine = e, metric = m, panel = RULES, value = 100 * colMeans(MAT[[e]][[m]]))))))
# Tie: the point estimates are rank_summary.tsv's (mvp_panel_rank.R: the same rule, other code).
x <- merge(OBS, RS[stratum == "pooled", .(engine, panel, metric, rs = value)], by = c("engine", "panel", "metric"))
stopifnot(nrow(x) == nrow(OBS), max(abs(x$value - x$rs)) < 1e-9)
set.seed(S1_SEED)
BOOT <- array(NA_real_, c(S1_B, NE, NM, N_PANEL))
for (b in seq_len(S1_B)) {
    i <- sample.int(NS, NS, replace = TRUE)
    for (ei in seq_len(NE)) for (mi in seq_len(NM)) BOOT[b, ei, mi, ] <- colMeans(MAT[[ei]][[mi]][i, , drop = FALSE])
}
set.seed(S1_SEED)
PERM <- array(NA_real_, c(S1_B, NE, NM, N_PANEL))
BLK <- rep(seq_len(NS), each = N_PANEL); ROWS <- rep(seq_len(NS), N_PANEL)
for (b in seq_len(S1_B)) {
    o <- order(BLK, runif(NS * N_PANEL))
    P <- matrix((o - 1L) %% N_PANEL + 1L, nrow = NS, byrow = TRUE)   # row i: a random order of the six sets
    ix <- cbind(ROWS, as.vector(P))                                 # set s of replicate i takes slot P[i, s]
    for (ei in seq_len(NE)) for (mi in seq_len(NM)) PERM[b, ei, mi, ] <- colMeans(matrix(MAT[[ei]][[mi]][ix], nrow = NS))
}
TA <- OBS[, {
    ei <- match(engine, ENGS); mi <- match(metric, names(A_MET)); si <- match(panel, RULES)
    bt <- 100 * BOOT[, ei, mi, si]; pm <- 100 * PERM[, ei, mi, si]; ch <- A_CHANCE[[metric]]
    kk <- sum(abs(pm - ch) >= abs(value - ch) - 1e-9)
    .(lo = quantile(bt, 0.025, names = FALSE), hi = quantile(bt, 0.975, names = FALSE), chance = ch,
      p_perm = (kk + 1) / (S1_B + 1), k_perm = kk, null_mean = mean(pm))
}, by = .(engine, metric, panel, value)]
stopifnot(nrow(TA) == NE * NM * N_PANEL, max(abs(TA$null_mean - TA$chance)) < 0.5)   # null centred on chance
# Check against the first Table S1 (ad hoc python, median of 3). It counted a tie fully on the
# inclusive side, so a value may move by half a replicate for every tied unit whose tie straddles
# that cutoff for that set; compared in counts (out of 600), the quoted % being rounded.
OLD_PY <- data.table(panel  = c("2/3 methods", "2/3 methods", "2/3 methods", "2/3 methods", "1/3 methods"),
                     metric = c("bottom1", "bottom2", "top1", "top2", "bottom1"),
                     old_pct = c(1.83, 11.17, 4.33, 18.83, 3.17),
                     old_lo = c(0.83, 8.67, NA, NA, 1.83), old_hi = c(3.00, 13.67, NA, NA, 4.67))
in_region <- function(m, r) switch(m, bottom1 = r >= N_PANEL, bottom2 = r >= N_PANEL - 1L,
                                   top1 = r <= 1L, top2 = r <= 2L)
T3 <- RX[engine == "median of 3" & tie_n > 1L]
OC <- rbindlist(lapply(seq_len(nrow(OLD_PY)), function(i) {
    m <- OLD_PY$metric[i]; s <- OLD_PY$panel[i]
    st <- T3[panel == s, .(inside = any(in_region(m, r)) && !all(in_region(m, r))), by = seed][inside == TRUE, seed]
    nw <- TA[engine == "median of 3" & metric == m & panel == s]
    data.table(OLD_PY[i], new_pct = nw$value, new_lo = nw$lo, new_hi = nw$hi,
               old_count = round(OLD_PY$old_pct[i] * NREP / 100), new_count = nw$value * NREP / 100,
               tied_units = length(st), tied_seeds = paste(st, collapse = " "))
}))
OC[, within := abs(new_count - old_count) <= 0.5 * tied_units + 1e-9]
message("Table S1 (a) vs the first version (ad hoc python; median of 3; counts out of ", NREP, "):")
print(OC[, .(set = set_lab(panel), metric, old_pct, new_pct = round(new_pct, 3), old_count, new_count, tied_units,
             tied_seeds, within, old_ci = fifelse(is.na(old_lo), "", sprintf("%.2f-%.2f", old_lo, old_hi)),
             new_ci = sprintf("%.2f-%.2f", new_lo, new_hi))])
if (!all(OC$within)) stop("Table S1 (a): a value moved by more than its tied units allow (table above)")
# (b) delta accuracy per set, from the per-replicate table; tied to the gallery's C_delta_stats
DB <- OFF[, .(n = .N, median = median(delta), min = min(delta), max = max(delta),
              ge0_pct = 100 * mean(delta >= 0), gt0_pct = 100 * mean(delta > 0), zero = sum(delta == 0)),
          by = .(panel = as.character(label))]
x <- merge(DB, CS[arch == "pooled", .(panel = label, cs_n = n, cs_median = median, reach_pct, win_pct)], by = "panel")
stopifnot(nrow(x) == N_PANEL, all(x$n == x$cs_n), max(abs(x$median - x$cs_median)) < 1e-12,
          max(abs(x$ge0_pct - x$reach_pct)) < 1e-9, max(abs(x$gt0_pct - x$win_pct)) < 1e-9)
SPREAD <- DB[, .(abs_lo = min(abs(median)), abs_hi = max(abs(median)), ge0_lo = min(ge0_pct), ge0_hi = max(ge0_pct))]
half_up <- function(v, d) floor(v * 10^d + 0.5) / 10^d
if (!isTRUE(all.equal(c(half_up(SPREAD$abs_lo, 3), half_up(SPREAD$abs_hi, 3), half_up(SPREAD$ge0_lo, 0),
                        half_up(SPREAD$ge0_hi, 0)), c(0.026, 0.040, 26, 38))))
    stop("Table S1 (b): the spread across sets no longer reads 0.026-0.040 and 26-38 % (the manuscript's numbers)")
message(sprintf("Table S1 (b): median delta %.3f to %.3f over the six sets; delta >= 0 in %.1f-%.1f %% (zero deltas: %d)",
                -SPREAD$abs_hi, -SPREAD$abs_lo, SPREAD$ge0_lo, SPREAD$ge0_hi, sum(DB$zero)))
# (c) the 2/3 set over strata, from the S8 values; (d) the variance partition; (e) the S8 tests
TC <- ALLS[panel == "2/3 methods", .(n_strata = .N, n_per_stratum = n[1], least_lo = min(least_pct),
                                     least_hi = max(least_pct), last2_lo = min(last2_pct), last2_hi = max(last2_pct)),
           by = kind]
stopifnot(nrow(TC) == 3L, ALLS[panel == "2/3 methods", uniqueN(n), by = kind][, all(V1 == 1L)])
VP_LAB <- c(replicate = "replicate", panel = "candidate set",
            `replicate x panel (residual)` = "replicate × candidate set (residual)")
stopifnot(setequal(RVP$source, names(VP_LAB)), uniqueN(RVP$n_rep) == 1L)
TS1 <- rbindlist(list(
    TA[, .(section = "a", set = set_lab(panel), engine = ENG_LAB[engine], stratum = "pooled",
           metric = paste("% replicates", A_MET[metric]), value, lo, hi, chance, p = p_perm, n = NREP)],
    melt(DB[, .(section = "b", set = set_lab(panel), engine = "median of 3", stratum = "pooled", n,
                `median delta accuracy` = median, `min delta accuracy` = min, `max delta accuracy` = max,
                `% replicates delta >= 0` = ge0_pct)],
         id.vars = c("section", "set", "engine", "stratum", "n"), variable.name = "metric", value.name = "value"),
    data.table(section = "b", set = "spread over the six sets", engine = "median of 3", stratum = "pooled",
               metric = c("|median delta accuracy|, min", "|median delta accuracy|, max",
                          "% replicates delta >= 0, min", "% replicates delta >= 0, max"),
               value = unlist(SPREAD[, .(abs_lo, abs_hi, ge0_lo, ge0_hi)], use.names = FALSE)),
    melt(TC[, .(section = "c", set = "2/3 set", engine = "median of 3", stratum = kind, n = n_strata,
                `% replicates least accurate, min` = least_lo, `% replicates least accurate, max` = least_hi,
                `% replicates among the two least accurate, min` = last2_lo,
                `% replicates among the two least accurate, max` = last2_hi)],
         id.vars = c("section", "set", "engine", "stratum", "n"), variable.name = "metric", value.name = "value"),
    RVP[, .(section = "d", set = "all six sets", engine = "median of 3", stratum = "pooled", n = n_rep,
            metric = paste("% of the variance of accuracy explained by", VP_LAB[source]), value = share_pct)],
    RVP[1, .(section = "d", set = "all six sets", engine = "median of 3", stratum = "pooled", n = n_rep,
             metric = c("spread of the six set medians (max - min)", "median spread of the six sets within a replicate"),
             value = c(panel_median_range, within_rep_range_median))],
    KWT[, .(section = "e", set = set_lab(panel), engine = "median of 3", stratum = factor, n,
           metric = "Kruskal-Wallis H", value = H, df, p, p_holm)]), use.names = TRUE, fill = TRUE)
TS1[, metric := as.character(metric)]
fwrite(TS1, file.path(OUT, "Table_S1.tsv"), sep = "\t")
message(sprintf("  Table_S1.tsv %d rows", nrow(TS1)))
num(TA[, .(rule = panel, stratum = "pooled", engine = ENG_LAB[engine], metric = paste0(metric, "_pct"), value, lo, hi,
           expected = chance, n = NREP, p_perm)],
    sprintf("Table S1 (a) within-replicate rank: %% (95%% bootstrap CI, %s resamples) and permutation P (%s) -- the text quotes these",
            format(S1_B, big.mark = ","), format(S1_B, big.mark = ",")),
    sprintf("figures_ssclines_offset13/rank/rank_long.tsv (ties split evenly; seed %d)", S1_SEED))
num(OC[, .(rule = panel, stratum = "pooled", engine = "median of 3", metric = paste0(metric, "_pct_first_version"),
           value = old_pct, lo = old_lo, hi = old_hi, n = tied_units)],
    "tie rule: the first Table S1 (ad hoc python, ties counted on the inclusive side); n = tied units inside the cutoff",
    "manuscript Table S1 v1 (values as quoted 2026-10-14)")
num(rbind(melt(DB[, .(rule = panel, stratum = "pooled", n, median_delta = median, min_delta = min, max_delta = max,
                      delta_ge0_pct = ge0_pct)], id.vars = c("rule", "stratum", "n"),
               variable.name = "metric", value.name = "value"),
          data.table(rule = "all six sets", stratum = "pooled", n = NREP,
                     metric = c("abs_median_delta_min", "abs_median_delta_max", "delta_ge0_pct_min", "delta_ge0_pct_max"),
                     value = unlist(SPREAD, use.names = FALSE)))[, engine := "median of 3"],
    "Table S1 (b) delta accuracy vs causal loci (delta >= 0 = C_delta_stats reach_pct)",
    "figures_ssclines_gallery13/offset_delta_per_seed.tsv (= C_delta_stats.tsv)")
num(melt(TC[, .(rule = "2/3 methods", stratum = kind, n = n_strata, least_pct_min = least_lo, least_pct_max = least_hi,
                last_two_pct_min = last2_lo, last_two_pct_max = last2_hi)],
         id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value")[, engine := "median of 3"],
    "Table S1 (c) 2/3 set: range over strata", "Figure S8 values (rank_long.tsv + detection_per_seed.tsv)")

# =============================================================================
# Data S2: every TSV a display item is drawn from, at its path under benchmarks/mvp_eval/, and a
# manifest mapping file -> display items (user, 2026-10-14)
# =============================================================================
message("=== Data S2")
rel_eval <- function(f) { stopifnot(startsWith(f, paste0(EVAL, "/"))); substring(f, nchar(EVAL) + 2L) }
DS2 <- rbindlist(list(
    list(file.path(GAL, "detection_per_seed.tsv"),
         "Figure 2a; Figure S4; Figure S5; Figure S3 (check at the 0.25% cutoff); Figure S8 (design cell and demographic scenario of each replicate)",
         "per replicate x candidate set: SNPs identified (`n`), causal (`n_causal`) and linked (`n_linked`) SNPs; true positives `TP` = causal + linked, false positives `FP` = neutral SNPs; design factors `arch`, `pleio`, `regime`, `block`; `label` / `panel` = candidate set (`1/3 methods`, `2/3 methods`, `3/3 methods` = the 1/3, 2/3 and 3/3 sets)",
         "benchmarks/mvp_j16_gallery.R"),
    list(file.path(RNK, "rank_long.tsv"),
         "Figure 2b; Figure S7; Figure S8; Table S1 (a, c, e)",
         "per replicate x offset engine x candidate set: accuracy (-Kendall tau), accuracy of the causal loci (`oracle`), `delta`, within-replicate `rank` (tied sets get the average rank; every count splits it evenly), `empty` = no SNP identified (ranked last); `label` / `panel` = candidate set (`1/3 methods`, `2/3 methods`, `3/3 methods` = the 1/3, 2/3 and 3/3 sets)",
         "benchmarks/mvp_panel_rank.R"),
    list(file.path(RNK, "rank_summary.tsv"),
         "order of the candidate sets in Figure 2b, Figure S2, Figure S3, Figure S6, Figure S7, Figure S8; Table S1 (a, check)",
         "per engine x candidate set x stratum: % replicates at the top / bottom k ranks (ties split evenly) and mean rank, 95% bootstrap CI over replicates (2,000 resamples); `label` / `panel` = candidate set (`1/3 methods`, `2/3 methods`, `3/3 methods` = the 1/3, 2/3 and 3/3 sets)",
         "benchmarks/mvp_panel_rank.R"),
    list(file.path(GAL, "offset_delta_per_seed.tsv"),
         "Figure S6; Table S1 (b)",
         "per replicate x candidate set: offset accuracy (-Kendall tau, median over 100 landscape gardens and three offset engines), accuracy of the causal loci (`oracle`), `delta`; `label` / `panel` = candidate set (`1/3 methods`, `2/3 methods`, `3/3 methods` = the 1/3, 2/3 and 3/3 sets)",
         "benchmarks/mvp_j16_gallery.R"),
    list(file.path(GAL, "C_delta_stats.tsv"),
         "Table S1 (b)",
         "delta summary per candidate set x architecture: median, bootstrap CI, IQR, `win_pct` (delta > 0), `reach_pct` (delta >= 0), Wilcoxon P; `label` / `panel` = candidate set (`1/3 methods`, `2/3 methods`, `3/3 methods` = the 1/3, 2/3 and 3/3 sets)",
         "benchmarks/mvp_j16_gallery.R"),
    list(file.path(GAL, "B4_aucpr_per_seed.tsv"),
         "Figure S2",
         "per replicate x method: AUC-PR of the causal loci (neutral SNPs the only false positives, linked SNPs excluded); RDA without structure correction",
         "benchmarks/mvp_j16_gallery.R"),
    list(file.path(REM, "calls_per_method.tsv"),
         "Figure S3 (single methods)",
         "per replicate x method x cutoff (`rung`; RDA at twice each): SNPs identified (`n_called`), causal (`tp_any_causal`), linked (`expected_linked`) and neutral (`fp_background`) SNPs",
         "benchmarks/mvp_remeasure_600.R"),
    list(file.path(REM, "calls_combine.tsv"),
         "Figure S3 (agreement sets)",
         "per replicate x agreement set (`union` = 1/3 set, `ge2` = 2/3 set, `all3` = 3/3 set) x cutoff x window (kb): the same counts; Figure S3 uses the 5 kb window",
         "benchmarks/mvp_remeasure_600.R"),
    list(file.path(EVAL, "remeasure600/rdaunc/rank_metrics.tsv"),
         "Figure S1 (RDA without structure correction)",
         "per replicate x method: AUC-PR and R-precision; RDA without structure correction (its LFMM and EMMAX rows equal the corrected arm's)",
         "benchmarks/mvp_remeasure_600.R"),
    list(file.path(EVAL, "remeasure600/onefit/rank_metrics.tsv"),
         "Figure S1 (RDA with structure correction; LFMM)",
         "per replicate x method: AUC-PR and R-precision; RDA with structure correction",
         "benchmarks/mvp_remeasure_600.R"),
    list(file.path(EVAL, "remeasure600/onefit/covariates.tsv"),
         "Figure S1 (genetic architecture of each replicate)",
         "per replicate: genetic architecture (`arch_level`), design factors and population-structure covariates",
         "benchmarks/mvp_remeasure_600.R"),
    list(file.path(RNK, "reliability_variance.tsv"),
         "Table S1 (d)",
         "sums of squares of accuracy (median of the three engines) by replicate, candidate set (`panel`) and residual; replicates with all six sets scored",
         "benchmarks/mvp_rank_reliability.R")))
setnames(DS2, c("src", "display_items", "description", "written_by"))
stopifnot(all(file.exists(DS2$src)))
DS2[, file := vapply(src, rel_eval, "")]
DS2[, c("rows", "columns") := {
    d <- lapply(src, function(f) fread(f, colClasses = "character"))
    list(vapply(d, nrow, 1L), vapply(d, ncol, 1L))
}]
for (it in c("Figure 2a", "Figure 2b", paste0("Figure S", 1:8)))
    if (!any(grepl(it, DS2$display_items, fixed = TRUE))) stop("Data S2: no source listed for ", it)
for (sec in letters[1:5])
    if (!any(grepl(paste0("Table S1 \\([^)]*", sec), DS2$display_items))) stop("Data S2: no source for Table S1 (", sec, ")")
MAN <- DS2[, .(file, display_items, rows, columns, description, written_by)]
STAGE <- file.path(tempdir(), "Data_S2")
dir.create(STAGE)
for (i in seq_len(nrow(DS2))) {
    dir.create(file.path(STAGE, dirname(DS2$file[i])), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(DS2$src[i], file.path(STAGE, DS2$file[i])))
}
fwrite(MAN, file.path(STAGE, "manifest.tsv"), sep = "\t")
fwrite(MAN, file.path(OUT, "Data_S2_manifest.tsv"), sep = "\t")
DS2_ZIP <- "Data_S2.zip"
old <- setwd(STAGE)
unlink(file.path(OUT, DS2_ZIP))
stopifnot(utils::zip(file.path(OUT, DS2_ZIP), c("manifest.tsv", DS2$file), flags = "-9Xq") == 0L)
setwd(old)
message(sprintf("  %s: %d source tables + manifest.tsv (%.1f MB)", DS2_ZIP, nrow(DS2),
                file.size(file.path(OUT, DS2_ZIP)) / 1e6))

# =============================================================================
# Captions (all explanation lives here, not on the figures)
# =============================================================================
nP <- function(rule) as.integer(CS[label == rule & arch == "pooled", n])
EMPTY_CLAUSE <- if (nrow(EMPTY)) paste0(
    paste(sprintf("The %s is empty in %d replicate%s", set_lab(EMPTY$label), EMPTY$N,
                  ifelse(EMPTY$N == 1, "", "s")), collapse = "; "), ".") else ""
SETS_DEF <- paste0(
    "Candidate sets: the SNPs identified by each of three single methods (LFMM, EMMAX and RDA without structure ",
    "correction) and three agreement sets built from them. The cutoff, fixed for all replicates: each method keeps ",
    "the 0.25% of tested SNPs with the lowest p-value for each of the two environmental predictors (LFMM and EMMAX ",
    "test the predictors separately and their two lists are pooled; RDA tests both in one multivariate p-value and ",
    "keeps 0.5%), so each method identifies about 0.5% of the tested SNPs; ",
    sprintf("the cutoff equals the median proportion of causal loci among tested SNPs over the %d replicates (%.2f%%). ", NREP, CS_MED),
    "1/3 set: every SNP identified by at least one method; 2/3 and 3/3 sets: SNPs within 5 kb of SNPs identified ",
    "by at least two / all three methods (each method keeps its own SNPs, so the 3/3 set can be larger than any ",
    "single-method set).")
SENS_DEF <- paste0("Sensitivity is the number of causal or linked SNPs a set identifies; precision is the proportion ",
                   "of causal or linked SNPs among those identified. ")
TPFP_DEF <- paste0(
    "True positives are causal or linked SNPs (causal loci and the SNPs linked to them); false positives are neutral ",
    "SNPs (unlinked to any causal locus). ")
PLANE_DEF <- paste0(
    "Both axes use a square-root scale with evenly spaced ticks; each plot has its own false-positive axis. ",
    "Outlines (bags) are the convex hull of each candidate set's median and the half of its replicates closest to it ",
    "(distance on the square-root axes, each scaled by its interquartile range); labelled points are the medians (true ",
    "and false positives medianed separately). ",
    SETS_DEF, " ", EMPTY_CLAUSE, if (nzchar(EMPTY_CLAUSE)) " " else "",
    "The causal loci themselves are not drawn: by construction they contain no linked SNPs.")
DETECTION_DEF <- paste0(TPFP_DEF, PLANE_DEF)
ACC_DEF <- paste0(
    "Accuracy is −Kendall τ between the predicted genetic offset and the fitness of the 100 ",
    "source populations in a common garden, median over the 100 landscape gardens, for each of three ",
    "offset engines (Gradient Forest, the LFMM2 geometric offset and RDA without structure correction).")
ORDER_DEF <- paste0(
    "Candidate sets are ordered by the percentage of replicates in which they are the least accurate under the ",
    "median of the three engines, fewest first (ties broken by the percentage in which they are among the two least ",
    "accurate).")
RANK_DEF <- paste0(
    "Within each replicate the six candidate sets are ranked by offset accuracy (1 = most accurate), ",
    "separately for each offset engine and for the median of the three. ", ACC_DEF,
    " The causal loci are a constant reference within a replicate, so ranks on accuracy equal ranks on accuracy relative to the causal loci. ",
    sprintf("Under random ordering every set would sit at each rank in %.1f%% of replicates. ", 100 / N_PANEL),
    sprintf("Tied sets split the tied ranks evenly: a set tied for ranks k and k + 1 counts one half at each (%d of %s replicate × engine units carry a tie, %s). ",
            N_TIED, format(N_UNITS, big.mark = ","),
            if (N_TIE3 == 0L) "all between two sets" else sprintf("%d of them among three or more sets", N_TIE3)),
    if (nrow(EMPTY)) "An empty set has no offset and is ranked last. " else "",
    ORDER_DEF)
ACC3_DEF <- paste0(
    "Accuracy is −Kendall τ between the predicted genetic offset and the fitness of the 100 source ",
    "populations in a common garden, median over the 100 landscape gardens and then over three offset ",
    "engines (Gradient Forest, the LFMM2 geometric offset and RDA without structure correction). ",
    "Half-violins: density; boxes: median and interquartile range; dots: replicates. ",
    "The grey track below the line is the same accuracy computed from the causal loci of each replicate (reference). ",
    if (nrow(EMPTY)) sprintf("Empty sets have no offset and are omitted (%s). ",
                             paste(sprintf("%s: n = %d of %d", set_lab(EMPTY$label),
                                           vapply(EMPTY$label, nP, 1L), NREP), collapse = "; ")) else "",
    "Candidate sets ordered as in Figure 2b.")
AUC_DEF <- paste0(
    "AUC-PR: the area under the precision–recall curve of the SNPs ranked by p-value, with the causal loci as ",
    "positives and neutral SNPs (unlinked to any causal locus) as the only negatives; SNPs linked to a causal locus ",
    "are excluded. It is computed as the mean, over the causal loci, of the proportion of causal loci among the SNPs ",
    "ranked up to each of them (causal loci never retrieved count 0). ")
fmt_p <- function(p) if (p < 1e-4) sprintf("%.1e", p) else formatC(p, format = "g", digits = 2)
kw_text <- function(f) {
    K <- KWT[factor == f]
    sprintf("across the %d %ss (df = %d): %s", K$groups[1], f, K$df[1],
            paste(sprintf("%s H = %.1f, P = %s (adjusted %s)", set_lab(K$panel), K$H,
                          vapply(K$p, fmt_p, ""), vapply(K$p_holm, fmt_p, "")), collapse = "; "))
}
EM3 <- CM[empty > 0L][order(match(label, PANEL_ORDER), cutoff)]
S3_EMPTY <- if (nrow(EM3)) paste(EM3[, .(t = paste(sprintf("%d replicate%s at %s%%", empty, ifelse(empty == 1L, "", "s"),
                                                            cutoff), collapse = ", ")), by = label][
                                     , sprintf("%s: %s", set_lab(label), t)], collapse = "; ") else "none"
CAP <- list(
    Fig2_simulation = paste0(
        sprintf("Detection and offset accuracy of the six candidate sets over %d simulated replicates. ", NREP),
        sprintf("(a) Detection of adaptive loci by genetic architecture (n = %d replicates per level): true against false positives of each replicate. ", NREP / 3),
        TPFP_DEF, SENS_DEF, PLANE_DEF, median_only(DET, "arch"),
        sprintf(" (b) Rank of each candidate set among the six within a replicate, per offset engine (n = %d replicates); bars: percentage of replicates at each rank. ", NREP),
        RANK_DEF),
    S_rda_correction = paste0(
        "Why the candidate sets use RDA without structure correction. ",
        sprintf("Ranking of the causal loci by RDA with and without structure correction, by genetic architecture (n = %d replicates per level), LFMM for reference. ", NREP / 3),
        sprintf("Without correction RDA ranks the causal loci better in %s of %s replicates (median AUC-PR %s vs %s, %s-fold). ",
                RDA_FACTS$wins, RDA_FACTS$n, RDA_FACTS$unc, RDA_FACTS$cor, RDA_FACTS$ratio),
        sprintf("On this landscape population structure follows the climate gradients: the conditioning principal components are collinear with the predictors (largest variance inflation factor ≥ 10 in %s of %s replicates), and the corrected model keeps a median %s of the uncorrected model's adjusted R². ",
                RDA_FACTS$vif10, RDA_FACTS$n, sub(" %", "%", RDA_FACTS$kept, fixed = TRUE)),
        AUC_DEF,
        "With correction: partial RDA conditioned on the first K principal components (K = number of ancestral populations, at least 2); ",
        "without correction: the same model with no conditioning. ",
        "Half-violins: density; boxes: median and interquartile range; dots: replicates. ",
        "Each architecture has its own x-axis range; compare within an architecture."),
    B4_aucpr_methods = paste0(
        sprintf("Ranking of the causal loci by each single method, by genetic architecture (n = %d replicates per level): LFMM, RDA without structure correction and EMMAX. ", NREP / 3),
        sprintf("Median AUC-PR over the %d replicates: %s. ", NREP,
                paste(sprintf("%s %.3f", set_lab(SINGLE), AUM[arch == "pooled"][match(SINGLE, label), m]), collapse = ", ")),
        AUC_DEF,
        "Half-violins: density; boxes: median and interquartile range; dots: replicates. ",
        "Each architecture has its own x-axis range; compare within an architecture. ",
        "Methods ordered as the candidate sets in Figure 2b."),
    A5_cutoff_sensitivity = paste0(
        "Detection at five cutoffs. ",
        sprintf("Median true positives (causal or linked SNPs; left) and false positives (neutral SNPs; right) of the six candidate sets over the %d replicates, ", NREP),
        "at cutoffs of 0.1, 0.25, 0.5, 1 and 2% of the tested SNPs per method and predictor (RDA without structure correction, which tests both predictors in one p-value, at twice each), 5 kb agreement window. ",
        "Dashed line: the 0.25% cutoff used throughout (the medians there are those of Figure 2a). ",
        sprintf("The order of the candidate sets by median true positives is the same at every cutoff from 0.25%% to 2%% (%s); at 0.1%% it is %s. ",
                S3_ORDER, S3_ORDER_LOW),
        sprintf("Medians include empty sets (%s). ", S3_EMPTY),
        SENS_DEF,
        "Both count axes use a square-root scale with evenly spaced ticks; the cutoff axis is logarithmic."),
    D5_detection_cells = paste0(
        sprintf("Detection (as in Figure 2a) for every design cell: genetic architecture (rows) × pleiotropy × selection regime (columns), n = %d replicates per cell. ", NREP / 12),
        "equal-S: both environmental axes under equal selection; unequal-S: selection on one axis 8-fold weaker. ",
        DETECTION_DEF, median_only(DET, "panel")),
    D6_detection_demography = paste0(
        sprintf("Detection (as in Figure 2a) for the five demographic scenarios of the simulation deposit, named by their deme-size (N) and migration (m) patterns (n = %d replicates each). ", NREP / 5),
        DETECTION_DEF, median_only(DET, "block")),
    C2_offset_distribution = paste0(
        sprintf("Genetic-offset accuracy per replicate for each candidate set and for the causal loci, by genetic architecture (n = %d per level). ", NREP / 3),
        ACC3_DEF),
    R2_rank_by_architecture = paste0(
        sprintf("Rank distribution of Figure 2b split by genetic architecture (rows; n = %d replicates per level) and offset engine (columns). ", NREP / 3),
        RANK_DEF),
    R3_least_accurate_strata = paste0(
        "Percentage of replicates in which each candidate set (LFMM, EMMAX, RDA without structure correction and the three agreement sets) is the least accurate of the six (median of the three offset engines; ranks as in Figure 2b, tied sets split evenly), ",
        sprintf("(a) by design cell — genetic architecture × pleiotropy (column groups) × selection regime (columns), n = %d replicates per cell — and (b) by demographic scenario, n = %d each. ",
                NREP / 12, NREP / 5),
        sprintf("Colour: white at chance (%.1f%%), teal below, red above. ", CHANCE),
        sprintf("The 2/3 set is least accurate in at most %.1f%% of the replicates of any stratum (%s: %s, %d of %d replicates). ",
                MX23$least_pct, MX23$kind, MX23$stratum, as.integer(round(MX23$least_pct * MX23$n / 100)), MX23$n),
        "Kruskal–Wallis tests of each set's within-replicate rank (median of the three engines; average rank for tied sets), Holm-adjusted over the 12 tests, ",
        kw_text("design cell"), "; ", kw_text("demographic scenario"), ". ",
        "Candidate sets ordered as in Figure 2b."))

# Display items in manuscript order; A1 and R1 are parts of Figure 2 (own files, no own caption).
DISPLAY <- c("Fig2_simulation", "S_rda_correction", "B4_aucpr_methods", "A5_cutoff_sensitivity",
             "D5_detection_cells", "D6_detection_demography", "C2_offset_distribution",
             "R2_rank_by_architecture", "R3_least_accurate_strata")
PARTS <- c(A1_detection_plane = "Fig2_simulation", R1_rank_distribution = "Fig2_simulation")
role_of <- function(s) FIGS[[s]]$role
stopifnot(setequal(names(CAP), DISPLAY), setequal(names(FIGS), c(DISPLAY, names(PARTS))),
          all(vapply(DISPLAY, role_of, "") %in% c("main", "supp")), all(vapply(names(PARTS), role_of, "") == "part"))
TAG <- c(Fig2_simulation = "Figure 2", S_rda_correction = "Figure S1", B4_aucpr_methods = "Figure S2",
         A5_cutoff_sensitivity = "Figure S3", D5_detection_cells = "Figure S4",
         D6_detection_demography = "Figure S5", C2_offset_distribution = "Figure S6",
         R2_rank_by_architecture = "Figure S7", R3_least_accurate_strata = "Figure S8")
TITLE <- c(Fig2_simulation = "Detection and offset rank of the six candidate sets",
           S_rda_correction = "Why RDA without structure correction: ranking of the causal loci with vs without correction",
           B4_aucpr_methods = "Ranking of the causal loci by each single method (AUC-PR)",
           A5_cutoff_sensitivity = "Detection at five cutoffs",
           D5_detection_cells = "Detection by design cell",
           D6_detection_demography = "Detection by demographic scenario",
           C2_offset_distribution = "Offset accuracy, candidate sets and causal loci",
           R2_rank_by_architecture = "Offset rank by genetic architecture",
           R3_least_accurate_strata = "Least accurate set by design cell and demographic scenario")
stopifnot(setequal(names(TAG), DISPLAY), setequal(names(TITLE), DISPLAY))
file_line <- function(s) {
    parts <- names(PARTS)[PARTS == s]
    paste0(sprintf("`%s.svg` / `.png` · %.0f × %.0f mm", s, FIGS[[s]]$w, FIGS[[s]]$h),
           if (length(parts)) paste0(" · parts a and b also as ",
                                     paste(sprintf("`%s` (%.0f × %.0f mm)", parts, vapply(FIGS[parts], `[[`, 0, "w"),
                                                   vapply(FIGS[parts], `[[`, 0, "h")), collapse = " and ")) else "",
           sprintf(" · source: `%s`", FIGS[[s]]$source))
}

# =============================================================================
# Table S1: caption + Table_S1.md (the tables for the manuscript; Table_S1.tsv holds every value)
# =============================================================================
f2  <- function(x) formatC(x, format = "f", digits = 2)
fd3 <- function(x) sub("^-", "−", formatC(x, format = "f", digits = 3))      # signed, typographic minus
p_txt <- function(p) ifelse(p < 1e-4, "<0.0001", formatC(p, format = "g", digits = 2))
md_table <- function(D) {
    D <- as.data.table(D)[, lapply(.SD, function(v) trimws(as.character(v)))]
    c(paste0("| ", paste(names(D), collapse = " | "), " |"),
      paste0("|", paste(rep("---", ncol(D)), collapse = "|"), "|"),
      apply(D, 1, function(v) paste0("| ", paste(v, collapse = " | "), " |")))
}
TS1_TITLE <- "Within-replicate rank of the candidate sets, their distance to the causal loci, and dependence on the design"
TS1_CAP <- list(
    a = paste0(sprintf("(a) Percentage of the %d replicates in which a candidate set is the least accurate, among the two least accurate, the most accurate and among the two most accurate of the six, ", NREP),
               "under the median of the three offset engines (Gradient Forest, the LFMM2 geometric offset and RDA without structure correction) and under each engine. ",
               "Candidate sets as in Figure 2: LFMM, EMMAX, RDA without structure correction and the three agreement sets. ",
               "Tied sets split the tied ranks evenly; an empty set is ranked last. ",
               sprintf("95%% CI: percentile bootstrap over replicates (%s resamples). P: two-sided permutation test against chance (16.7%% for least and most accurate, 33.3%% for the two least and two most accurate), candidate-set labels shuffled within each replicate (%s permutations); P = (k + 1) / (%s + 1), k = permutations at least as far from chance as observed. Seed %d for both.",
                       format(S1_B, big.mark = ","), format(S1_B, big.mark = ","), format(S1_B, big.mark = ","), S1_SEED)),
    b = paste0("(b) Δ accuracy = accuracy of the candidate set − accuracy of the causal loci of the same replicate (−Kendall τ, median over the 100 landscape gardens and the three engines): ",
               "median, minimum and maximum over replicates, and the percentage of replicates with Δ ≥ 0, i.e. in which the set matched or exceeded the causal loci. ",
               sprintf("n = replicates with an offset (%s). ", if (nrow(EMPTY)) paste(sprintf("the %s is empty in %d", set_lab(EMPTY$label), EMPTY$N), collapse = "; ") else "all"),
               "Δ ≥ 0 is column `reach_pct` of `C_delta_stats.tsv`; its `win_pct` counts Δ > 0, the same numbers here because no replicate has Δ = 0 exactly."),
    c = "(c) Range of the 2/3 set's percentages over the 3 genetic architectures, the 5 demographic scenarios and the 12 design cells (median of the three engines; values in Figure S8 and `numbers.tsv`).",
    d = sprintf("(d) Sums of squares of accuracy (median of the three engines) by replicate and by candidate set, over the %d replicates in which all six sets have an offset; spread = range of the six set medians, and the median over replicates of the range of the six sets within a replicate (−Kendall τ).",
                RVP$n_rep[1]),
    e = "(e) Kruskal–Wallis tests of Figure S8: each set's within-replicate rank (median of the three engines; average rank for tied sets) across the 12 design cells and across the 5 demographic scenarios, Holm-adjusted over the 12 tests.")
TA_W <- dcast(TA[, .(engine = factor(engine, levels = ENGS), panel = factor(panel, levels = PANEL_ORDER), metric,
                     cell = sprintf("%s [%s–%s]", f2(value), f2(lo), f2(hi)), P = p_txt(p_perm))],
              engine + panel ~ metric, value.var = c("cell", "P"))
setorder(TA_W, engine, panel)
TAB <- list(
    a = TA_W[, .(Engine = ENG_LAB[as.character(engine)], `Candidate set` = set_lab(panel),
                 `Least accurate, % [95% CI]` = cell_bottom1, `P (least)` = P_bottom1,
                 `Two least accurate, % [95% CI]` = cell_bottom2, `P (two least)` = P_bottom2,
                 `Most accurate, % [95% CI]` = cell_top1, `P (most)` = P_top1,
                 `Two most accurate, % [95% CI]` = cell_top2, `P (two most)` = P_top2)],
    b = rbind(DB[match(PANEL_ORDER, panel), .(`Candidate set` = set_lab(panel), n = as.character(n),
                                              `Median Δ accuracy` = fd3(median), `Minimum Δ` = fd3(min),
                                              `Maximum Δ` = fd3(max), `Δ ≥ 0, % of replicates` = formatC(ge0_pct, format = "f", digits = 1))],
              data.table(`Candidate set` = "Spread over the six sets", n = "",
                         `Median Δ accuracy` = sprintf("%s to %s (magnitude %s–%s)", fd3(-SPREAD$abs_hi), fd3(-SPREAD$abs_lo),
                                                       fd3(SPREAD$abs_lo), fd3(SPREAD$abs_hi)),
                         `Minimum Δ` = "", `Maximum Δ` = "",
                         `Δ ≥ 0, % of replicates` = sprintf("%.1f–%.1f", SPREAD$ge0_lo, SPREAD$ge0_hi))),
    c = TC[, .(Stratification = kind, `Strata (replicates each)` = sprintf("%d (%d)", n_strata, n_per_stratum),
               `2/3 set least accurate, % (min–max)` = sprintf("%s–%s", f2(least_lo), f2(least_hi)),
               `2/3 set among the two least accurate, % (min–max)` = sprintf("%s–%s", f2(last2_lo), f2(last2_hi)))],
    d = data.table(Quantity = c(paste0("Variance of accuracy explained by ", VP_LAB[RVP$source], ", %"),
                                "Spread of the six set medians (−Kendall τ)",
                                "Median spread of the six sets within a replicate (−Kendall τ)", "Replicates"),
                   Value = c(f2(RVP$share_pct), fd3(RVP$panel_median_range[1]), fd3(RVP$within_rep_range_median[1]),
                             as.character(RVP$n_rep[1]))),
    e = KWT[, .(`Candidate set` = set_lab(panel), Stratification = sprintf("%s (%d)", factor, groups), n = as.character(n),
               H = f2(H), df = as.character(df), P = vapply(p, fmt_p, ""), `P (Holm)` = vapply(p_holm, fmt_p, ""))])
TAB_HEAD <- c(a = "(a) Within-replicate rank of the candidate sets", b = "(b) Offset accuracy relative to the causal loci",
              c = "(c) The 2/3 set across strata", d = "(d) Variance of accuracy", e = "(e) Dependence of the rank on design cell and demographic scenario")
TS1_MD <- c(sprintf("# Table S1 — %s", TS1_TITLE), "",
            sprintf("SS-Clines, %d replicates. Values: `Table_S1.tsv` (long format) and `numbers.tsv`; sources in Data S2.", NREP), "",
            unlist(lapply(names(TAB), function(k) c(sprintf("## %s", TAB_HEAD[[k]]), "", md_table(TAB[[k]]), "", TS1_CAP[[k]], ""))))
writeLines(TS1_MD, file.path(OUT, "Table_S1.md"))

md <- unlist(lapply(DISPLAY, function(s) c(sprintf("## %s — %s", TAG[[s]], TITLE[[s]]), "", file_line(s), "", CAP[[s]], "")))
md <- c(md, sprintf("## Table S1 — %s", TS1_TITLE), "",
        "`Table_S1.md` (tables) / `Table_S1.tsv` (every value, long format) · source: `rank_long.tsv`, `rank_summary.tsv`, `offset_delta_per_seed.tsv`, `C_delta_stats.tsv`, `reliability_variance.tsv`, `detection_per_seed.tsv`",
        "", paste(unlist(TS1_CAP), collapse = " "), "")
CAPTIONS <- c("# Simulation display items — draft captions", "",
              sprintf("SS-Clines corpus, %d replicates. Figure 2, Figures S1–S8 and Table S1; their source tables are Data S2 (`Data_S2.zip`, manifest `Data_S2_manifest.tsv`). Arial 7 pt, drawn at print size. Numbers: `numbers.tsv`; for Figure S1 also `figures_ssclines_rda/numbers.tsv` and `RDA_CORRECTION_BRIEF.md`.", NREP),
              "", md)
writeLines(CAPTIONS, file.path(OUT, "captions.md"))

NUMS <- rbindlist(NUM, use.names = TRUE, fill = TRUE)
setcolorder(NUMS, intersect(c("figure", "rule", "stratum", "engine", "metric", "value", "lo", "hi",
                              "expected", "n", "p_perm", "source"), names(NUMS)))
fwrite(NUMS, file.path(OUT, "numbers.tsv"), sep = "\t")
message(sprintf("  numbers.tsv %d rows", nrow(NUMS)))

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
    "each candidate set's accuracy on each half, centred it on the replicate's mean over its six candidate sets, and ",
    "correlated the two halves (Spearman–Brown-corrected to 100 gardens). We also measured agreement of the order of ",
    "the candidate sets across the three offset engines with Kendall's W, against a permutation null in which each ",
    "engine's values of the six sets were shuffled within the replicate (%s permutations). Replicates with an empty ",
    "set were excluded (%d of %d analysed)."),
    M3R$n_splits, format(KW$n_perm, big.mark = ","), M3R$n_rep, NREP)
TXT_RESULTS <- sprintf(paste0(
    "Mean accuracy barely differed among the candidate sets (medians within %s; the candidate set explained %s%% of the ",
    "variance, the replicate %s%%), but within a replicate the sets differed reproducibly: the split-half correlation ",
    "was %s–%s per engine (%s for the median of the three engines at 100 gardens), the least accurate set was ",
    "the same in both halves in %s–%s%% of replicates (chance %s%%), and the three offset engines agreed on the ",
    "order of the sets (Kendall's W = %s; %s under the null; P ≤ %s). Which set was best changed between replicates ",
    "(median range of the six sets within a replicate %s), which is why the pooled distributions overlap ",
    "while the ranks do not."),
    f_(RVP[1, panel_median_range], 3), f_(RVP[source == "panel", share_pct], 2),
    f_(RVP[source == "replicate", share_pct], 1),
    f_(min(RSH[engine != "median of 3", r_half]), 2), f_(max(RSH[engine != "median of 3", r_half]), 2),
    f_(M3R$r_full, 2), f_(min(RSH$worst_same_pct), 0), f_(max(RSH$worst_same_pct), 0), f_(100 / N_PANEL, 0),
    f_(KW$value, 2), f_(KW$null_mean, 2), formatC(KW$p_perm, format = "g", digits = 1),
    f_(RVP[1, within_rep_range_median], 3))
TXT_CAVEAT <- paste0(
    "Not covered by these tests: whether the 2/3 set ranks mid-field because it is built from the other ",
    "sets' SNPs. Present it as a hedge against the worst case, not as a way to select better loci.")

# Reproduce: every script between the offset13 candidate sets and this page. Paths relative to the
# repo root /mnt/data/eugene/ADAPTOGENE; journal = work/journal/18_ssclines_onefit_rule.html.
REPRO <- data.table(
    stage = c(rep("this report (driver: bash work/p5_figures.sh all)", 5), rep("upstream (commands in the journal)", 6)),
    script = c("benchmarks/mvp_j16_gallery.R", "benchmarks/mvp_panel_rank.R", "benchmarks/mvp_rank_ties.R",
               "benchmarks/mvp_rank_reliability.R", "benchmarks/mvp_ms_sim_figures.R",
               "benchmarks/mvp_build_snp_sets.R",
               "work/p4b_configs.sh + work/p4b_sweep.sh, then benchmarks/eval_offset_lind.R + benchmarks/mvp_panel_tables.R",
               "benchmarks/mvp_remeasure_600.R (SHARE_MULT=RDA=2)",
               "benchmarks/mvp_oracle_stats.R (+ MVP_SIZE_PANELS=1 on block b1)",
               "benchmarks/mvp_offset_engines.R",
               "benchmarks/mvp_ms_rda_correction.R"),
    output = c("benchmarks/mvp_eval/figures_ssclines_gallery13/ — per-replicate detection and offset tables, AUC-PR",
               "benchmarks/mvp_eval/figures_ssclines_offset13/rank/ — rank_long, rank_summary",
               "the tie convention of the within-replicate rank (tied sets split the tied ranks evenly), sourced by mvp_panel_rank.R and mvp_ms_sim_figures.R",
               "benchmarks/mvp_eval/figures_ssclines_offset13/rank/reliability_*.tsv",
               "benchmarks/mvp_eval/figures_ssclines_ms/ — figures, captions.md, Table S1, Data S2, numbers.tsv, README.md, this page",
               "benchmarks/mvp_eval/offset13/snp_sets_summary.tsv + the six candidate sets per replicate",
               "benchmarks/mvp_eval/offset13/ — garden_performance.tsv, phase1_seed_medians_solo.tsv, panel_pr_recomputed.tsv",
               "benchmarks/mvp_eval/remeasure600/rdaunc_rda2x/ — detection at every cutoff and window (Figure S3)",
               "benchmarks/mvp_eval/figures_ssclines_offset13/ (+ size_b1/) — delta vs causal loci, size curve",
               "benchmarks/mvp_eval/figures_ssclines_offset13/ — backward gardens, RDA engine corrected vs uncorrected",
               "benchmarks/mvp_eval/figures_ssclines_rda/ — S_rda_correction (Figure S1), numbers.tsv, RDA_CORRECTION_BRIEF.md"),
    journal = c("Step 19", "Steps 18, 22", "Step 22", "Step 21", "Steps 19–21, 23–24", "Step 13", "Steps 14–15",
                "Steps 14, 16", "Step 17", "Step 17", "Steps 12, 16, 23"))
REPRO_NOTE <- c(
    "Repo root on lab (host fidel): /mnt/data/eugene/ADAPTOGENE. Image cline-go:latest (93bd6025); docker via `nix shell nixpkgs#docker-client -c docker`.",
    "Corpus: MVP_ARM=primary_ssclines, MVP_ADDED = the five SS-Clines block tags, MVP_N_EXPECT=600 (benchmarks/mvp_arm.R). Every arm input is passed explicitly — the scripts' defaults name older arms.",
    sprintf("Journal: work/journal/18_ssclines_onefit_rule.html (Steps 12–24). Logs: work/p5_logs/ (the 2026-10-03 logs in work/p5_logs/_2026-10-03/). Seeds: jitter 16, bootstrap 16 / 20261003, splits and permutations 20261003; Table S1 bootstrap and permutations %d (%s each).",
            S1_SEED, format(S1_B, big.mark = ",")),
    "Reproduction check 2026-10-03: re-running the rank step reproduced all 8 rank / reliability tables byte for byte.",
    "Tie convention 2026-10-14: tied sets split the tied ranks evenly (benchmarks/mvp_rank_ties.R). The rank rebuild left rank_long, rank_gaps, rank_gap_summary and winrate byte-identical and changed 137 of the 1,512 rows of rank_summary, all of them through the 25 tied units.",
    "Data keys: numbers.tsv and the source tables in Data S2 call the agreement sets `1/3 methods`, `2/3 methods`, `3/3 methods`; figures, captions and Table S1 call them the 1/3, 2/3 and 3/3 sets.",
    "Superseded sets (read-only): figures_ssclines_ms/_2026-10-03/ (before the 2026-10-14 numbering, vocabulary and tie convention) and figures_ssclines_ms/_offset12/ (journal-16 candidate sets). Dropped 2026-10-14: D5_offset_cells, kept in figures_ssclines_ms/_removed/.")
README <- c("# Simulation results — report README",
            "", sprintf("SS-Clines, %d replicates. Contains every display item (Figure 2, Figures S1–S8, Table S1, Data S2), caption and number needed to write the simulation part of the manuscript's Results (and the matching Methods). Open `index.html` for the figures; numbers in `numbers.tsv`; captions in `captions.md`; Table S1 in `Table_S1.md`; source tables in `Data_S2.zip`.", NREP),
            "", "## Text-only result: the within-replicate ranking is signal, not noise", "",
            "**Methods (draft).** ", TXT_METHODS, "", "**Results (draft).** ", TXT_RESULTS, "", paste0("*", TXT_CAVEAT, "*"),
            "", "## Reproduce", "", paste0("- ", REPRO_NOTE), "",
            "| stage | script | output | journal 18 |", "|---|---|---|---|",
            REPRO[, sprintf("| %s | `%s` | %s | %s |", stage, script, output, journal)])
writeLines(README, file.path(OUT, "README.md"))

# =============================================================================
# Download-all zip + index.html
# =============================================================================
FILES <- c(DISPLAY, names(PARTS))
files <- c(paste0(FILES, ".svg"), paste0(FILES, ".png"), "captions.md", "numbers.tsv", "README.md",
           "Table_S1.md", "Table_S1.tsv")
stopifnot(all(file.exists(file.path(OUT, files))))
ZIP <- "sim_figures_ms.zip"
old <- setwd(OUT)
unlink(ZIP)
stopifnot(utils::zip(ZIP, files, flags = "-9Xq") == 0L)
setwd(old)

esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
                     gsub(">", "&gt;", x, fixed = TRUE) }
# `code` spans of the markdown text -> <code> on the page
code_html <- function(x) gsub("`([^`]*)`", "<code>\\1</code>", esc(x))
kb <- function(f) { b <- file.size(file.path(OUT, f)); if (b > 1e6) sprintf("%.1f MB", b / 1e6) else sprintf("%.0f kB", b / 1e3) }
part_links <- function(s) {
    parts <- names(PARTS)[PARTS == s]
    if (!length(parts)) return("")
    paste(vapply(seq_along(parts), function(i) sprintf(
        '<a class="btn" href="%s.svg" download>part %s SVG · %s</a><a class="btn" href="%s.png" download>part %s PNG · %s</a>',
        parts[i], letters[i], kb(paste0(parts[i], ".svg")), parts[i], letters[i], kb(paste0(parts[i], ".png"))), ""),
        collapse = "")
}
card <- function(s, cls = "") sprintf(paste0(
    '<figure class="card %s" id="%s">',
    '<figcaption class="head"><span class="tag">%s</span><span class="ttl">%s</span>',
    '<span class="dim">%.0f × %.0f mm</span></figcaption>',
    '<div class="img"><img src="%s.svg" alt="%s" style="--w:%.0f" loading="lazy"></div>',
    '<p class="cap">%s</p>',
    '<div class="dl"><a class="btn" href="%s.svg" download>SVG · %s</a>',
    '<a class="btn" href="%s.png" download>PNG 600 dpi · %s</a>%s</div></figure>'),
    cls, s, esc(TAG[[s]]), esc(TITLE[[s]]), FIGS[[s]]$w, FIGS[[s]]$h, s, esc(TITLE[[s]]),
    FIGS[[s]]$w, esc(CAP[[s]]), s, kb(paste0(s, ".svg")), s, kb(paste0(s, ".png")), part_links(s))
html_table <- function(D, cls = "repro") {
    D <- as.data.table(D)
    paste0('<div style="overflow-x:auto"><table class="', cls, '"><thead><tr>',
           paste(sprintf("<th>%s</th>", esc(names(D))), collapse = ""), "</tr></thead><tbody>",
           paste(apply(D, 1, function(v) paste0("<tr>", paste(sprintf("<td>%s</td>", esc(trimws(as.character(v)))), collapse = ""),
                                                "</tr>")), collapse = ""),
           "</tbody></table></div>")
}
TS1_HTML <- unlist(lapply(names(TAB), function(k) c(sprintf('<h3>%s</h3>', esc(TAB_HEAD[[k]])), html_table(TAB[[k]]),
                                                   sprintf('<p class="cap">%s</p>', code_html(TS1_CAP[[k]])))))
HEADER_NOTE <- sprintf("Editable SVG (Arial 7 pt, print size) and PNG 600 dpi; captions are drafts. Candidate sets: the cutoff is 0.25%% per method and predictor (RDA 0.5%%), RDA without structure correction, 5 kb agreement window. Built %s.",
                       format(Sys.time(), "%Y-%m-%d %H:%M"))
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
'h3{font-size:14px;margin:14px 0 6px}',
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
sprintf('<header><div><h1>Simulation figures — SS-Clines, %d replicates</h1><p>%s</p></div>', NREP, esc(HEADER_NOTE)),
sprintf(paste0('<div class="dl"><button class="btn" id="ps" type="button" aria-pressed="false">Show print size</button>',
               '<a class="btn" href="captions.md" download>captions.md</a><a class="btn" href="numbers.tsv" download>numbers.tsv</a>',
               '<a class="btn" href="Table_S1.md" download>Table_S1.md</a><a class="btn" href="README.md" download>README.md</a>',
               '<a class="btn" href="%s" download>Data S2 (.zip, %s)</a>',
               '<a class="btn primary" href="%s" download>Download all (.zip, %s)</a></div></header>'),
        DS2_ZIP, kb(DS2_ZIP), ZIP, kb(ZIP)),
'<main><h2>Figure 2</h2>',
'<p class="note">(a) Detection; (b) offset accuracy, rank of each candidate set within a replicate. 180 mm wide; parts a and b also as separate files.</p>',
'<section class="stack">', card("Fig2_simulation"), '</section>',
'<h2>Supplementary figures</h2>',
'<section class="stack">', vapply(setdiff(DISPLAY, "Fig2_simulation"), card, ""), '</section>',
sprintf('<h2>Table S1 — %s</h2>', esc(TS1_TITLE)),
'<section class="card">', TS1_HTML,
'<div class="dl"><a class="btn" href="Table_S1.md" download>Table_S1.md</a><a class="btn" href="Table_S1.tsv" download>Table_S1.tsv (every value)</a></div></section>',
'<h2>Data S2 — source tables</h2>',
'<section class="card"><p class="note" style="margin:0 0 8px">Every table a figure or Table S1 is drawn from, at its path under <code>benchmarks/mvp_eval/</code>, with <code>manifest.tsv</code>.</p>',
html_table(MAN[, .(file, `display items` = display_items, rows, columns, description, `written by` = written_by)]),
sprintf('<div class="dl"><a class="btn" href="%s" download>Data S2 (.zip, %s)</a><a class="btn" href="Data_S2_manifest.tsv" download>manifest</a></div></section>', DS2_ZIP, kb(DS2_ZIP)),
'<h2>Text-only result — the within-replicate ranking is signal, not noise</h2>',
'<section class="card"><p class="note" style="margin:0 0 8px">No figure (user decision 2026-10-03): these tests go into the manuscript text. Draft wording; numbers from <code>figures_ssclines_offset13/rank/reliability_*.tsv</code>.</p>',
sprintf('<p><strong>Methods.</strong> %s</p>', esc(TXT_METHODS)),
sprintf('<p><strong>Results.</strong> %s</p>', esc(TXT_RESULTS)),
sprintf('<p class="note"><em>%s</em></p>', esc(TXT_CAVEAT)), '</section>',
'<h2>Reproduce</h2><section class="card">',
paste0('<ul class="note" style="margin:0 0 10px;padding-left:18px">', paste(sprintf('<li>%s</li>', code_html(REPRO_NOTE)), collapse = ""), '</ul>'),
'<div style="overflow-x:auto"><table class="repro"><thead><tr><th>stage</th><th>script</th><th>output</th><th>journal 18</th></tr></thead><tbody>',
paste(REPRO[, sprintf('<tr><td>%s</td><td><code>%s</code></td><td>%s</td><td>%s</td></tr>', esc(stage), esc(script), esc(output), esc(journal))], collapse = ""),
'</tbody></table></div></section></main>',
# Every image is shown at one common scale (width in mm / 180 of the row), so relative figure
# sizes stay true; "print size" shows each at its physical width (CSS mm) to judge 7 pt text.
'<script>document.getElementById("ps").addEventListener("click",function(){var on=document.body.classList.toggle("print");this.setAttribute("aria-pressed",on);this.textContent=on?"Fit to page":"Show print size";});</script>',
'</body></html>')
writeLines(html, file.path(OUT, "index.html"))
# every file the page links to exists
links <- unique(unlist(regmatches(html, gregexpr('(?<=src=")[^"#]+|(?<=href=")[^"#]+', html, perl = TRUE))))
stopifnot(length(links) > 0L, all(file.exists(file.path(OUT, links))))

# =============================================================================
# Vocabulary guard over everything a reader sees (figure text was checked at each save_ms)
# =============================================================================
# the page as read: no style or script blocks, no code spans, no tags (attributes go with them)
page_text <- function(h) {
    h <- paste(h, collapse = "\n")
    h <- gsub("(?s)<(style|script)>.*?</\\1>", " ", h, perl = TRUE)
    h <- gsub("(?s)<code>.*?</code>", " ", h, perl = TRUE)
    gsub("<[^>]+>", " ", h, perl = TRUE)
}
n_txt <- c(captions.md = vocab_check(CAPTIONS, "captions.md"),
           Table_S1.md = vocab_check(TS1_MD, "Table_S1.md"),
           Table_S1.tsv = vocab_check(unlist(TS1[, .(set, engine, stratum, metric)]), "Table_S1.tsv"),
           Data_S2_manifest = vocab_check(unlist(MAN[, .(display_items, description)]), "Data_S2_manifest.tsv"),
           README.md = vocab_check(README, "README.md"),
           index.html = vocab_check(page_text(html), "index.html"),
           titles = vocab_check(c(TAG, TITLE, TS1_TITLE, unlist(TAB_HEAD)), "titles and tags"))
message("vocabulary: no banned wording in ", paste(sprintf("%s (%d strings)", names(n_txt), n_txt), collapse = ", "))
# Every caption that names RDA introduces it as "RDA without structure correction" (user, 2026-10-14).
RDA_CAPS <- c(unlist(CAP), table_s1 = paste(unlist(TS1_CAP), collapse = " "))
no_rda_name <- names(RDA_CAPS)[grepl("\\bRDA\\b", RDA_CAPS) & !grepl("RDA without structure correction", RDA_CAPS, fixed = TRUE)]
if (length(no_rda_name)) stop("caption(s) name RDA without 'RDA without structure correction': ", paste(no_rda_name, collapse = ", "))
message(sprintf("RDA named 'RDA without structure correction' in all %d captions that mention it", sum(grepl("\\bRDA\\b", RDA_CAPS))))

message(sprintf("OK: %d display figures + %d parts (SVG + PNG), Table S1, Data S2, captions.md, numbers.tsv, index.html, %s -> %s",
                length(DISPLAY), length(PARTS), ZIP, OUT))
