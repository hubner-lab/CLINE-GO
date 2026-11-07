# [landscape figs copy of work/ld_window_lab/scratch_figs/sim_figures_ld.R, made by work/mvp_land/scratch/make_figs_land.py]
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
SIZ  <- NA_character_   # [landscape figs] size section dropped
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

# ---- [LD-lab scratch copy] fixed-number checks -> warnings, logged to SOFT_LOG (one line each)
SOFT_LOG <- Sys.getenv("SOFT_LOG", "/tmp/soft_asserts.tsv")
soft_note <- function(where, msg) {
    warning("SOFT ASSERT [", where, "]: ", msg, call. = FALSE)
    cat(where, "\t", gsub("[\t\n]+", " ", msg), "\n", sep = "", file = SOFT_LOG, append = TRUE)
}
soft_stopifnot <- function(.where, ...) {
    pf <- parent.frame(); ex <- as.list(substitute(list(...)))[-1L]
    for (e in ex) {
        ok <- tryCatch(isTRUE(all(eval(e, pf))), error = function(err) FALSE)
        if (!ok) soft_note(.where, paste(deparse(e), collapse = " "))
    }
    invisible(TRUE)
}
soft_stop <- function(.where, ...) soft_note(.where, paste0(...))
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
message("[LD-lab] sets without an offset (< 3 SNPs): ", nrow(off_missing), "; empty: ", sum(EMPTY$N))
stopifnot(nrow(off_missing) == DET[n < 3, .N],
          nrow(merge(off_missing, DET[n < 3, .(seed, label = as.character(label))],
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
message("[LD-lab] Figure 2a medians, new (TPm/FPm) vs quoted (TP/FP):"); print(fx); print(EMPTY)
soft_stopifnot("FIX_DET Fig2a medians + empty 3/3", nrow(fx) == N_PANEL, all(fx$TPm == fx$TP), all(fx$FPm == fx$FP),
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
if (!identical(ORD_OLD$panel, PANEL_ORDER)) soft_stop("ORD tie rule", "the tie rule changed the order of the candidate sets")
# [2026-10-18] the order of Figure 2b, as in the 2026-10-14 build (display work must not move it)
soft_stopifnot("PANEL_ORDER Fig2b", identical(PANEL_ORDER, c("2/3 methods", "1/3 methods", "LFMM", "3/3 methods", "RDA", "EMMAX")))
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
    soft_stopifnot("A1 y range", pp$y$continuous_range[1] > 0, y_lo < A1_LOW)
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
# Within-replicate rank (R1 main, R2 supp) -- the V4 view of mvp_panel_rank.R
# =============================================================================
message("=== rank")
# Empty panels in the rank table are exactly the empty detection panels, one rank per engine.
stopifnot(setequal(RL[empty == TRUE, unique(paste(seed, panel))],
                   DET[n < 3, paste(seed, label)]))   # [LD-lab] < 3 SNPs: no offset, ranked last
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

ARCH_WRAP <- c("oligogenic" = "oligogenic", "moderately polygenic" = "moderately\npolygenic",
               "highly polygenic" = "highly\npolygenic")
# =============================================================================
# Figure S2: AUC-PR of the three single methods (gallery code B4), drawn as S1 (user, 2026-10-14)
# =============================================================================
message("=== S2 AUC-PR per method")
RNUM <- setnames(data.table(k_ = character(), value = character()), "k_", "key")  # [landscape figs] no S1 brief: tie logged as soft assert
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
        soft_stop("S2 AUC-PR", "S2: pooled AUC-PR of ", meth, " = ", v, " differs from figures_ssclines_rda/numbers.tsv or ", AU_FIX[[meth]])
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
# [landscape figs] numbers + figure list (the full report's captions / Data S2 / zip are not built)
# =============================================================================
NUMS <- rbindlist(NUM, use.names = TRUE, fill = TRUE)
setcolorder(NUMS, intersect(c("figure", "rule", "stratum", "engine", "metric", "value", "lo", "hi",
                              "expected", "n", "p_perm", "source"), names(NUMS)))
fwrite(NUMS, file.path(OUT, "numbers.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(FIGS), function(k) data.table(stem = k, role = FIGS[[k]]$role,
       w_mm = FIGS[[k]]$w, h_mm = FIGS[[k]]$h, source = FIGS[[k]]$source))), file.path(OUT, "figures.tsv"), sep = "\t")
message(sprintf("=== landscape figs done: %d figures, numbers.tsv %d rows -> %s", length(FIGS), nrow(NUMS), OUT))
