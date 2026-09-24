#!/usr/bin/env Rscript
# =============================================================================
# mvp_j16_gallery.R -- labelled presentation variants for the pooled SS-Clines corpus
# (600 replicates), for the user to choose from before journal 16 is written.
#
# FAMILIES (codes are the contract with work/journal/16g_simulation_variants.Rmd):
#   A  detection, representation   A1 TP-vs-FP plane, contours   A2 same, raw points
#                                  A3 paired arrows LFMM -> 2/3  A4 composition (means)
#   B  detection, one scalar       B1 Pareto dominance           B2 precision + yield
#                                  B3 F1 on causal loci          B4 AUC-PR (single methods)
#                                  B5 net discoveries TP - FP
#   C  offset, representation      C1 estimation (dots + median + bootstrap CI)
#      (delta vs causal loci)      C2 raincloud   C3 rank share   C4 ECDF   C5 reach tile
#   D  split                       D1 pooled  D2 genic level  D3 x regime  D4 x pleiotropy
#                                  D5 all 12 cells  D6 demography block
#                                  (each: an offset view Dk_offset + a detection view Dk_detection)
#
# DEFINITIONS, stated once:
#   detection, per seed x rule, fixed operating points (LFMM top 100, RDA top 100,
#     EMMAX p < 1e-4, 5 kb agreement window) -- offset12_ssclines_pooled/panel_pr_recomputed.tsv
#     TP = causal + linked_neutral hits, FP = background_neutral hits (LG 11-20) = n - TP.
#     Empty panels are kept (n = 0): 3/3 on 67 seeds, EMMAX on 8.
#   offset, per seed x panel: accuracy = -tau, median over the 100 landscape gardens
#     (upstream), then median over the 3 working engines (GFoffset, LFMM2offset,
#     RDA-uncorrected); delta = panel - causal loci, paired within seed. Exactly the
#     mvp_oracle_stats.R aggregation; regression-checked against journal 15's
#     delta_per_seed.tsv below. 3/3 has 514 seeds, EMMAX 536 (panels < 3 SNPs dropped).
#
# EXCLUSIONS in the A-plane: causal loci (no linked hits by construction -- TP 7/43/460
#   at FP 0, so it would sit BELOW every GEA panel at oligo/moderate and read as "2/3
#   beats the truth set") and the random panel (drawn from background only: TP 0).
#   Both appear in A4 as labelled reference bars.
#
#   FIG_OUT     default benchmarks/mvp_eval/figures_ssclines_j16_gallery
#   OFFSET_DIR  default offset12_ssclines_pooled
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT via benchmarks/mvp_arm.R (mandatory, 600)
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggdist); library(ggrepel)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset12_ssclines_pooled")
J15  <- Sys.getenv("J15_DIR", file.path(EVAL, "figures_ssclines_pooled_j15"))
DET  <- Sys.getenv("DET_DIR", file.path(EVAL, "detection600"))
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_j16_gallery"))
BOOT_SEED <- 16L
N_BOOT    <- 2000L
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ vocabulary
ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus ships the typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
BLOCK_LAB <- c(ssclines_nvar_mvar      = "N variable, m variable",
               ssclines_ncline_ns      = "N cline north-south",
               ssclines_nequal_mconst  = "N equal, m constant",
               ssclines_ncline_ctredge = "N cline centre-edge",
               ssclines_nequal_mbreaks = "N equal, m breaks")
RULES  <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
ORACLE <- "causal loci"
RANDOM <- "random, size-matched"
OFF_ORDER <- c(RULES, RANDOM)

# Rule colours as mvp_dist_panel.R: sage shades = combinations, amber = single methods.
MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49", sage = "#66a182",
           navy = "#2e4057", grey = "#8d96a3")
sh <- function(b, n) c(b, grDevices::colorRampPalette(c(b, "white"))(6)[2:(1 + (n - 1) %/% 2 + (n - 1) %% 2)],
                          grDevices::colorRampPalette(c(b, "black"))(6)[2:(1 + (n - 1) %/% 2)])[seq_len(n)]
g3 <- sh(MINOU[["sage"]], 3); a3 <- sh(MINOU[["amber"]], 3)
RULE_COL <- c("2/3 methods" = g3[1], "3/3 methods" = g3[2], "1/3 methods" = g3[3],
              "RDA" = a3[1], "EMMAX" = a3[2], "LFMM" = a3[3],
              "causal loci" = MINOU[["navy"]], "random, size-matched" = MINOU[["grey"]])
ZERO <- geom_vline(xintercept = 0, linetype = "dashed", colour = CLINEGO_THRESHOLD, linewidth = 0.4)

# ------------------------------------------------------------------ helpers
MANIFEST <- list()
save_fig <- function(code, stem, p, w, h, what) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h)
    MANIFEST[[length(MANIFEST) + 1L]] <<- data.table(code = code, stem = stem, what = what)
    message(sprintf("  %-12s %s", code, stem))
}
emit <- function(dt, stem) {
    if ("seed" %in% names(dt))
        stopifnot(all(dt$seed %in% PRIM$seed), uniqueN(dt$seed) <= mvp_n_expect())
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
}
# The RNG is seeded ONCE (below, before the first summary), not inside boot_ci(): a per-call
# set.seed() would hand every group the same resampling stream, so CIs of adjacent rows would
# move together -- the wrong property for a figure whose question is "are these rules
# distinguishable?". Seeded once, the stream is still fully reproducible because the script's
# evaluation order is fixed.
boot_ci <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3L) return(c(NA_real_, NA_real_))
    m <- apply(matrix(sample(x, length(x) * N_BOOT, replace = TRUE), nrow = N_BOOT), 1, median)
    unname(quantile(m, c(0.025, 0.975)))
}
p_two <- function(d) {
    d <- d[is.finite(d)]
    if (length(d) < 6L || all(d == 0)) return(NA_real_)
    tryCatch(stats::wilcox.test(d, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
}
summ_delta <- function(D, by) D[, {
    ci <- boot_ci(delta)
    .(n = .N, median = median(delta), ci_lo = ci[1], ci_hi = ci[2],
      q25 = quantile(delta, 0.25), q75 = quantile(delta, 0.75),
      win_pct = 100 * mean(delta > 0), reach_pct = 100 * mean(delta >= 0),
      p_wilcox_two_sided = p_two(delta))
}, by = by]
facet_of <- function(rows = NULL, cols = NULL, scales = "fixed") {
    if (is.null(rows) && is.null(cols)) return(NULL)
    facet_grid(rows = if (length(rows)) vars(!!!rlang::syms(rows)) else NULL,
               cols = if (length(cols)) vars(!!!rlang::syms(cols)) else NULL,
               scales = scales)
}
ARCH_FACET <- facet_of(cols = "arch")

# ------------------------------------------------------------------ manifest
man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
PRIM[, arch   := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]
# "unequal-S" first: grepl("equal-S") matches both (mvp_pleiotropy_report.R:94).
PRIM[, regime := factor(fifelse(grepl("unequal-S", architecture), "unequal-S", "equal-S"),
                        levels = c("equal-S", "unequal-S"))]
PRIM[, pleio  := factor(fifelse(ispleiotropy == 1L, "pleiotropy", "no pleiotropy"),
                        levels = c("no pleiotropy", "pleiotropy"))]
PRIM[, block  := factor(BLOCK_LAB[added], levels = BLOCK_LAB)]
stopifnot(!anyNA(PRIM$arch), !anyNA(PRIM$block),
          nrow(PRIM[, .N, by = .(arch, pleio, regime)]) == 12L,
          PRIM[, .N, by = .(arch, pleio, regime)][, all(N == mvp_n_expect() / 12L)])
COV <- PRIM[, .(seed, arch, regime, pleio, block)]
set.seed(BOOT_SEED)

# =============================================================================
# DATA 1 -- offset delta vs causal loci, per seed
# =============================================================================
WORKING <- c("GFoffset", "LFMM2offset", "RDA-uncorrected")
SET_LAB <- c(gea_best = "2/3 methods", gea_union = "1/3 methods", gea_strict = "3/3 methods",
             gea_lfmm_only = "LFMM", gea_rda_only = "RDA", gea_emmax_only = "EMMAX",
             adaptive = ORACLE, random_matched = RANDOM)
acc <- fread(file.path(EVAL, OFF, "phase1_seed_medians_solo.tsv"), colClasses = c(seed = "character"))
acc <- mvp_require_rows(acc[seed %in% PRIM$seed & method_label %in% WORKING &
                            marker_set %in% names(SET_LAB)], "offset accuracy rows")
acc[, accuracy := -tau]                       # never abs()
SA <- acc[, .(accuracy = median(accuracy), n_methods = .N), by = .(seed, marker_set)]
SA[, label := SET_LAB[marker_set]]
W  <- dcast(SA, seed ~ label, value.var = "accuracy")
stopifnot(ORACLE %in% names(W), !anyNA(W[[ORACLE]]), nrow(W) == nrow(PRIM))
OFFD <- rbindlist(lapply(OFF_ORDER, function(p)
    data.table(seed = W$seed, label = p, accuracy = W[[p]], oracle = W[[ORACLE]])))[is.finite(accuracy)]
OFFD[, delta := accuracy - oracle]
OFFD <- merge(OFFD, COV, by = "seed")
OFFD[, label := factor(label, levels = OFF_ORDER)]

# Regression tie to journal 15: same seeds, same labels, same deltas.
j15 <- fread(file.path(J15, "delta_per_seed.tsv"), colClasses = c(seed = "character"))
chk <- merge(OFFD[, .(seed, label = as.character(label), delta)],
             j15[, .(seed, label, delta_j15 = delta)], by = c("seed", "label"), all = TRUE)
stopifnot(!anyNA(chk$delta), !anyNA(chk$delta_j15),
          max(abs(chk$delta - chk$delta_j15)) < 1e-9)
message(sprintf("regression vs journal 15: %d seed x panel deltas, max |diff| = %.1e",
                nrow(chk), max(abs(chk$delta - chk$delta_j15))))
emit(OFFD[, .(seed, label, arch, regime, pleio, block, accuracy, oracle, delta)],
     "offset_delta_per_seed")

# =============================================================================
# DATA 2 -- detection composition, per seed
# =============================================================================
DSET <- c(best = "2/3 methods", union = "1/3 methods", intersect3 = "3/3 methods",
          solo_lfmm = "LFMM", solo_rda = "RDA", solo_emmax = "EMMAX",
          truth = ORACLE, rand_best1 = RANDOM)
pr <- fread(file.path(EVAL, OFF, "panel_pr_recomputed.tsv"), colClasses = c(seed = "character"))
pr <- mvp_require_rows(pr[seed %in% PRIM$seed & set %in% names(DSET)], "panel composition rows")
pr[, label := factor(DSET[set], levels = c(RULES, ORACLE, RANDOM))]
pr[, `:=`(TP = n_causal + n_linked, FP = n - n_causal - n_linked)]
pr[, `:=`(prec_cl = fifelse(n > 0, TP / n, NA_real_),
          prec_c  = fifelse(n > 0, n_causal / n, NA_real_),
          rec_c   = n_causal / causal_total)]
pr[, F1 := fifelse(n_causal > 0, 2 * prec_c * rec_c / (prec_c + rec_c), 0)]
pr[, net := TP - FP]
pr <- merge(pr, COV, by = "seed")
stopifnot(all(pr$FP >= 0), pr[, .N, by = label][, all(N == nrow(PRIM))],
          !anyDuplicated(pr[, .(seed, label)]))
DETR <- pr[label %in% RULES]
DETR[, label := factor(as.character(label), levels = RULES)]
emit(pr[, .(seed, label, arch, regime, pleio, block, n, n_causal, n_linked, TP, FP,
            causal_total, prec_cl, prec_c, rec_c, F1, net)], "detection_per_seed")

det_stats <- function(D, by) D[, .(
    n_seeds = .N, n_empty = sum(n == 0),
    TP_median = median(TP), TP_q25 = quantile(TP, .25), TP_q75 = quantile(TP, .75),
    FP_median = median(FP), FP_q25 = quantile(FP, .25), FP_q75 = quantile(FP, .75),
    prec_cl_median = median(prec_cl, na.rm = TRUE), F1_median = median(F1),
    net_median = median(net)), by = by]
DS <- rbind(det_stats(pr, c("label", "arch")),
            det_stats(pr, "label")[, arch := "pooled"], use.names = TRUE)
emit(DS, "detection_stats")

message("=== A: detection representation")
# =============================================================================
# A -- the TP-vs-FP plane
# =============================================================================
MED <- function(D, by) D[, .(TP = median(TP), FP = median(FP)), by = by]
# Iso-precision guides: precision = TP / (TP + FP)  ->  TP = FP * p / (1 - p). Straight
# lines through the origin on a sqrt-sqrt plot too, so the sqrt axes keep them readable.
iso_lines <- function(D) {
    fmax <- max(D$FP); tmax <- max(D$TP)
    rbindlist(lapply(c(0.5, 0.75, 0.9), function(p) {
        k <- p / (1 - p); f <- seq(0, fmax, length.out = 60)
        data.table(p = p, FP = f, TP = f * k)[TP <= tmax]
    }))
}
# A 2-D density needs spread in both (sqrt) dimensions; 3/3 has FP = 0 on most seeds
# and EMMAX sits on a handful of integers. Those get the median marker only.
contour_ok <- function(D, by) D[, .(ok = .N >= 10L && IQR(sqrt(FP)) > 0 && IQR(sqrt(TP)) > 0),
                                by = by]

plane_plot <- function(D, rows = NULL, cols = "arch", title, show_contours = TRUE,
                       show_points = FALSE) {
    by <- c("label", rows, cols)
    M  <- MED(D, by)
    ok <- contour_ok(D, by)
    DC <- merge(D, ok[ok == TRUE], by = by)
    ISO <- iso_lines(D)
    n_facets <- prod(vapply(c(rows, cols), function(v) uniqueN(D[[v]]), 1L))
    p <- ggplot(D, aes(FP, TP)) +
        geom_line(data = ISO, aes(group = p), linetype = "dotted", colour = CLINEGO_COL$fg,
                  linewidth = 0.3, alpha = 0.6)
    # Label the guides only on small grids; on D5 (12 facets) / D6 (5) the same three labels
    # would repeat in every small panel. The caption names the three levels instead.
    if (n_facets <= 3L)
        p <- p + geom_text(data = ISO[, .SD[.N], by = p], aes(label = sprintf("precision %.2f", p)),
                           size = 2.3, hjust = 1.05, vjust = -0.4, colour = CLINEGO_COL$fg, alpha = 0.7)
    if (show_points)
        p <- p + geom_point(aes(colour = label), size = 0.45, alpha = 0.22)
    if (show_contours && nrow(DC))
        p <- p + geom_density_2d(data = DC, aes(colour = label), contour_var = "ndensity",
                                 breaks = c(0.25, 0.6), linewidth = 0.45)
    p + geom_point(data = M, aes(fill = label), shape = 21, size = 2.8, colour = "black",
                   stroke = 0.4) +
        geom_label_repel(data = M, aes(label = label, colour = label), size = 2.4,
                         label.size = 0.15, label.padding = 0.12, min.segment.length = 0,
                         seed = BOOT_SEED, show.legend = FALSE, max.overlaps = Inf) +
        facet_of(rows, cols) +
        scale_x_sqrt() + scale_y_sqrt() +
        scale_colour_manual(values = RULE_COL, guide = "none") +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        labs(x = "false positives per seed (background loci, sqrt scale)",
             y = "true hits per seed (causal + linked, sqrt scale)", title = title) +
        theme_clinego()
}
emit(MED(DETR, c("label", "arch")), "A_plane_medians")
save_fig("A1", "A1_plane_contours",
         plane_plot(DETR, title = "A1  detection plane, density contours + median"),
         13, 4.8, "TP-vs-FP plane: per-rule density contours (25%/60% of peak) + median, by genic level")
save_fig("A2", "A2_plane_points",
         plane_plot(DETR, title = "A2  detection plane, every seed + median",
                    show_contours = FALSE, show_points = TRUE),
         13, 4.8, "TP-vs-FP plane: every seed as a point + median, by genic level")

# A3 -- paired arrows, LFMM -> 2/3, one per seed
AR <- dcast(DETR[label %in% c("LFMM", "2/3 methods")], seed + arch ~ label, value.var = c("TP", "FP"))
setnames(AR, c("TP_2/3 methods", "FP_2/3 methods"), c("TP_best", "FP_best"))
AR[, `:=`(dTP = TP_best - TP_LFMM, dFP = FP_best - FP_LFMM)]
AR[, direction := fcase(dTP >= 0 & dFP <= 0 & (dTP > 0 | dFP < 0), "2/3 dominates",
                        dTP <= 0 & dFP >= 0 & (dTP < 0 | dFP > 0), "LFMM dominates",
                        dTP == 0 & dFP == 0, "identical",
                        default = "trade-off")]
AR[, direction := factor(direction, levels = c("2/3 dominates", "trade-off", "identical",
                                                "LFMM dominates"))]
emit(AR, "A3_arrows_per_seed")
emit(AR[, .N, by = .(arch, direction)][, pct := 100 * N / sum(N), by = arch][order(arch, direction)],
     "A3_arrow_direction_share")
DIR_COL <- c("2/3 dominates" = CLINEGO_RETAINED, "trade-off" = MINOU[["grey"]],
             "identical" = MINOU[["navy"]], "LFMM dominates" = CLINEGO_REMOVED)
ARM <- AR[, .(FP_LFMM = median(FP_LFMM), TP_LFMM = median(TP_LFMM),
              FP_best = median(FP_best), TP_best = median(TP_best)), by = arch]
pA3 <- ggplot(AR) +
    geom_segment(aes(x = FP_LFMM, y = TP_LFMM, xend = FP_best, yend = TP_best, colour = direction),
                 linewidth = 0.25, alpha = 0.45,
                 arrow = grid::arrow(length = grid::unit(0.06, "cm"), type = "closed")) +
    geom_segment(data = ARM, aes(x = FP_LFMM, y = TP_LFMM, xend = FP_best, yend = TP_best),
                 linewidth = 1.1, colour = "black",
                 arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")) +
    ARCH_FACET + scale_x_sqrt() + scale_y_sqrt() +
    scale_colour_manual(values = DIR_COL, name = NULL) +
    labs(x = "false positives per seed (sqrt scale)", y = "true hits per seed (sqrt scale)",
         title = "A3  per-seed move from LFMM alone to the 2/3 rule") +
    theme_clinego() + theme(legend.position = "bottom")
save_fig("A3", "A3_arrows_lfmm_to_best", pA3, 13, 5.2,
         "one arrow per seed from LFMM-alone to the 2/3 panel; black = median move")

# A4 -- composition from per-seed MEANS (medians of components do not sum; filed defect)
COMP <- pr[, .(causal = mean(n_causal), linked = mean(n_linked), background = mean(FP),
               total_mean = mean(n), total_median = median(n),
               total_q25 = quantile(n, .25), total_q75 = quantile(n, .75)), by = .(label, arch)]
stopifnot(COMP[, max(abs(causal + linked + background - total_mean))] < 1e-9)
emit(COMP, "A4_composition_means")
CL <- melt(COMP, id.vars = c("label", "arch"), measure.vars = c("causal", "linked", "background"),
           variable.name = "class", value.name = "markers")
CL[, class := factor(class, levels = c("background", "linked", "causal"))]
CLASS_COL <- c(causal = CLINEGO_RETAINED, linked = MINOU[["sage"]], background = CLINEGO_REMOVED)
A4_ORDER <- c(RULES, ORACLE, RANDOM)
pA4 <- ggplot(CL, aes(markers, label)) +
    geom_col(aes(fill = class), width = 0.7, alpha = 0.85) +
    geom_linerange(data = COMP, aes(xmin = total_q25, xmax = total_q75, y = label),
                   inherit.aes = FALSE, linewidth = 0.4, colour = CLINEGO_COL$fg) +
    geom_point(data = COMP, aes(x = total_median, y = label), inherit.aes = FALSE,
               shape = 124, size = 3, colour = CLINEGO_COL$fg) +
    facet_of(cols = "arch", scales = "free_x") +
    scale_y_discrete(limits = rev(A4_ORDER)) +
    scale_fill_manual(values = CLASS_COL, name = NULL, breaks = c("causal", "linked", "background")) +
    labs(x = "markers per seed (bars: mean composition; tick + line: median and IQR of the total)",
         y = NULL, title = "A4  panel composition, causal / linked / background") +
    theme_clinego() + theme(legend.position = "bottom")
save_fig("A4", "A4_composition", pA4, 13, 4.8,
         "stacked mean composition per rule incl. causal loci and random references; median/IQR of panel size")

message("=== B: detection scalar")
# =============================================================================
# B -- one number per seed
# =============================================================================
dist_plot <- function(D, x, xlab, title, order, scales = "free_x") {
    ggplot(D, aes(.data[[x]], label)) +
        geom_boxplot(aes(fill = label), outlier.shape = NA, alpha = 0.55, linewidth = 0.3, width = 0.62) +
        geom_jitter(height = 0.16, width = 0, size = 0.45, alpha = 0.28, colour = CLINEGO_COL$fg) +
        facet_of(cols = "arch", scales = scales) +
        scale_y_discrete(limits = rev(order)) +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        labs(x = xlab, y = NULL, title = title) + theme_clinego()
}

# B1 -- Pareto dominance. Tie rule declared BEFORE computing: A dominates B when
# TP_A >= TP_B and FP_A <= FP_B with at least one strict. Three shares per ordered pair.
PW <- dcast(DETR, seed + arch ~ label, value.var = c("TP", "FP"), sep = "|")
DOM <- rbindlist(lapply(RULES, function(a) rbindlist(lapply(setdiff(RULES, a), function(b) {
    ta <- PW[[paste0("TP|", a)]]; fa <- PW[[paste0("FP|", a)]]
    tb <- PW[[paste0("TP|", b)]]; fb <- PW[[paste0("FP|", b)]]
    a_dom <- ta >= tb & fa <= fb & (ta > tb | fa < fb)
    b_dom <- tb >= ta & fb <= fa & (tb > ta | fb < fa)
    data.table(arch = PW$arch, row = a, col = b, a_dom, b_dom)
}))))
DOMS <- rbind(DOM[, .(n = .N, pct_row = 100 * mean(a_dom), pct_col = 100 * mean(b_dom)),
                  by = .(arch, row, col)],
              DOM[, .(n = .N, pct_row = 100 * mean(a_dom), pct_col = 100 * mean(b_dom)),
                  by = .(row, col)][, arch := "pooled"], use.names = TRUE)
DOMS[, pct_neither := 100 - pct_row - pct_col]
DOMS[, arch := factor(arch, levels = c(ARCH_LABELS, "pooled"))]
DOMS[, `:=`(row = factor(row, levels = RULES), col = factor(col, levels = RULES))]
emit(DOMS, "B1_dominance")
pB1 <- ggplot(DOMS, aes(col, row)) +
    geom_tile(aes(fill = pct_row - pct_col), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f/%.0f/%.0f", pct_row, pct_col, pct_neither)), size = 2.2) +
    facet_wrap(~ arch, nrow = 1) +
    scale_y_discrete(limits = rev(RULES)) +
    scale_fill_gradient2(low = CLINEGO_REMOVED, mid = "white", high = CLINEGO_RETAINED,
                         midpoint = 0, limits = c(-100, 100),
                         name = "row dominates\nminus col dominates (%)") +
    labs(x = "compared against (column)", y = "rule (row)",
         title = "B1  Pareto dominance: % seeds row dominates / col dominates / neither") +
    theme_clinego() + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                            legend.position = "bottom")
save_fig("B1", "B1_dominance", pB1, 15, 5.6,
         "pairwise per-seed Pareto dominance (>= TP and <= FP, one strict), three shares per cell")

# B2 -- precision (causal+linked)/called beside absolute true hits
B2 <- rbind(DETR[, .(seed, label, arch, metric = "precision: true hits / called", value = prec_cl)],
            DETR[, .(seed, label, arch, metric = "true hits (causal + linked)", value = as.numeric(TP))])
B2 <- B2[is.finite(value)]
emit(B2, "B2_precision_and_yield")
pB2 <- ggplot(B2, aes(value, label)) +
    geom_boxplot(aes(fill = label), outlier.shape = NA, alpha = 0.55, linewidth = 0.3, width = 0.62) +
    geom_jitter(height = 0.16, width = 0, size = 0.4, alpha = 0.25, colour = CLINEGO_COL$fg) +
    facet_grid(metric ~ arch, scales = "free_x") +
    scale_y_discrete(limits = rev(RULES)) +
    scale_fill_manual(values = RULE_COL, guide = "none") +
    labs(x = NULL, y = NULL, title = "B2  precision and yield, side by side") + theme_clinego()
save_fig("B2", "B2_precision_and_yield", pB2, 13, 6.4,
         "two distributions per rule: precision (causal+linked)/called and the absolute true-hit count")

save_fig("B3", "B3_f1_causal",
         dist_plot(DETR, "F1", "F1 against causal loci (linked hits count as misses)",
                   "B3  F1 on causal loci", RULES),
         13, 4.2, "F1 = harmonic mean of causal precision and causal recall, per seed")
save_fig("B5", "B5_net_discoveries",
         dist_plot(DETR, "net", "true hits minus false positives per seed",
                   "B5  net discoveries (TP - FP)", RULES) + ZERO,
         13, 4.2, "TP - FP per seed; rewards volume, favours the union")

# B4 -- AUC-PR, single methods only (no combination rule has a p-value)
auc <- fread(file.path(DET, "aucpr_per_seed.tsv"), colClasses = c(seed = "character"))
auc <- mvp_require_rows(auc[seed %in% PRIM$seed], "AUC-PR rows")
stopifnot(nrow(auc) == 3L * nrow(PRIM), !anyNA(auc$aucpr))
auc <- merge(auc[, .(seed, label = method, aucpr, n_testable)], COV, by = "seed")
auc[, label := factor(label, levels = c("LFMM", "RDA", "EMMAX"))]
emit(auc, "B4_aucpr_per_seed")
emit(rbind(auc[, .(n = .N, median = median(aucpr), q25 = quantile(aucpr, .25),
                   q75 = quantile(aucpr, .75)), by = .(label, arch)],
           auc[, .(n = .N, median = median(aucpr), q25 = quantile(aucpr, .25),
                   q75 = quantile(aucpr, .75)), by = label][, arch := "pooled"], use.names = TRUE),
     "B4_aucpr_stats")
save_fig("B4", "B4_aucpr_single_methods",
         dist_plot(auc, "aucpr", "AUC-PR (positives = causal loci, FP = background only)",
                   "B4  AUC-PR, single methods (threshold-free)", c("LFMM", "RDA", "EMMAX")),
         13, 3.4, "threshold-free AUC-PR per single method; RDA ranks on one multivariate column")

message("=== C: offset representation")
# =============================================================================
# C -- delta vs causal loci
# =============================================================================
est_plot <- function(D, rows = NULL, cols = "arch", title) {
    S <- summ_delta(D, c("label", rows, cols))
    ggplot(D, aes(delta, label)) + ZERO +
        geom_jitter(aes(colour = label), height = 0.2, width = 0, size = 0.45, alpha = 0.25) +
        geom_linerange(data = S, aes(xmin = ci_lo, xmax = ci_hi, y = label), inherit.aes = FALSE,
                       linewidth = 1.5, colour = "black") +
        geom_point(data = S, aes(x = median, y = label), inherit.aes = FALSE, size = 2.1,
                   shape = 21, fill = "white", colour = "black", stroke = 0.6) +
        facet_of(rows, cols) +
        scale_y_discrete(limits = rev(OFF_ORDER)) +
        scale_colour_manual(values = RULE_COL, guide = "none") +
        labs(x = "accuracy of the panel minus accuracy of the causal loci (per seed)", y = NULL,
             title = title) + theme_clinego()
}
CS <- rbind(summ_delta(OFFD, c("label", "arch")),
            summ_delta(OFFD, "label")[, arch := "pooled"], use.names = TRUE)
emit(CS, "C_delta_stats")
# Regression tie to journal 15's delta_by_arch.tsv medians (same corpus, same aggregation).
dba <- fread(file.path(J15, "delta_by_arch.tsv"))
tie <- merge(CS[, .(label = as.character(label), arch = as.character(arch), median, n)],
             dba[, .(label, arch = arch_lab, median_j15 = median_delta, n_j15 = n)],
             by = c("label", "arch"))
stopifnot(nrow(tie) == length(OFF_ORDER) * 4L, max(abs(tie$median - tie$median_j15)) < 1e-9,
          all(tie$n == tie$n_j15))
message(sprintf("regression vs journal 15 delta_by_arch: %d medians, max |diff| = %.1e",
                nrow(tie), max(abs(tie$median - tie$median_j15))))

save_fig("C1", "C1_estimation", est_plot(OFFD, title = "C1  delta vs causal loci: seeds, median, 95% CI"),
         13, 4.6, "per-seed deltas + median with 95% bootstrap CI; zero = matches the causal loci")

# C2 -- raincloud
RC <- copy(OFFD)[, ypos := as.numeric(factor(label, levels = rev(OFF_ORDER)))]
pC2 <- ggplot(RC, aes(x = delta, y = ypos)) + ZERO +
    stat_halfeye(aes(fill = label, group = label), orientation = "horizontal", adjust = 0.8,
                 height = 0.55, justification = -0.18, .width = 0, point_colour = NA,
                 slab_alpha = 0.65) +
    geom_boxplot(aes(group = label), orientation = "y", width = 0.14, outlier.shape = NA,
                 alpha = 0.6, linewidth = 0.3) +
    geom_point(aes(y = ypos - 0.22, colour = label), position = position_jitter(height = 0.07, width = 0,
               seed = BOOT_SEED), size = 0.35, alpha = 0.3) +
    ARCH_FACET +
    scale_y_continuous(breaks = seq_along(OFF_ORDER), labels = rev(OFF_ORDER)) +
    scale_fill_manual(values = RULE_COL, guide = "none") +
    scale_colour_manual(values = RULE_COL, guide = "none") +
    labs(x = "accuracy of the panel minus accuracy of the causal loci (per seed)", y = NULL,
         title = "C2  delta vs causal loci, raincloud") + theme_clinego()
save_fig("C2", "C2_raincloud", pC2, 13, 5.2, "density slab + box + seeds per rule")

# C3 -- within-seed rank among the six obtainable rules
RK <- merge(CJ(seed = PRIM$seed, label = RULES), OFFD[, .(seed, label = as.character(label), accuracy)],
            by = c("seed", "label"), all.x = TRUE)
RK[, rank := frank(-accuracy, ties.method = "min", na.last = "keep"), by = seed]
RANK_LAB <- c("1st", "2nd", "3rd", "4th", "5th", "6th", "no panel (< 3 SNPs)")
RK[, rank_lab := factor(fifelse(is.na(rank), RANK_LAB[7], RANK_LAB[pmin(rank, 6L)]), levels = rev(RANK_LAB))]
RK <- merge(RK, COV[, .(seed, arch)], by = "seed")
RKS <- rbind(RK[, .N, by = .(label, arch, rank_lab)][, pct := 100 * N / sum(N), by = .(label, arch)],
             RK[, .N, by = .(label, rank_lab)][, pct := 100 * N / sum(N), by = label][, arch := "pooled"],
             use.names = TRUE)
RKS[, arch := factor(arch, levels = c(ARCH_LABELS, "pooled"))]
emit(RK[, .(seed, label, arch, accuracy, rank, rank_lab)], "C3_rank_per_seed")
emit(RKS, "C3_rank_share")
RANK_COL <- setNames(c(grDevices::colorRampPalette(c(CLINEGO_RETAINED, "#EFEFEF", CLINEGO_REMOVED))(6),
                       MINOU[["grey"]]), RANK_LAB)
pC3 <- ggplot(RKS, aes(pct, label, fill = rank_lab)) +
    geom_col(width = 0.72) +
    facet_wrap(~ arch, nrow = 1) +
    scale_y_discrete(limits = rev(RULES)) +
    scale_fill_manual(values = RANK_COL, breaks = RANK_LAB, name = "within-seed rank") +
    labs(x = "% of seeds", y = NULL, title = "C3  how often each rule ranks 1st ... 6th on the same seed") +
    theme_clinego() + theme(legend.position = "bottom")
save_fig("C3", "C3_rank_share", pC3, 14, 4.6,
         "share of seeds at each within-seed accuracy rank among the six obtainable rules")

# C4 -- ECDF
pC4 <- ggplot(OFFD, aes(delta, colour = label, linetype = label == RANDOM)) + ZERO +
    stat_ecdf(linewidth = 0.55) + ARCH_FACET +
    scale_colour_manual(values = RULE_COL, name = NULL) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
    labs(x = "accuracy of the panel minus accuracy of the causal loci (per seed)",
         y = "cumulative share of seeds", title = "C4  cumulative distribution of delta") +
    theme_clinego() + theme(legend.position = "bottom")
save_fig("C4", "C4_ecdf", pC4, 13, 4.8, "ECDF of per-seed delta; curve further right = better")

# C5 -- reach tile (% seeds at or above the causal loci)
CS5 <- copy(CS)[, arch := factor(arch, levels = c(ARCH_LABELS, "pooled"))]
pC5 <- ggplot(CS5, aes(arch, label)) +
    geom_tile(aes(fill = reach_pct), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f%%", reach_pct)), size = 3) +
    scale_y_discrete(limits = rev(OFF_ORDER)) +
    scale_fill_gradient2(low = CLINEGO_REMOVED, mid = "white", high = CLINEGO_RETAINED,
                         midpoint = 50, limits = c(0, 100), name = "% seeds\ndelta >= 0") +
    labs(x = NULL, y = NULL, title = "C5  share of seeds where the panel reaches the causal loci") +
    theme_clinego()
save_fig("C5", "C5_reach_tile", pC5, 7.5, 4.6, "% seeds with delta >= 0 per rule x genic level")

message("=== D: splits")
# =============================================================================
# D -- the same two views under six splits
# =============================================================================
SPLITS <- list(
    D1 = list(name = "pooled",                         rows = NULL,     cols = NULL,                w = 7,  h = 4.6),
    D2 = list(name = "genic level",                    rows = NULL,     cols = "arch",              w = 13, h = 4.6),
    D3 = list(name = "genic level x selection regime", rows = "regime", cols = "arch",              w = 13, h = 7.6),
    D4 = list(name = "genic level x pleiotropy",       rows = "pleio",  cols = "arch",              w = 13, h = 7.6),
    D5 = list(name = "all 12 cells",                   rows = "arch",   cols = c("pleio", "regime"), w = 17, h = 10),
    D6 = list(name = "demography block",               rows = NULL,     cols = "block",             w = 19, h = 4.8))
for (code in names(SPLITS)) {
    s <- SPLITS[[code]]
    emit(summ_delta(OFFD, c("label", s$rows, s$cols)), paste0(code, "_offset_stats"))
    save_fig(paste0(code, "-off"), paste0(code, "_offset"),
             est_plot(OFFD, s$rows, s$cols, sprintf("%s  offset, split: %s", code, s$name)),
             s$w, s$h, sprintf("C1 view (delta vs causal loci) split by %s", s$name))
    emit(det_stats(DETR, c("label", s$rows, s$cols)), paste0(code, "_detection_stats"))
    save_fig(paste0(code, "-det"), paste0(code, "_detection"),
             plane_plot(DETR, s$rows, s$cols, sprintf("%s  detection, split: %s", code, s$name)),
             s$w, s$h, sprintf("A1 view (TP-vs-FP plane) split by %s", s$name))
}

M <- rbindlist(MANIFEST)
fwrite(M, file.path(OUT, "gallery_manifest.tsv"), sep = "\t")
message(sprintf("OK: %d figures, arm [%s] -> %s", nrow(M), mvp_arm_label(), OUT))
