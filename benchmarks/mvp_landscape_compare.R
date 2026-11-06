#!/usr/bin/env Rscript
# =============================================================================
# mvp_landscape_compare.R -- SS-Clines | SS-Mtn | Est-Clines side by side (plan nifty-prancing-grove,
# Phase 9). Plotting and tables only: every number is read from a landscape's own scored tables,
# written by the same scoring code for all three (SS-Clines: the published LD-window arm; SS-Mtn and
# Est-Clines: work/mvp_land/post.sh + figs.sh, whose copies tie to SS-Clines byte for byte).
#
# Per landscape (600 2-trait replicates each):
#   detection   <gallery>/detection_per_seed.tsv   option B (causal or within the replicate's LD window)
#   AUC-PR      <gallery>/B4_aucpr_per_seed.tsv    RDA without structure correction
#               <remeasure>/onefit_rda2x_ldb/rank_metrics.tsv   RDA with structure correction
#   offset      <gallery>/offset_delta_per_seed.tsv  accuracy minus the causal loci's (median of 3 engines)
#   rank        <rank>/rank_long.tsv               within-replicate rank of the six sets, median of 3
#               engines; ties split evenly (benchmarks/mvp_rank_ties.R)
#   LD window   work/mvp_land/ld_window_by_replicate.tsv  (clumping_distance_bp, the window the sets used)
#   covariates  benchmarks/mvp_seeds.tsv
# GATE (stops): each landscape's medians recomputed here equal the ones its own figure set published
# (<ms>/numbers.tsv: Figure 2a detection, S2 AUC-PR, Figure 2b least accurate; <gallery>/C_delta_stats.tsv
# offset) to 1e-12, and every table holds exactly the landscape's 600 seeds.
# Tests: Kruskal-Wallis across the three landscapes, set by set (pooled replicates), Holm within each
# metric family (the statistics fixed 2026-10-26). The pooled 1,800 appear only as a labelled extra row.
# Outputs (OUT_DIR, default benchmarks/mvp_eval/landscape_compare/):
#   summary.tsv     landscape x family x item x stratum: n, median, q25, q75 (+ lo/hi for least accurate)
#   tests.tsv       family x item: Kruskal-Wallis H, df, p, Holm p, per-landscape medians
#   per_seed.tsv    the long per-replicate table every number and figure is drawn from
#   counts.tsv      empty and under-3 sets per landscape
#   gate.tsv        the reproduction check, one row per compared number
#   L1..L7 figures (SVG + PNG, ms_style)
# Usage: Rscript benchmarks/mvp_landscape_compare.R   (env PIPELINE_ROOT, OUT_DIR)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OUT  <- Sys.getenv("OUT_DIR", file.path(EVAL, "landscape_compare"))
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_rank_ties.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TOL <- 1e-12

LANDS <- data.table(
    land  = c("SS-Clines", "SS-Mtn", "Est-Clines"),
    lkey  = c("ssclines", "ssmtn", "estclines"),
    gal   = c("figures_ld_gallery14", "figures_land_ssmtn_gallery", "figures_land_estclines_gallery"),
    rank  = c("figures_ld_offset14/rank", "figures_land_ssmtn_offset15_ssmtn/rank",
              "figures_land_estclines_offset15_estclines/rank"),
    ms    = c("figures_ld_ms", "figures_land_ssmtn_ms", "figures_land_estclines_ms"),
    rem   = c("remeasure600", "remeasure600_ssmtn", "remeasure600_estclines"),
    seeds = c("work/ld_phase3_2026-10-26/scripts/seeds600.txt", "work/mvp_land/seeds_ssmtn.txt",
              "work/mvp_land/seeds_estclines.txt"))
LAND_LEV  <- LANDS$land
LAND_COL  <- setNames(c(MINOU[["navy"]], MINOU[["teal"]], MINOU[["amber"]]), LAND_LEV)
SETS      <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
SET_LAB   <- c("2/3 methods" = "2/3 set", "1/3 methods" = "1/3 set", "3/3 methods" = "3/3 set",
               LFMM = "LFMM", RDA = "RDA", EMMAX = "EMMAX")
ARCH_LEV  <- c("oligogenic", "moderately polygenic", "highly polygenic")
ARCH_OF   <- c("oliogenic" = "oligogenic", "mod-polygenic" = "moderately polygenic",
               "highly-polygenic" = "highly polygenic")
rd <- function(...) fread(file.path(...), colClasses = c(seed = "character"))
man <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
LDW <- fread(file.path(ROOT, "work/mvp_land/ld_window_by_replicate.tsv"), colClasses = c(seed = "character"))

# ------------------------------------------------------------------ per-seed long table
PS <- list(); GATE <- list(); CNT <- list()
gate <- function(L, what, mine, pub) {
    x <- merge(mine, pub, by = setdiff(intersect(names(mine), names(pub)), c("v_mine", "v_pub")))
    if (nrow(x) != nrow(mine) || nrow(x) == 0L) stop(L, " ", what, ": ", nrow(mine), " rows recomputed, ", nrow(x), " matched")
    x[, `:=`(landscape = L, what = what, abs_diff = abs(v_mine - v_pub))]
    if (any(!is.finite(x$abs_diff)) || max(x$abs_diff) > TOL)
        stop(L, " ", what, ": recomputed medians differ from the published ones (max |diff| ", max(x$abs_diff), ")")
    GATE[[length(GATE) + 1L]] <<- x[, .(n_compared = .N, max_abs_diff = max(abs_diff)), by = .(landscape, what)]
}
# na.rm: precision is undefined for an empty set (the published medians skip those replicates too)
med_by <- function(D, val) rbind(D[, .(v_mine = median(get(val), na.rm = TRUE)), by = .(label, arch)],
                                 D[, .(v_mine = median(get(val), na.rm = TRUE)), by = .(label)][, arch := "pooled"])
for (i in seq_len(nrow(LANDS))) {
    L <- LANDS$land[i]; G <- file.path(EVAL, LANDS$gal[i]); R <- file.path(EVAL, LANDS$rank[i])
    M <- file.path(EVAL, LANDS$ms[i]); RM <- file.path(EVAL, LANDS$rem[i])
    SEEDS <- trimws(readLines(file.path(ROOT, LANDS$seeds[i]))); SEEDS <- SEEDS[nzchar(SEEDS)]
    stopifnot(length(SEEDS) == 600L, !anyDuplicated(SEEDS))
    own <- function(D, what) {
        if (!setequal(unique(D$seed), SEEDS)) stop(L, " ", what, ": seeds differ from ", LANDS$seeds[i])
        D
    }
    NUM <- fread(file.path(M, "numbers.tsv"), colClasses = "character")
    # detection (option B)
    DET <- own(rd(G, "detection_per_seed.tsv")[label %in% SETS], "detection")
    for (v in c("TP", "FP", "prec_cl")) {
        pub <- NUM[figure == "Figure 2a / S4 / S5 detection" & metric == paste0(v, "_median"),
                   .(label = rule, arch = stratum, v_pub = as.numeric(value))]
        gate(L, paste("detection", v), med_by(DET, v), pub)
    }
    PS[[length(PS) + 1L]] <- melt(DET[, .(landscape = L, seed, item = label, arch, TP, FP, precision = prec_cl)],
                                  measure.vars = c("TP", "FP", "precision"), variable.name = "family")[!is.na(value)]
    CNT[[L]] <- DET[, .(landscape = L, empty = sum(n == 0), under3 = sum(n < 3), n = .N), by = .(item = label)]
    # AUC-PR (RDA uncorrected) + RDA corrected
    AU <- own(rd(G, "B4_aucpr_per_seed.tsv"), "AUC-PR")
    gate(L, "AUC-PR", med_by(AU, "aucpr"),
         NUM[figure == "Figure S2 AUC-PR per method (RDA without structure correction)" & metric == "median",
             .(label = rule, arch = stratum, v_pub = as.numeric(value))])
    RC <- own(rd(RM, "onefit_rda2x_ldb", "rank_metrics.tsv")[method == "RDA"], "RDA corrected")
    RC <- merge(RC, AU[label == "RDA", .(seed, arch)], by = "seed")
    PS[[length(PS) + 1L]] <- rbind(AU[, .(landscape = L, seed, item = label, arch, family = "AUC-PR", value = aucpr)],
                                   RC[, .(landscape = L, seed, item = "RDA corrected", arch, family = "AUC-PR", value = aucpr)])
    # offset accuracy delta (sets under 3 SNPs have no offset: 3/3 rows may be missing)
    OD <- rd(G, "offset_delta_per_seed.tsv")[label %in% SETS]
    stopifnot(all(OD$seed %in% SEEDS))
    CD <- rd(G, "C_delta_stats.tsv")
    gate(L, "offset delta", med_by(OD, "delta"), CD[, .(label, arch, v_pub = median)])
    PS[[length(PS) + 1L]] <- OD[, .(landscape = L, seed, item = label, arch, family = "offset delta", value = delta)]
    # within-replicate rank, median of 3 engines: rank and least-accurate weight (ties split evenly)
    RL <- own(rd(R, "rank_long.tsv")[engine == "median of 3"], "rank")
    RX <- rank_split(RL, "seed")
    IND <- rank_indicators(RX, "seed", n_set = length(SETS))
    IND <- merge(IND, RL[, .(seed, panel, arch = arch_lab)], by = c("seed", "panel"))
    least_pub <- NUM[metric == "bottom1" & engine == "Median of the three",
                     .(label = rule, arch = stratum, v_pub = as.numeric(value))]
    gate(L, "least accurate %", rbind(IND[, .(v_mine = 100 * mean(bottom1)), by = .(label = panel, arch)],
                                      IND[, .(v_mine = 100 * mean(bottom1)), by = .(label = panel)][, arch := "pooled"]),
         least_pub)
    PS[[length(PS) + 1L]] <- rbind(IND[, .(landscape = L, seed, item = panel, arch, family = "rank", value = mean_rank)],
                                   IND[, .(landscape = L, seed, item = panel, arch, family = "least accurate", value = 100 * bottom1)])
    # LD window and covariates
    LW <- own(LDW[landscape == LANDS$lkey[i]], "LD window")
    CV <- own(man[seed %in% SEEDS], "covariates")
    CV[, arch := unname(ARCH_OF[arch_level])]
    PS[[length(PS) + 1L]] <- rbind(
        merge(LW[, .(seed, value = clumping_distance_bp)], CV[, .(seed, arch)], by = "seed")[
            , .(landscape = L, seed, item = "LD window (bp)", arch, family = "LD window", value)],
        melt(CV[, .(landscape = L, seed, arch, r2_pc1_temp, r2_pc1_sal, meanFst, final_LA, n_causal_maf01)],
             measure.vars = c("r2_pc1_temp", "r2_pc1_sal", "meanFst", "final_LA", "n_causal_maf01"),
             variable.name = "item")[, family := "covariate"][, .(landscape, seed, item = as.character(item), arch, family, value)])
    message(sprintf("%s: %d seeds, gate passed", L, length(SEEDS)))
}
PS <- rbindlist(PS, use.names = TRUE)
PS[, `:=`(landscape = factor(landscape, levels = LAND_LEV), family = as.character(family),
          value = as.numeric(value))]
stopifnot(!anyNA(PS$value), !anyNA(PS$arch), all(PS$arch %in% ARCH_LEV))
fwrite(PS, file.path(OUT, "per_seed.tsv"), sep = "\t")
GATE <- rbindlist(GATE); fwrite(GATE, file.path(OUT, "gate.tsv"), sep = "\t")
CNT <- rbindlist(CNT); fwrite(CNT, file.path(OUT, "counts.tsv"), sep = "\t")
message(sprintf("GATE PASSED: %d recomputed medians equal the published ones (max |diff| %.1e)",
                sum(GATE$n_compared), max(GATE$max_abs_diff)))

# ------------------------------------------------------------------ summary + tests
summ <- function(D) D[, .(n = .N, median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75)),
                      by = .(landscape, family, item, stratum)]
POOL <- copy(PS)[, landscape := "pooled 1,800 (extra row)"]
SUMM <- rbind(summ(PS[, .(landscape = as.character(landscape), family, item, value, stratum = "pooled")]),
              summ(PS[, .(landscape = as.character(landscape), family, item, value, stratum = arch)]),
              summ(POOL[, .(landscape, family, item, value, stratum = "pooled")]))
# least accurate: share (mean of the per-seed weight) with a bootstrap CI over replicates
set.seed(20261106L)
boot <- function(x, B = 2000L) unname(quantile(colMeans(matrix(sample(x, length(x) * B, TRUE), ncol = B)), c(.025, .975)))
LA <- PS[family == "least accurate", {ci <- boot(value); .(share = mean(value), lo = ci[1], hi = ci[2], n = .N)},
         by = .(landscape, item)]
SUMM <- merge(SUMM, LA[, .(landscape = as.character(landscape), family = "least accurate", item, stratum = "pooled",
                           share = share, lo, hi)], by = c("landscape", "family", "item", "stratum"), all.x = TRUE)
setorder(SUMM, family, item, stratum, landscape)
fwrite(SUMM, file.path(OUT, "summary.tsv"), sep = "\t")

TESTS <- PS[, {
    k <- kruskal.test(value ~ landscape)
    md <- tapply(value, landscape, median)
    .(n = .N, H = unname(k$statistic), df = unname(k$parameter), p = k$p.value,
      median_ss_clines = md[["SS-Clines"]], median_ss_mtn = md[["SS-Mtn"]], median_est_clines = md[["Est-Clines"]])
}, by = .(family, item)]
TESTS[, p_holm := p.adjust(p, "holm"), by = family]
setcolorder(TESTS, c("family", "item", "n", "H", "df", "p", "p_holm"))
fwrite(TESTS, file.path(OUT, "tests.tsv"), sep = "\t")

# ------------------------------------------------------------------ figures
SET_COL <- c("2/3 methods" = MINOU[["teal"]], "1/3 methods" = MINOU[["sage"]], "3/3 methods" = MINOU[["navy"]],
             LFMM = MINOU[["red"]], RDA = MINOU[["amber"]], EMMAX = MS_ORANGE)
setf <- function(x) factor(x, levels = SETS, labels = SET_LAB[SETS])
save <- function(stem, p, w, h) { ms_save(file.path(OUT, stem), p, w, h); message("  ", stem) }

# L1 detection plane: median TP vs FP with IQR, per landscape x architecture
D1 <- PS[family %in% c("TP", "FP"), .(m = median(value), lo = quantile(value, .25), hi = quantile(value, .75)),
         by = .(landscape, arch, item, family)]
D1 <- dcast(D1, landscape + arch + item ~ family, value.var = c("m", "lo", "hi"))
D1[, `:=`(arch = factor(arch, levels = ARCH_LEV), set = setf(item))]
p1 <- ggplot(D1, aes(m_FP, m_TP, colour = item)) +
    geom_errorbar(aes(ymin = lo_TP, ymax = hi_TP), width = 0, linewidth = 0.3, alpha = 0.6) +
    geom_errorbar(aes(xmin = lo_FP, xmax = hi_FP), orientation = "y", width = 0, linewidth = 0.3, alpha = 0.6) +
    geom_point(size = 1.4) +
    facet_grid(rows = vars(arch), cols = vars(landscape), scales = "free") +
    scale_colour_manual(values = SET_COL, labels = SET_LAB, name = NULL) +
    labs(x = "False positives (neutral SNPs), median and IQR", y = "True positives (causal or linked SNPs)") + theme_ms()
save("L1_detection_plane", p1, 180, 120)

# L2 rank distribution (median of 3 engines), share at each rank, per landscape
RK <- rbindlist(lapply(LAND_LEV, function(L) {
    RL <- rd(EVAL, LANDS[land == L, rank], "rank_long.tsv")[engine == "median of 3"]
    rank_split(RL, "seed")[, .(share = 100 * sum(w) / uniqueN(RL$seed)), by = .(panel, r)][, landscape := L]
}))
RK[, `:=`(landscape = factor(landscape, levels = LAND_LEV), set = factor(SET_LAB[panel], levels = rev(SET_LAB[SETS])),
          r = factor(r, levels = 6:1))]
RANK_COL <- setNames(colorRampPalette(c(MINOU[["teal"]], "#e8eef0", MINOU[["red"]]))(6), 1:6)
p2 <- ggplot(RK, aes(share, set, fill = r)) + geom_col(width = 0.7) + facet_wrap(~landscape, nrow = 1) +
    scale_fill_manual(values = RANK_COL, breaks = 1:6, name = "rank (1 = most accurate)") +
    labs(x = "% replicates at each rank (median of the three offset engines)", y = NULL) + theme_ms() +
    theme(legend.position = "bottom")
save("L2_rank_distribution", p2, 180, 70)

# L3 least accurate share with bootstrap CI, by set, landscapes side by side
LA[, `:=`(set = setf(item))]
p3 <- ggplot(LA, aes(share, set, colour = landscape)) +
    geom_vline(xintercept = 100 / 6, linetype = "dashed", colour = MS_REF, linewidth = 0.35) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.25, linewidth = 0.4) +
    scale_colour_manual(values = LAND_COL, name = NULL) + scale_y_discrete(limits = rev) +
    labs(x = "% replicates in which the set is the least accurate (95 % bootstrap CI; dashed = 1 in 6)", y = NULL) +
    theme_ms() + theme(legend.position = "bottom")
save("L3_least_accurate", p3, 120, 80)

# L4 AUC-PR by method x architecture, landscapes side by side (RDA with and without correction)
A4 <- PS[family == "AUC-PR"][, `:=`(arch = factor(arch, levels = ARCH_LEV),
          method = factor(item, levels = c("LFMM", "EMMAX", "RDA", "RDA corrected"),
                          labels = c("LFMM", "EMMAX", "RDA uncorrected", "RDA corrected")))]
p4 <- ggplot(A4, aes(method, value, fill = landscape)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, width = 0.7, position = position_dodge(width = 0.8)) +
    facet_wrap(~arch, nrow = 1, scales = "free_y") + scale_fill_manual(values = LAND_COL, name = NULL) +
    labs(x = NULL, y = "AUC-PR (causal loci vs neutral SNPs)") + theme_ms() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
save("L4_aucpr", p4, 180, 75)

# L5 offset accuracy relative to the causal loci, by set (median, IQR), landscapes side by side
O5 <- PS[family == "offset delta", .(m = median(value), lo = quantile(value, .25), hi = quantile(value, .75)),
         by = .(landscape, item)][, set := setf(item)]
p5 <- ggplot(O5, aes(m, set, colour = landscape)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = MS_REF, linewidth = 0.35) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.25, linewidth = 0.4) +
    scale_colour_manual(values = LAND_COL, name = NULL) + scale_y_discrete(limits = rev) +
    labs(x = "Δ accuracy (candidate set − causal loci), median and IQR", y = NULL) +
    theme_ms() + theme(legend.position = "bottom")
save("L5_offset_delta", p5, 120, 80)

# L6 LD window per replicate
L6 <- PS[family == "LD window"]
p6 <- ggplot(L6, aes(value, fill = landscape)) + geom_histogram(bins = 40, colour = NA, alpha = 0.85) +
    facet_wrap(~landscape, ncol = 1, scales = "free_y") + scale_x_log10() +
    scale_fill_manual(values = LAND_COL, guide = "none") +
    labs(x = "LD window used by the sets (bp, log scale)", y = "replicates") + theme_ms()
save("L6_ld_window", p6, 90, 90)

# L7 realised covariates
C7 <- PS[family == "covariate"]
p7 <- ggplot(C7, aes(landscape, value, fill = landscape)) + geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    facet_wrap(~item, nrow = 1, scales = "free_y") + scale_fill_manual(values = LAND_COL, guide = "none") +
    labs(x = NULL, y = NULL) + theme_ms() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save("L7_covariates", p7, 180, 60)

message(sprintf("OK: %d per-seed rows, %d summary rows, %d tests -> %s", nrow(PS), nrow(SUMM), nrow(TESTS), OUT))
