#!/usr/bin/env Rscript
# =============================================================================
# mvp_lind_compare.R -- put CLINE-GO's garden-sweep numbers next to the published ones.
#
# Every published value below is quoted from Lind & Lotterhos 2025 (Mol Ecol Resour 25(4):e14008)
# and is hard-coded here with its source string, so the comparison table can never drift from what
# the paper actually says.
#
# Four comparisons, all on their metric (Kendall's tau, negative = good):
#   1. tau vs degree of local adaptation, against their two published anchors
#   2. marker-set comparison: their three win rates and their "<3% median gain"
#   3. climate-novelty decline curve
#   4. method ranking within strata
#
# Usage:  Rscript mvp_lind_compare.R [--outdir=DIR]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args   <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt    <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
OUTDIR <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset09"))

# ---- published values, verbatim ---------------------------------------------
PUB <- data.table(
    quantity = c("tau at LA ~ 0.5 (best methods)",
                 "tau at LA ~ 0.3 (best methods)",
                 "median performance gain, adaptive vs all/neutral",
                 "% of models where adaptive beats neutral",
                 "% of models where adaptive beats all",
                 "% of models where all beats neutral"),
    published = c("-0.6", "-0.2", "< 3%", "68%", "67%", "74%"),
    source = "Lind & Lotterhos 2025, Mol Ecol Resour 25(4):e14008")

GP <- fread(file.path(OUTDIR, "garden_performance.tsv"), colClasses = c("seed" = "character"))
LAND <- GP[garden_type == "landscape" & control == FALSE]
NOV  <- GP[garden_type == "novelty"   & control == FALSE]

# Performance in their orientation: they invert the tau axis so that "higher = better".
LAND[, perf := -tau]; NOV[, perf := -tau]

# ---- 1. tau vs local adaptation ---------------------------------------------
TL <- LAND[, .(tau_med = median(tau, na.rm = TRUE),
               tau_q25 = quantile(tau, 0.25, na.rm = TRUE),
               tau_q75 = quantile(tau, 0.75, na.rm = TRUE),
               n_gardens = uniqueN(garden_id)),
           by = .(seed, final_LA, arch, method_label, marker_set)]
fwrite(TL, file.path(OUTDIR, "tau_vs_LA.tsv"), sep = "\t")

anchors <- data.table(final_LA = c(0.3, 0.5), tau = c(-0.2, -0.6))
p1 <- ggplot(TL[marker_set %in% c("adaptive", "all", "neutral")],
             aes(final_LA, tau_med, colour = method_label, shape = marker_set)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(size = 2.4, alpha = 0.85) +
    geom_point(data = anchors, aes(final_LA, tau), inherit.aes = FALSE,
               shape = 4, size = 5, stroke = 1.4, colour = "black") +
    labs(x = "degree of local adaptation (final_LA, deposit)",
         y = "median Kendall's tau across landscape gardens",
         title = "CLINE-GO offsets on the Lind & Lotterhos scale",
         subtitle = "black crosses = their published anchors (-0.2 at LA 0.3, -0.6 at LA 0.5)") +
    theme_minimal(base_size = 11)
ggsave(file.path(OUTDIR, "tau_vs_LA.png"), p1, width = 9, height = 5.5, dpi = 150)
ggsave(file.path(OUTDIR, "tau_vs_LA.svg"), p1, width = 9, height = 5.5)

# ---- 2. marker-set comparison ------------------------------------------------
# Their unit of comparison is a model: here (seed, method, garden). Two marker sets are compared
# on the same unit, so the contrast is strictly within-garden and within-method.
W <- dcast(LAND[marker_set %in% c("all", "adaptive", "neutral")],
           seed + method_label + garden_id ~ marker_set, value.var = "perf")
W <- W[complete.cases(W)]

pair_stats <- function(a, b, lbl) {
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]; b <- b[ok]
    data.table(comparison = lbl, n_models = length(a),
               pct_first_wins = 100 * mean(a > b),
               median_gain_pct = 100 * median((a - b) / abs(b)),
               median_perf_first = median(a), median_perf_second = median(b))
}
MS <- rbindlist(list(
    pair_stats(W$adaptive, W$neutral, "adaptive > neutral"),
    pair_stats(W$adaptive, W$all,     "adaptive > all"),
    pair_stats(W$all,      W$neutral, "all > neutral")))
fwrite(MS, file.path(OUTDIR, "marker_set_comparison.tsv"), sep = "\t")

# ---- 2b. OUR axis: GEA-derived panels, which their design never tested -------------------
# Their "adaptive" panel is the true QTN set, i.e. an oracle. The practical question is how
# close a GEA-identified panel gets to it, across the precision-recall span:
#   gea_union  (any of LFMM/RDA/EMMAX)  -> maximum recall
#   gea_best   (>= 2 of 3, 5 kb window) -> the journal-07 anchor scheme
#   gea_strict (all 3 agree)            -> maximum precision
WG <- dcast(LAND, seed + method_label + garden_id ~ marker_set, value.var = "perf")
WG <- WG[complete.cases(WG[, .(all, adaptive, neutral)])]
gea_rows <- list()
for (g in intersect(c("gea_union", "gea_best", "gea_strict", "random_matched"), names(WG))) {
    for (ref in intersect(c("adaptive", "all", "neutral", "random_matched"), names(WG))) {
        if (g == ref) next
        ok <- is.finite(WG[[g]]) & is.finite(WG[[ref]])
        if (!any(ok)) next
        gea_rows[[length(gea_rows) + 1L]] <- data.table(
            panel = g, reference = ref, n_models = sum(ok),
            pct_panel_wins = 100 * mean(WG[[g]][ok] > WG[[ref]][ok]),
            median_gain_pct = 100 * median((WG[[g]][ok] - WG[[ref]][ok]) / abs(WG[[ref]][ok])))
    }
}
if (length(gea_rows)) {
    GEA <- rbindlist(gea_rows)
    fwrite(GEA, file.path(OUTDIR, "gea_panels_vs_refs.tsv"), sep = "\t")
    message("\n=== GEA-derived panels vs the reference panels (their design has no equivalent) ===")
    print(GEA[, .(panel, reference, n_models,
                  pct_panel_wins = round(pct_panel_wins, 1),
                  median_gain_pct = round(median_gain_pct, 1))])
}

MSm <- rbindlist(lapply(unique(W$method_label), function(m) {
    w <- W[method_label == m]
    r <- rbindlist(list(pair_stats(w$adaptive, w$neutral, "adaptive > neutral"),
                        pair_stats(w$adaptive, w$all,     "adaptive > all"),
                        pair_stats(w$all,      w$neutral, "all > neutral")))
    r[, method_label := m][]
}))
fwrite(MSm, file.path(OUTDIR, "marker_set_comparison_by_method.tsv"), sep = "\t")

# ---- 3. climate-novelty decline ---------------------------------------------
NV <- NOV[, .(tau_med = median(tau, na.rm = TRUE), n = .N,
              offset_cv_med = median(offset_cv, na.rm = TRUE)),
          by = .(dev, method_label, marker_set)]
fwrite(NV, file.path(OUTDIR, "novelty_curve.tsv"), sep = "\t")

# ---- 3b. the novelty curve read against BOTH baselines -----------------------
# Their deviation grid starts at dev = 0, which is a single synthetic garden sitting exactly at
# the climate CENTRE. Every deme is roughly equidistant from it, so the predicted offsets bunch
# together and ranking them is mostly noise -- a floor effect, not a method failure. Reading the
# curve against that point therefore exaggerates any apparent "improvement" with novelty.
#
# The honest reference is the median over the 30 WITHIN-LANDSCAPE gardens: real, varied
# environments in which the demes genuinely differ. Both readings are reported, per method and
# panel, together with the offset dispersion that decides which one is meaningful.
LAND_BASE <- LAND[, .(tau_landscape_med = median(tau, na.rm = TRUE),
                      offset_cv_landscape = median(offset_cv, na.rm = TRUE)),
                  by = .(method_label, marker_set)]
DEV0 <- NV[dev == 0, .(method_label, marker_set, tau_dev0 = tau_med,
                       offset_cv_dev0 = offset_cv_med)]
NVB <- merge(merge(NV[dev > 0], DEV0, by = c("method_label", "marker_set")),
             LAND_BASE, by = c("method_label", "marker_set"))
NVB[, `:=`(delta_vs_dev0      = tau_med - tau_dev0,
           delta_vs_landscape = tau_med - tau_landscape_med)]
setorder(NVB, method_label, marker_set, dev)
fwrite(NVB, file.path(OUTDIR, "novelty_baselines.tsv"), sep = "\t")

message("\n=== novelty curve against BOTH baselines (adaptive panel) ===")
print(NVB[marker_set == "adaptive",
          .(method_label, dev,
            tau = round(tau_med, 3),
            vs_dev0 = round(delta_vs_dev0, 3),
            vs_landscape = round(delta_vs_landscape, 3),
            offset_cv = round(offset_cv_med, 2))])
message("\n=== offset dispersion: dev0 vs landscape vs far novelty (is dev0 a floor?) ===")
print(rbind(
    LAND_BASE[, .(garden_type = "landscape (median of 30)", method_label, marker_set,
                  offset_cv = round(offset_cv_landscape, 3))],
    DEV0[, .(garden_type = "novelty dev = 0 (climate centre)", method_label, marker_set,
             offset_cv = round(offset_cv_dev0, 3))],
    NV[dev == max(dev), .(garden_type = "novelty dev = 5.88", method_label, marker_set,
                          offset_cv = round(offset_cv_med, 3))]
)[marker_set == "adaptive"])

p2 <- ggplot(NV[marker_set %in% c("adaptive", "all", "neutral")],
             aes(dev, tau_med, colour = method_label, linetype = marker_set)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    labs(x = "climate novelty (SD from the climate centre, their deviation grid)",
         y = "median Kendall's tau", title = "Performance vs climate novelty",
         subtitle = "their finding: performance declines as gardens become more novel") +
    theme_minimal(base_size = 11)
ggsave(file.path(OUTDIR, "novelty_curve.png"), p2, width = 9, height = 5.5, dpi = 150)
ggsave(file.path(OUTDIR, "novelty_curve.svg"), p2, width = 9, height = 5.5)

# ---- 4. method ranking -------------------------------------------------------
MR <- LAND[, .(tau_med = median(tau, na.rm = TRUE), n_gardens = .N),
           by = .(arch, method_label, marker_set)]
setorder(MR, arch, marker_set, tau_med)
fwrite(MR, file.path(OUTDIR, "method_rank.tsv"), sep = "\t")

# ---- 5. stratification by every simulation property we hold ------------------
# The deposit's own per-seed metadata (benchmarks/mvp_seeds.tsv) is the only legitimate
# stratifier here: it was computed by the authors, not derived by us from the same data we
# are scoring. Anything we derived ourselves would risk circularity.
MAN <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
             colClasses = c("seed" = "character"))
LAND <- merge(LAND, MAN[, .(seed, architecture, arch_level, demog_level, n_traits,
                            ispleiotropy, k_best, meanFst, r2_pc1_temp, r2_pc1_sal,
                            n_snps, n_causal_maf01, final_LA)],
              by = "seed", all.x = TRUE, suffixes = c("", "_man"))

strat_tbl <- function(by_col) {
    x <- LAND[, .(n_seeds = uniqueN(seed), n_models = .N,
                  tau_med = round(median(tau, na.rm = TRUE), 3)),
              by = c(by_col, "method_label", "marker_set")]
    setnames(x, by_col, "stratum_value")
    x[, stratum := by_col][]
}
STRATA <- rbindlist(lapply(c("arch_level", "n_traits", "ispleiotropy", "architecture"),
                           strat_tbl), fill = TRUE)
setcolorder(STRATA, c("stratum", "stratum_value", "method_label", "marker_set"))
fwrite(STRATA, file.path(OUTDIR, "strata_tau.tsv"), sep = "\t")

# Continuous simulation properties: correlate per-seed performance with each, per method and
# panel. Spearman, because n = 11-14 seeds and nothing here is expected to be linear.
PER_SEED <- LAND[, .(tau_med = median(tau, na.rm = TRUE)),
                 by = .(seed, method_label, marker_set, meanFst, r2_pc1_temp, final_LA,
                        n_causal_maf01, n_snps, k_best)]
cont_rows <- list()
for (v in c("meanFst", "r2_pc1_temp", "final_LA", "n_causal_maf01", "n_snps", "k_best")) {
    r <- PER_SEED[, {
        ok <- is.finite(get(v)) & is.finite(tau_med)
        if (sum(ok) >= 5L) .(property = v, n_seeds = sum(ok),
                             rho = round(cor(get(v)[ok], tau_med[ok], method = "spearman"), 3))
        else .(property = v, n_seeds = sum(ok), rho = NA_real_)
    }, by = .(method_label, marker_set)]
    cont_rows[[length(cont_rows) + 1L]] <- r
}
CONT <- rbindlist(cont_rows)
fwrite(CONT, file.path(OUTDIR, "property_correlations.tsv"), sep = "\t")

message("\n=== tau by architecture (rows) x method, GEA-best panel ===")
print(dcast(STRATA[stratum == "arch_level" & marker_set == "gea_best"],
            stratum_value ~ method_label, value.var = "tau_med"))
message("\n=== tau by architecture, QTN-oracle panel ===")
print(dcast(STRATA[stratum == "arch_level" & marker_set == "adaptive"],
            stratum_value ~ method_label, value.var = "tau_med"))
message("\n=== Spearman(seed property, tau) -- negative rho = property makes offset BETTER ===")
print(dcast(CONT[marker_set == "gea_best"], property ~ method_label, value.var = "rho"))

# ---- headline table ----------------------------------------------------------
best_lo <- TL[marker_set == "adaptive" & final_LA < 0.42, .(v = median(tau_med)), by = method_label][
    order(v)][1]
best_hi <- TL[marker_set == "adaptive" & final_LA >= 0.42, .(v = median(tau_med)), by = method_label][
    order(v)][1]
OURS <- c(
    sprintf("%.3f (%s)", best_hi$v, best_hi$method_label),
    sprintf("%.3f (%s)", best_lo$v, best_lo$method_label),
    sprintf("%.1f%% / %.1f%%", MS[comparison == "adaptive > neutral", median_gain_pct],
                               MS[comparison == "adaptive > all", median_gain_pct]),
    sprintf("%.1f%%", MS[comparison == "adaptive > neutral", pct_first_wins]),
    sprintf("%.1f%%", MS[comparison == "adaptive > all",     pct_first_wins]),
    sprintf("%.1f%%", MS[comparison == "all > neutral",      pct_first_wins]))
HEAD <- copy(PUB)[, ours := OURS]
fwrite(HEAD, file.path(OUTDIR, "headline_vs_published.tsv"), sep = "\t")

message("=== CLINE-GO vs Lind & Lotterhos 2025 ===")
print(HEAD[, .(quantity, published, ours)])
message("\n=== marker-set comparison, per method ===")
print(MSm[, .(method_label, comparison, n_models, pct_first_wins = round(pct_first_wins, 1),
              median_gain_pct = round(median_gain_pct, 1))])
message("\n=== median tau by method and marker set (landscape gardens) ===")
print(dcast(LAND[, .(tau = median(tau, na.rm = TRUE)), by = .(method_label, marker_set)],
            marker_set ~ method_label, value.var = "tau"))
message("\nWrote tables and figures to ", OUTDIR)
