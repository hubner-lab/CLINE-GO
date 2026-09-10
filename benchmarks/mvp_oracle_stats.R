#!/usr/bin/env Rscript
# =============================================================================
# mvp_oracle_stats.R -- every number the "reaches the oracle" figure plots.
#
# WHY THIS EXISTS. The old headline plotted ABSOLUTE accuracy (-Kendall tau), and
# ~0.62 of that scale is a floor every marker panel gets for free: random and
# neutral markers score 0.62-0.63 because genomic offset tracks isolation-by-
# environment genome-wide. That floor is not an artefact -- it is Lind &
# Lotterhos 2025's own Q3 result -- but it leaves every bar between 0.62 and 0.75
# and nothing separable. Worse, the ranking it produced put the 2/3 panel ABOVE
# the true causal loci, which is not the claim this paper needs.
#
# So the causal-loci panel becomes the REFERENCE AT ZERO and every panel is
# reported as its paired difference from it:
#
#     delta = accuracy(panel) - accuracy(true QTNs),  paired within replicate
#
# Success = REACHING the line in every genetic architecture. Not exceeding it.
#
# WHAT IS AND IS NOT CLAIMED (read before touching the tests).
# A one-sided test of "panel >= oracle" against a point null is a test that the
# panel is significantly BETTER than the oracle -- exactly the claim being
# retired. So:
#   * 2/3 panel  -> DESCRIPTIVE only (median delta, IQR, win%). No p-value in its
#                   favour. A real equivalence claim needs a measured margin
#                   delta, and none exists (rand_best2/rand_best3 were built but
#                   never run through the garden arm).
#   * panels that fall below -> INFERENTIAL. One-sided paired Wilcoxon that the
#                   panel is WORSE than the oracle. That is an ordinary
#                   superiority test of the oracle and is legitimate.
# The paper's point is the second bullet: single GEA methods demonstrably fall
# below the causal loci in at least one architecture; the 2/3 panel does not.
#
# Reads only tables already on disk. Fits nothing that needs the pipeline.
#
# Usage:
#   Rscript /pipeline/benchmarks/mvp_oracle_stats.R
#   OFFSET_DIR=offset11 FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_oracle ...
# =============================================================================

suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_oracle"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rd <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) stop("MISSING: ", f)
    fread(f, colClasses = c("seed" = "character"))
}
emit <- function(dt, stem) {
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    message(sprintf("  wrote %s.tsv (%d rows)", stem, nrow(dt)))
    invisible(dt)
}

# RDA-corrected is pathological (2.65% of its gardens predict BACKWARDS, and on
# 5 replicates it returns byte-identical accuracy for every marker panel). It is
# excluded from every aggregate and appears only in the method panel.
DROP_METHOD <- "RDA-corrected"
WORKING     <- c("GFoffset", "LFMM2offset", "RDA-uncorrected")
ORACLE      <- "causal loci"

ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus ships the typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")

PANELS <- data.table(
    marker_set = c("gea_best", "gea_strict", "gea_union", "gea_rda_only",
                   "gea_emmax_only", "gea_lfmm_only", "all", "random_matched",
                   "neutral", "adaptive"),
    label      = c("2/3 methods", "3/3 methods", "1/3 methods", "RDA",
                   "EMMAX", "LFMM", "all loci", "random, size-matched",
                   "neutral", ORACLE),
    role       = c("combined", "combined", "combined", "single",
                   "single", "single", "reference", "reference",
                   "reference", "oracle"))
PANEL_ORDER <- PANELS$label[c(1, 2, 3, 4, 5, 6, 7, 8, 9)]   # oracle is the line, not a bar

# ---------------------------------------------------------------- shared data
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- seeds[arm == "primary"]
PRIM[, arch_lab := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]
stopifnot(nrow(PRIM) == 90L, !any(is.na(PRIM$arch_lab)))
message(sprintf("primary replicates: %d (%s)", nrow(PRIM),
                paste(sprintf("%s=%d", levels(PRIM$arch_lab), table(PRIM$arch_lab)),
                      collapse = ", ")))

acc <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
acc[, accuracy := -tau]                      # NEVER abs(): an anti-predicting model
                                             # must stay negative, not become accurate
acc <- merge(acc, PRIM[, .(seed, arch_lab)], by = "seed")           # drops the 2 controls
acc <- merge(acc, PANELS[, .(marker_set, label, role)], by = "marker_set")
message(sprintf("accuracy rows: %d (%d seeds, %d offset methods, %d panels)",
                nrow(acc), uniqueN(acc$seed), uniqueN(acc$method_label),
                uniqueN(acc$label)))

# =============================================================================
# delta = panel - oracle, paired within replicate
# =============================================================================
# Aggregation, stated once and used everywhere downstream: median over gardens
# within (seed, method) -- already done upstream in phase1_seed_medians_solo --
# then median over the three WORKING offset methods within seed. A replicate
# carrying more methods must not outweigh one carrying fewer.
seed_acc <- acc[method_label %in% WORKING,
                .(accuracy = median(accuracy), n_methods = .N),
                by = .(seed, arch_lab, label, role)]

wide <- dcast(seed_acc, seed + arch_lab ~ label, value.var = "accuracy")
stopifnot(ORACLE %in% names(wide))

delta_long <- rbindlist(lapply(PANEL_ORDER, function(p) {
    if (!p %in% names(wide)) { message("  panel absent, skipped: ", p); return(NULL) }
    d <- wide[[p]] - wide[[ORACLE]]
    data.table(seed = wide$seed, arch_lab = wide$arch_lab, label = p, delta = d)
}))[is.finite(delta)]
delta_long <- merge(delta_long, PANELS[, .(label, role)], by = "label")
delta_long[, label := factor(label, levels = PANEL_ORDER)]
emit(delta_long, "delta_per_seed")

# One-sided paired Wilcoxon, H1: the panel is WORSE than the oracle.
# Interpreted only where median(delta) < 0 -- see the header.
p_below <- function(d) {
    d <- d[is.finite(d)]
    if (length(d) < 6L || all(d == 0)) return(NA_real_)
    tryCatch(stats::wilcox.test(d, mu = 0, alternative = "less")$p.value,
             error = function(e) NA_real_)
}
p_two <- function(d) {
    d <- d[is.finite(d)]
    if (length(d) < 6L || all(d == 0)) return(NA_real_)
    tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA_real_)
}
summarise_delta <- function(dt, by) dt[, .(
    n              = .N,
    median_delta   = median(delta),
    q25            = quantile(delta, .25),
    q75            = quantile(delta, .75),
    win_pct        = 100 * mean(delta > 0),
    p_below_oracle = p_below(delta),
    p_two_sided    = p_two(delta)), by = by]

by_arch <- summarise_delta(delta_long, c("label", "role", "arch_lab"))
pooled  <- summarise_delta(delta_long, c("label", "role"))[, arch_lab := "pooled"]
BY_ARCH <- rbind(by_arch, pooled, use.names = TRUE)
setorder(BY_ARCH, label, arch_lab)
emit(BY_ARCH, "delta_by_arch")

message("\n=== delta from the causal-loci reference (median), by architecture ===")
print(dcast(BY_ARCH, label ~ arch_lab, value.var = "median_delta"))
message("\n=== p(panel is WORSE than the oracle) -- read only where delta < 0 ===")
print(dcast(BY_ARCH, label ~ arch_lab, value.var = "p_below_oracle"))

# Every cell must be complete. A short cell means seeds were silently dropped
# for a missing panel, which would make the strata non-comparable.
short <- BY_ARCH[arch_lab != "pooled" & n < 30]
if (nrow(short)) {
    message("\n!! INCOMPLETE CELLS (expected n=30 per architecture):")
    print(short[, .(label, arch_lab, n)])
}

# =============================================================================
# worst stratum -- the robustness claim, in delta units
# =============================================================================
WORST <- by_arch[, .(worst_arch     = arch_lab[which.min(median_delta)],
                     worst_delta    = min(median_delta),
                     best_delta     = max(median_delta),
                     spread_delta   = max(median_delta) - min(median_delta),
                     n_strata_below = sum(median_delta < 0)),
                 by = .(label, role)]
setorder(WORST, -worst_delta)
emit(WORST, "worst_stratum")
message("\n=== worst architecture, in delta-from-oracle units ===")
print(WORST)

# =============================================================================
# per offset method -- for the method panel (RDA-corrected INCLUDED here)
# =============================================================================
wide_m <- dcast(acc, seed + arch_lab + method_label ~ label, value.var = "accuracy")
delta_m <- rbindlist(lapply(PANEL_ORDER, function(p) {
    if (!p %in% names(wide_m)) return(NULL)
    data.table(seed = wide_m$seed, arch_lab = wide_m$arch_lab,
               method_label = wide_m$method_label,
               label = p, delta = wide_m[[p]] - wide_m[[ORACLE]])
}))[is.finite(delta)]
BY_METHOD <- delta_m[, .(n = .N, median_delta = median(delta),
                         q25 = quantile(delta, .25), q75 = quantile(delta, .75)),
                     by = .(label, method_label, arch_lab)]
BY_METHOD_POOL <- delta_m[, .(n = .N, median_delta = median(delta),
                              q25 = quantile(delta, .25), q75 = quantile(delta, .75)),
                          by = .(label, method_label)][, arch_lab := "pooled"]
emit(rbind(BY_METHOD, BY_METHOD_POOL, use.names = TRUE), "delta_by_method")
message("\n=== delta from oracle, per offset method (pooled over architecture) ===")
print(dcast(BY_METHOD_POOL, label ~ method_label, value.var = "median_delta"))

# How OFTEN each panel falls below the oracle, per offset method. The collapsed
# version (delta_by_arch) takes the median across the three working methods
# within a replicate first, which answers "does the panel work" but hides
# whether one method is carrying it. This keeps the method axis open.
BELOW_METHOD <- delta_m[method_label %in% WORKING,
                        .(n = .N, below_pct = 100 * mean(delta < 0)),
                        by = .(label, method_label, arch_lab)]
BELOW_METHOD_POOL <- delta_m[method_label %in% WORKING,
                             .(n = .N, below_pct = 100 * mean(delta < 0)),
                             by = .(label, method_label)][, arch_lab := "pooled"]
emit(rbind(BELOW_METHOD, BELOW_METHOD_POOL, use.names = TRUE), "below_by_method")
message("\n=== % of replicates below the oracle, per offset method (pooled) ===")
print(dcast(BELOW_METHOD_POOL, label ~ method_label, value.var = "below_pct"))

# Absolute accuracy per method too -- the method panel needs to show that
# RDA-corrected is bad in level, not merely relative to the oracle.
ABS_METHOD <- acc[, .(n = .N, accuracy = median(accuracy),
                      q25 = quantile(accuracy, .25), q75 = quantile(accuracy, .75)),
                  by = .(label, method_label, arch_lab)]
emit(ABS_METHOD, "accuracy_by_method")

message("\nAll tables written to ", OUT)
