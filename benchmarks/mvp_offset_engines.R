#!/usr/bin/env Rscript
# =============================================================================
# mvp_offset_engines.R -- RDA offset with vs without structure correction, measured.
#
# WHY THIS EXISTS. Every aggregate in mvp_oracle_stats.R / mvp_block_report.R drops
# RDA-corrected and medians only GFoffset, LFMM2offset and RDA-uncorrected. The
# reason given in their headers ("2.65 % of its gardens predict BACKWARDS, and on 5
# replicates it returns byte-identical accuracy for every marker panel") was measured
# on an older arm, and the Phase 4a gate (2026-10-01) kept BOTH RDA engines to "see
# which is better". So the exclusion is re-measured here on the arm being reported,
# and the comparison is made directly instead of inherited:
#
#   backward_gardens.tsv   per garden_type x engine x panel: share of gardens with
#                          tau > 0 (offset and fitness move together = the model
#                          predicts backwards), and with tau NA
#   degenerate_seeds.tsv   per (seed, engine): number of distinct panel accuracies;
#                          a replicate whose every panel scores the same value has an
#                          engine that ignores the markers it is given
#   rda_engine_paired.tsv  per panel x architecture (+ pooled): accuracy(corrected) -
#                          accuracy(uncorrected), paired within replicate; median, IQR,
#                          % replicates where corrected is better, two-sided paired
#                          Wilcoxon. Accuracy = -tau, never abs(tau).
#
# Accuracy per (seed, panel, engine) is phase1_seed_medians_solo.tsv: median tau over
# the 100 LANDSCAPE gardens (mvp_panel_tables.R). The backward share is computed per
# garden from garden_performance.tsv, for both garden types.
#
# OFFSET_DIR and FIG_OUT have no defaults on purpose: older arms (offset11, offset12)
# still exist on disk, so a forgotten variable would give a complete, wrong-arm result.
#
# Reads only tables already on disk. Fits nothing.
#
# Usage:
#   OFFSET_DIR=offset13 FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_offset13 \
#   MVP_ARM=primary_ssclines MVP_ADDED=<5 block tags> MVP_N_EXPECT=600 MVP_N_PER_ARCH=200 \
#     Rscript /pipeline/benchmarks/mvp_offset_engines.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "")
OUT  <- Sys.getenv("FIG_OUT", "")
if (!nzchar(OFF) || !nzchar(OUT)) stop("set OFFSET_DIR and FIG_OUT explicitly (no defaults)")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

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

CORR   <- "RDA-corrected"
UNCORR <- "RDA-uncorrected"
ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus ships the typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")

seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- mvp_prim(seeds)
PRIM[, arch_lab := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]
stopifnot(nrow(PRIM) == mvp_n_expect(), !any(is.na(PRIM$arch_lab)))
message(sprintf("arm %s, OFFSET_DIR=%s: %d replicates", mvp_arm_label(), OFF, nrow(PRIM)))

# ------------------------------------------------------------ backward gardens
GP <- rd(EVAL, OFF, "garden_performance.tsv")
GP <- GP[control == FALSE & seed %in% PRIM$seed]
stopifnot(uniqueN(GP$seed) == nrow(PRIM), all(c(CORR, UNCORR) %in% GP$method_label))
BACK <- GP[, .(n_gardens = .N, n_seeds = uniqueN(seed),
               backward_pct = 100 * mean(tau > 0, na.rm = TRUE),
               na_pct = 100 * mean(is.na(tau))),
           by = .(garden_type, method_label, marker_set)]
BACK_ALL <- GP[, .(n_gardens = .N, n_seeds = uniqueN(seed),
                   backward_pct = 100 * mean(tau > 0, na.rm = TRUE),
                   na_pct = 100 * mean(is.na(tau))),
               by = .(garden_type, method_label)][, marker_set := "all panels"]
BACK <- rbind(BACK, BACK_ALL, use.names = TRUE)
setorder(BACK, garden_type, marker_set, method_label)
emit(BACK, "backward_gardens")
message("\n=== % gardens predicting backwards (tau > 0), all panels pooled ===")
print(dcast(BACK_ALL, method_label ~ garden_type, value.var = "backward_pct"))
rm(GP); invisible(gc())

# ------------------------------------------------------------ accuracy per seed
acc <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
acc <- merge(acc, PRIM[, .(seed, arch_lab)], by = "seed")
acc[, accuracy := -tau]                      # NEVER abs(): backwards must stay negative

# A replicate whose every panel scores the same accuracy under an engine: that engine
# is not responding to the markers. Counted over panels present for that seed.
DEG <- acc[, .(n_panels = .N, n_distinct = uniqueN(accuracy)),
           by = .(seed, arch_lab, method_label)]
DEG[, degenerate := n_panels >= 2L & n_distinct == 1L]
emit(DEG, "degenerate_seeds")
message("\n=== replicates where every panel gets one identical accuracy ===")
print(DEG[, .(n_seeds = .N, degenerate = sum(degenerate)), by = method_label])

# ------------------------------------------------------------ paired RDA engines
w <- dcast(acc[method_label %in% c(CORR, UNCORR)], seed + arch_lab + marker_set ~ method_label,
           value.var = "accuracy")
w <- w[is.finite(get(CORR)) & is.finite(get(UNCORR))]
w[, diff := get(CORR) - get(UNCORR)]
p_two <- function(d) {
    d <- d[is.finite(d)]
    if (length(d) < 6L || all(d == 0)) return(NA_real_)
    tryCatch(stats::wilcox.test(d, mu = 0)$p.value, error = function(e) NA_real_)
}
summ <- function(dt, by) dt[, .(
    n                  = .N,
    acc_corrected      = median(get(CORR)),
    acc_uncorrected    = median(get(UNCORR)),
    median_diff        = median(diff),
    q25                = quantile(diff, .25),
    q75                = quantile(diff, .75),
    corrected_better_pct = 100 * mean(diff > 0),
    p_two_sided        = p_two(diff)), by = by]
PAIR <- rbind(summ(w, c("marker_set", "arch_lab")),
              summ(w, "marker_set")[, arch_lab := "pooled"], use.names = TRUE)
setorder(PAIR, marker_set, arch_lab)
emit(PAIR, "rda_engine_paired")
message("\n=== accuracy(corrected) - accuracy(uncorrected), median, pooled ===")
print(PAIR[arch_lab == "pooled", .(marker_set, n, acc_corrected, acc_uncorrected,
                                   median_diff, corrected_better_pct, p_two_sided)])

message("\nAll tables written to ", OUT)
