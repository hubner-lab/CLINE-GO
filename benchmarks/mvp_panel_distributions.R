#!/usr/bin/env Rscript
# mvp_panel_distributions.R -- the marker-panel comparison shown as DISTRIBUTIONS
# across the 30 replicates, not as single medians.
#
# Why this exists: every replicate is an independent simulation with its own
# demography, architecture and SNP count, so panel size, precision, recall and
# offset accuracy all vary between them. A median alone hides whether a gap is
# consistent across replicates or driven by a handful. Each replicate is one
# point here, and the panel-to-panel tests are PAIRED WITHIN replicate --
# comparing panels within the same simulation, which is what the design supports.
#
# Panel names spelled out: the three GEA methods are LFMM, RDA and EMMAX, so
# "union" = LFMM u RDA u EMMAX, ">=2 agree" = at least two of those three, and
# "all 3 agree" = the intersection of all three.

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
# Which scored offset root to read. Defaults to offset09, the 32-replicate run, so an
# unparameterised call reproduces the frozen figure set; a new cohort passes OFFSET_DIR.
OFF <- Sys.getenv("OFFSET_DIR", "offset09")
OUT <- Sys.getenv("FIG_OUT", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/figures"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Recomputed directly from each panel file vs the replicate's truth table, keyed
# on chr:pos. NOT snp_sets_summary.tsv -- that file is overwritten by whichever
# --sets run wrote it last, so it covered 32 replicates for the solo panels but
# only 18 for the rest, which would compare panels on different replicate sets.
SS <- fread(file.path(PIPELINE_ROOT, file.path("benchmarks/mvp_eval", OFF, "panel_pr_recomputed.tsv")),
            colClasses = c(seed = "character"))
SS[, `:=`(n_snps = n, written = TRUE)]
S  <- fread(file.path(PIPELINE_ROOT, file.path("benchmarks/mvp_eval", OFF, "phase1_seed_medians_solo.tsv")),
            colClasses = c(seed = "character"))
LAB <- c(truth = "true QTNs\n(oracle)", intersect3 = "LFMM & RDA & EMMAX\n(all 3 agree)",
         best = "any 2 of\nLFMM/RDA/EMMAX", solo_emmax = "EMMAX\nonly",
         solo_rda = "RDA\nonly", solo_lfmm = "LFMM\nonly",
         union = "LFMM or RDA or EMMAX\n(union)",
         all = "ALL SNPs\n(no selection)", neutral_all = "neutral SNPs\n(reference)",
         rand_best1 = "random,\nsize-matched")
RL <- c(adaptive = "true QTNs\n(oracle)", gea_strict = "LFMM & RDA & EMMAX\n(all 3 agree)",
        gea_best = "any 2 of\nLFMM/RDA/EMMAX", gea_emmax_only = "EMMAX\nonly",
        gea_rda_only = "RDA\nonly", gea_lfmm_only = "LFMM\nonly",
        gea_union = "LFMM or RDA or EMMAX\n(union)",
        all = "ALL SNPs\n(no selection)", neutral = "neutral SNPs\n(reference)",
        random_matched = "random,\nsize-matched")
SS[, panel := LAB[set]]; SS <- SS[!is.na(panel) & written == TRUE]
S[,  panel := RL[marker_set]]; S <- S[!is.na(panel)]


A <- S[, .(accuracy = -median(tau)), by = .(seed, panel)]
D <- merge(SS[, .(seed, panel, n_snps, precision, recall)], A, by = c("seed", "panel"))

ord <- D[, .(m = median(accuracy)), by = panel][order(-m)]$panel
D[, panel := factor(panel, levels = ord)]

# ---- table: median [IQR] and full range, per replicate ----------------------
q <- function(x) sprintf("%.3f [%.3f-%.3f]", median(x), quantile(x, .25), quantile(x, .75))
tab <- D[, .(n_replicates = .N,
             markers = sprintf("%d [%d-%d]", as.integer(median(n_snps)),
                               as.integer(min(n_snps)), as.integer(max(n_snps))),
             precision = q(precision), recall = q(recall), accuracy = q(accuracy)),
         by = panel][order(match(panel, ord))]
print(tab)
fwrite(tab, file.path(OUT, "panel_distribution_table.tsv"), sep = "\t")

# ---- paired tests: every panel vs the best-performing one -------------------
best <- as.character(ord[1])
res <- rbindlist(lapply(setdiff(as.character(ord), best), function(p) {
    m <- merge(D[panel == best, .(seed, a = accuracy)], D[panel == p, .(seed, b = accuracy)], by = "seed")
    w <- suppressWarnings(wilcox.test(m$a, m$b, paired = TRUE))
    data.table(panel = p, n = nrow(m), median_diff = median(m$a - m$b),
               wins_for_best = sum(m$a > m$b), p_value = w$p.value)
}))
res[, p_holm := p.adjust(p_value, "holm")]
cat(sprintf("\n=== paired Wilcoxon vs '%s' (n=30 replicates) ===\n", gsub("\n", " ", best)))
print(res)
fwrite(res, file.path(OUT, "panel_paired_tests.tsv"), sep = "\t")

# ---- figure: distributions, one point per replicate ------------------------
mk <- function(var, title, sub, ylab, logy = FALSE) {
    p <- ggplot(D, aes(panel, get(var))) +
        geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey95",
                     colour = CLINEGO_COL$secondary, linewidth = 0.4) +
        geom_jitter(width = 0.13, height = 0, size = 1.5, alpha = 0.65, colour = CLINEGO_RETAINED) +
        labs(title = title, subtitle = sub, x = NULL, y = ylab) +
        theme_clinego() + theme(axis.text.x = element_text(size = 7.5))
    if (logy) p <- p + scale_y_log10()
    p
}
ggsave(file.path(OUT, "dist_accuracy.png"),
       mk("accuracy", "Prediction accuracy per replicate",
          "Each dot = one of 30 simulated replicates. Higher = better. Box = median and IQR.",
          "prediction accuracy (-tau)"), width = 10, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "dist_markers.png"),
       mk("n_snps", "Marker-set size varies a lot between replicates",
          "Log scale. Each dot = one replicate.", "markers in set", logy = TRUE),
       width = 10, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "dist_recall.png"),
       mk("recall", "Recall: share of true causal loci captured",
          "Each dot = one replicate.", "recall"), width = 10, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "dist_precision.png"),
       mk("precision", "Precision: share of the set that is truly causal",
          "Each dot = one replicate.", "precision"), width = 10, height = 5.5, dpi = 300, bg = "white")
message("\nwritten to ", OUT)
