#!/usr/bin/env Rscript
# plot_pregea_ladder.R — assemble one preGEA ladder (LFMM-K or EMMAX-#PC) into
# a wide ladder TSV + 4 plot families. Serves BOTH engines — the wide
# p-value table shape and the rung_stats long-table shape are identical
# across engines (Code Structure Rule: avoid redundancy).
#
# Rule 6: the QQ grid uses geom_scattermore (plot_pop_diff.R's idiom), NOT
# geom_point — N scales with SNP count even on the LD-pruned set. Do NOT
# copy plot_manhattan.R's geom_point QQ (an existing Rule 6 violation there).
# Rule 8: no annotate() commentary — only data-derived geom_hline/geom_vline
# reference lines. Rule 9: theme_clinego() + ADAPT_* constants throughout.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scattermore)
})

source("/pipeline/scripts/R/utils/theme_clinego.R")
source("/pipeline/scripts/R/utils/io_pvalues.R")

args <- commandArgs(trailingOnly = TRUE)
################################################################################
ENGINE          <- args[1]                   # "lfmm" | "emmax"
RUNG_PARAM      <- args[2]                   # "K" | "n_pcs"
RUNG_STATS_STR  <- args[3]                   # quoted, space-separated rung_stats.tsv paths
PVALUES_STR     <- args[4]                   # quoted, space-separated pvalues.tsv paths, SAME order
RUNG_VALUES_STR <- args[5]                   # comma-separated rung values, SAME order
DEFLATION_FLOOR <- args[6]                   # numeric string, or "NULL" (lfmm has none)
FDR             <- as.numeric(args[7])
OUT_LADDER_TSV  <- args[8]
PLOT_DIR        <- args[9]
INTER_DIR       <- args[10]
################################################################################

dir.create(dirname(OUT_LADDER_TSV), recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INTER_DIR, recursive = TRUE, showWarnings = FALSE)

deflation_floor <- if (identical(DEFLATION_FLOOR, "NULL")) NA_real_ else as.numeric(DEFLATION_FLOOR)

# Display labels — user-facing text only; ENGINE/RUNG_PARAM stay as the raw
# wildcard values everywhere else (paths, TSV columns) so nothing downstream
# (pregea_recommend.R, write_summary.R) needs to change.
ENGINE_LABEL <- c(lfmm = "LFMM", emmax = "EMMAX / GAPIT")[[ENGINE]]
RUNG_LABEL   <- c(K = "K", n_pcs = "#PCs")[[RUNG_PARAM]]

stats_paths  <- strsplit(trimws(RUNG_STATS_STR), "\\s+")[[1]]
pvals_paths  <- strsplit(trimws(PVALUES_STR), "\\s+")[[1]]
rung_values  <- strsplit(RUNG_VALUES_STR, ",")[[1]]
stopifnot(length(stats_paths) == length(rung_values), length(pvals_paths) == length(rung_values))

message("INFO: Assembling ", ENGINE, " ladder over ", RUNG_PARAM, " = ", RUNG_VALUES_STR)

################################################################################
# 1. Long stats table across all rungs -> wide ladder TSV
################################################################################
long_stats <- rbindlist(lapply(stats_paths, fread, sep = "\t", header = TRUE,
                               colClasses = "character"))
# Numeric metrics widened to columns; hist_shape/hist_spike0 stay as text.
ladder_wide <- dcast(long_stats, engine + rung_param + rung_value + trait ~ metric, value.var = "value")
num_cols <- intersect(c("n_tests", "lambda_gc", "lambda_gc_applied", "hits_bonf", "hits_qval",
                        "hist_flatness_ks", "frac_p_gt_half"), names(ladder_wide))
ladder_wide[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
ladder_wide[, rung_value_num := suppressWarnings(as.numeric(rung_value))]
setorder(ladder_wide, rung_value_num, trait)
fwrite(ladder_wide, OUT_LADDER_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote ladder table: ", OUT_LADDER_TSV, " (", nrow(ladder_wide), " rows)")

################################################################################
# 2. Load all rung p-value tables into one long data.table for the grids
################################################################################
pval_long_list <- vector("list", length(pvals_paths))
for (i in seq_along(pvals_paths)) {
    dt <- read_pvalues_tsv(pvals_paths[i])
    fixed_cols <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
    trait_cols <- setdiff(names(dt), fixed_cols)
    molten <- melt(dt, id.vars = "SNPID", measure.vars = trait_cols,
                   variable.name = "trait", value.name = "pvalue")
    molten[, rung_value := rung_values[i]]
    pval_long_list[[i]] <- molten
}
pval_long <- rbindlist(pval_long_list)
pval_long[, rung_value_num := suppressWarnings(as.numeric(rung_value))]
pval_long <- pval_long[!is.na(pvalue)]

rung_lab <- function(v) paste0(RUNG_LABEL, "=", v)
pval_long[, rung_label := factor(rung_lab(rung_value_num), levels = rung_lab(sort(unique(rung_value_num))))]

################################################################################
# 3. p-value histogram grid — THE selection signal (C.0/C.1): read shape, not lambda
################################################################################
g_hist <- ggplot(pval_long, aes(x = pvalue)) +
    geom_histogram(breaks = seq(0, 1, length.out = 41), fill = CLINEGO_NEUTRAL, color = NA) +
    facet_grid(rung_label ~ trait) +
    labs(x = "p-value", y = "Count",
        title = paste0(ENGINE_LABEL, " candidate p-value histogram grid"),
        subtitle = "Read the SHAPE per panel: flat + spike near 0 is the target; not lambda") +
    theme_clinego_grid()
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_pvalue_histogram_grid.png")), g_hist,
      width = 2.2 * uniqueN(pval_long$trait) + 2, height = 1.6 * uniqueN(pval_long$rung_value) + 1.5,
      dpi = 200, limitsize = FALSE)
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_pvalue_histogram_grid.svg")), g_hist,
      width = 2.2 * uniqueN(pval_long$trait) + 2, height = 1.6 * uniqueN(pval_long$rung_value) + 1.5,
      device = svglite::svglite, bg = "transparent", limitsize = FALSE)

################################################################################
# 4. QQ grid — geom_scattermore (Rule 6), not geom_point
################################################################################
pval_long[order(rung_value_num, trait, pvalue),
         expected := -log10(ppoints(.N)), by = .(rung_value_num, trait)]
pval_long[, observed := -log10(pvalue)]
g_qq <- ggplot(pval_long, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, color = CLINEGO_THRESHOLD, linetype = "dashed") +
    geom_scattermore(color = CLINEGO_NEUTRAL, pointsize = 3, pixels = c(600, 600)) +
    facet_grid(rung_label ~ trait) +
    labs(x = expression(-log[10](expected)), y = expression(-log[10](observed)),
        title = paste0(ENGINE_LABEL, " QQ grid")) +
    theme_clinego_grid()
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_qq_grid.png")), g_qq,
      width = 2.2 * uniqueN(pval_long$trait) + 2, height = 1.6 * uniqueN(pval_long$rung_value) + 1.5,
      dpi = 200, limitsize = FALSE)
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_qq_grid.svg")), g_qq,
      width = 2.2 * uniqueN(pval_long$trait) + 2, height = 1.6 * uniqueN(pval_long$rung_value) + 1.5,
      device = svglite::svglite, bg = "transparent", limitsize = FALSE)

################################################################################
# 5. lambda-vs-rung (pooled only — a single curve, not per-trait clutter)
################################################################################
pooled <- ladder_wide[trait == "__pooled__"]
g_lambda <- ggplot(pooled, aes(x = rung_value_num, y = lambda_gc)) +
    geom_hline(yintercept = 1, color = CLINEGO_THRESHOLD, linetype = "dashed")
if (!is.na(deflation_floor)) {
    g_lambda <- g_lambda + geom_hline(yintercept = deflation_floor, color = CLINEGO_REMOVED, linetype = "dotted")
}
g_lambda <- g_lambda +
    geom_line(color = CLINEGO_NEUTRAL) + geom_point(color = CLINEGO_NEUTRAL, size = 2) +
    labs(x = RUNG_LABEL, y = expression(lambda[GC]),
        title = paste0(ENGINE_LABEL, " genomic inflation factor vs ", RUNG_LABEL)) +
    theme_clinego()
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_lambda_vs_", RUNG_PARAM, ".png")), g_lambda,
      width = 7, height = 4.5, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_lambda_vs_", RUNG_PARAM, ".svg")), g_lambda,
      width = 7, height = 4.5, device = svglite::svglite, bg = "transparent")

################################################################################
# 6. hit-count vs K/#PCs at fixed FDR (pooled only, per C.2:315/321) — an
# OVER-CORRECTION detector, not a power measure: read the TREND, not the
# height. A steep DROP as K/#PCs rises means the model is absorbing real
# signal (over-correction onset); a stable plateau means the correction is
# doing its job. On small test datasets (e.g. SIMDATA) this line is often
# flat at 1 hit across the whole range — the metric needs WGS-scale SNP
# density to separate values; that is not a bug (see fct_help_notes.R).
################################################################################
g_hits <- ggplot(pooled, aes(x = rung_value_num, y = hits_qval)) +
    geom_line(color = CLINEGO_NEUTRAL) + geom_point(color = CLINEGO_NEUTRAL, size = 2) +
    labs(x = RUNG_LABEL, y = sprintf("Significant SNPs (FDR %.2g)", FDR),
        title = paste0(ENGINE_LABEL, " significant-hit count vs ", RUNG_LABEL)) +
    theme_clinego()
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_hits_vs_", RUNG_PARAM, ".png")), g_hits,
      width = 7, height = 4.5, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(ENGINE, "_hits_vs_", RUNG_PARAM, ".svg")), g_hits,
      width = 7, height = 4.5, device = svglite::svglite, bg = "transparent")

message("INFO: preGEA ", ENGINE, " ladder plots complete in ", PLOT_DIR)
