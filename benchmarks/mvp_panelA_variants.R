#!/usr/bin/env Rscript
# =============================================================================
# mvp_panelA_variants.R -- alternative normalisations for panel A, plus the
# metric-comparison figure. Additive only: it never touches the count panel
# already written by mvp_main_figure.R, so both can be compared side by side.
#
# WHY. Replicates differ 6.2x in SNP count (4,975-30,650) and 137x in causal
# loci (4-548), so a median of ABSOLUTE counts mixes replicates of very
# different size. The fix is not one obvious normalisation -- different ones
# crown different rules, which is itself the result this script measures:
#
#   volume metrics (counts, per-class recall)      -> 1/3 methods wins
#   ratio metrics (% causal, enrichment, discrim.) -> EMMAX wins
#   offset accuracy (what actually matters)        -> 2/3 methods wins
#
# Outputs:
#   A2_recovery_rates      per-class recovery rate, proper per-replicate
#                          denominators, log axis (scale-free companion to A1)
#   A_counts_and_rates     A1 (existing counts) + A2 stacked
#   M_metric_comparison    every rule under all six metrics, winner marked
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(cowplot)})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
           sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")
shades <- function(base, n) {
    if (n == 1L) return(base)
    c(base, grDevices::colorRampPalette(c(base, "white"))(6)[2:(1 + (n - 1) %/% 2 + (n - 1) %% 2)],
      grDevices::colorRampPalette(c(base, "black"))(6)[2:(1 + (n - 1) %/% 2)])[seq_len(n)]
}
green3 <- shades(MINOU[["sage"]], 3); amber3 <- shades(MINOU[["amber"]], 3)
CLASS_COL <- c("causal" = green3[3], "linked" = MINOU[["sage"]], "neutral" = MINOU[["red"]])
RULE_COL  <- c("2/3 methods" = green3[1], "3/3 methods" = green3[2], "1/3 methods" = green3[3],
               "RDA" = amber3[1], "EMMAX" = amber3[2], "LFMM" = amber3[3])

ARCH <- c(oliogenic = "oligogenic", `mod-polygenic` = "moderately polygenic",
          `highly-polygenic` = "highly polygenic")
LAB  <- c(best = "2/3 methods", intersect3 = "3/3 methods", union = "1/3 methods",
          solo_rda = "RDA", solo_emmax = "EMMAX", solo_lfmm = "LFMM")

man <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
man <- man[arm == "primary", .(seed, arch = factor(ARCH[arch_level], levels = ARCH))]
D <- fread(file.path(EVAL, OFF, "panel_pr_recomputed.tsv"), colClasses = c(seed = "character"))
D[, background := n - n_causal - n_linked]
D <- merge(D, man, by = "seed")

# Per-replicate denominators: the `all` panel IS the whole genome for that replicate.
den <- D[set == "all", .(seed, C_tot = n_causal, L_tot = n_linked,
                         B_tot = background, N_tot = n)]
M <- merge(D[set %in% names(LAB)], den, by = "seed")
M[, rule := factor(LAB[set], levels = names(RULE_COL))]

# Six metrics, one row per replicate x rule.
M[, `:=`(
    rate_causal  = fifelse(C_tot > 0, n_causal   / C_tot, NA_real_),
    rate_linked  = fifelse(L_tot > 0, n_linked   / L_tot, NA_real_),
    rate_neutral = fifelse(B_tot > 0, background / B_tot, NA_real_),
    usable       = n_causal + n_linked,
    pct_causal   = fifelse(n > 0, n_causal / n, NA_real_))]
M[, enrichment := fifelse(n > 0 & C_tot > 0, (n_causal / n) / (C_tot / N_tot), NA_real_)]
M[, discrimination := fifelse(rate_neutral > 0, rate_causal / rate_neutral, NA_real_)]

# =============================== A2: recovery rates ==========================
# Each class divided by how much of that class the replicate actually contained,
# so a 5,000-SNP replicate and a 30,000-SNP one are directly comparable.
a2 <- melt(M, id.vars = c("rule", "arch"),
           measure.vars = c("rate_causal", "rate_linked", "rate_neutral"),
           variable.name = "class", value.name = "rate")
a2[, class := factor(class, levels = c("rate_causal", "rate_linked", "rate_neutral"),
                     labels = names(CLASS_COL))]
a2s <- a2[!is.na(rate) & rate > 0, .(rate = median(rate), n_seeds = .N), by = .(rule, arch, class)]
# Order identical to A1: by usable markers pooled over architectures.
ord <- M[, .(u = median(usable)), by = rule][order(u)]$rule
a2s[, rule := factor(as.character(rule), levels = as.character(ord))]
fwrite(a2s, file.path(OUT, "A2_recovery_rates.tsv"), sep = "\t")

pA2 <- ggplot(a2s, aes(rule, rate, fill = class)) +
    geom_col(width = 0.68, position = position_dodge(width = 0.75)) +
    coord_flip() +
    facet_wrap(~ arch, nrow = 1) +
    scale_y_log10(labels = function(x) paste0(signif(100 * x, 2), "%")) +
    scale_fill_manual(values = CLASS_COL, name = NULL) +
    labs(tag = "A2", x = NULL, y = "markers of that class recovered (log scale)") +
    theme_clinego() +
    theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 13),
          strip.text = element_text(size = 8.5), axis.text.y = element_text(size = 8))
clinego_save_both(file.path(OUT, "A2_recovery_rates"), pA2, w = 11, h = 4.4)
message("  OK A2_recovery_rates")

# ============================ A1 + A2 stacked ================================
# A1 rebuilt here identically to mvp_main_figure.R so the pair can be shown as one
# object; the standalone A written by that script is left untouched.
a1s <- M[, .(causal = median(as.numeric(n_causal)), linked = median(as.numeric(n_linked)),
             neutral = median(as.numeric(background)), total = median(as.numeric(n)),
             usable = median(as.numeric(usable))), by = .(arch, rule)]
a1s[, rule := factor(as.character(rule), levels = as.character(ord))]
a1l <- melt(a1s, id.vars = c("arch", "rule", "total", "usable"),
            measure.vars = c("causal", "linked", "neutral"),
            variable.name = "class", value.name = "markers")
a1l[, class := factor(class, levels = names(CLASS_COL))]
fwrite(a1s, file.path(OUT, "A1_counts.tsv"), sep = "\t")

pA1 <- ggplot(a1l, aes(rule, markers, fill = class)) +
    geom_col(width = 0.7) +
    geom_text(data = a1s, aes(rule, total, label = as.integer(usable)),
              inherit.aes = FALSE, hjust = -0.25, size = 2.7, colour = CLINEGO_COL$fg) +
    coord_flip() + facet_wrap(~ arch, nrow = 1) +
    scale_fill_manual(values = CLASS_COL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(tag = "A1", x = NULL, y = "markers selected (median)") +
    theme_clinego() +
    theme(legend.position = "none", plot.tag = element_text(face = "bold", size = 13),
          strip.text = element_text(size = 8.5), axis.text.y = element_text(size = 8))
clinego_save_both(file.path(OUT, "A_counts_and_rates"),
                plot_grid(pA1, pA2, ncol = 1, rel_heights = c(1, 1.12)), w = 11, h = 8.6)
message("  OK A_counts_and_rates")

# ======================== M: the metric comparison ===========================
# Every rule under all six detection metrics, plus the accuracy it achieves.
# The point is that the winner MOVES: no detection metric picks the rule that
# predicts best.
acc <- fread(file.path(EVAL, OFF, "phase1_seed_medians_solo.tsv"), colClasses = c(seed = "character"))
acc <- acc[method_label != "RDA-corrected"]
ACCLAB <- c(gea_best = "2/3 methods", gea_strict = "3/3 methods", gea_union = "1/3 methods",
            gea_rda_only = "RDA", gea_emmax_only = "EMMAX", gea_lfmm_only = "LFMM")
acc <- acc[marker_set %in% names(ACCLAB)]
acc[, rule := ACCLAB[marker_set]]
acc <- merge(acc, man, by = "seed")
accs <- acc[, .(v = median(-tau)), by = .(rule, arch)]

mets <- rbindlist(list(
    M[, .(v = median(usable)),                        by = .(rule, arch)][, metric := "usable markers (count)"],
    M[, .(v = median(rate_causal,  na.rm = TRUE)),    by = .(rule, arch)][, metric := "causal recovery rate"],
    M[, .(v = median(pct_causal,   na.rm = TRUE)),    by = .(rule, arch)][, metric := "% of calls that are causal"],
    M[, .(v = median(enrichment,   na.rm = TRUE)),    by = .(rule, arch)][, metric := "fold enrichment for causal"],
    M[, .(v = median(discrimination, na.rm = TRUE)),  by = .(rule, arch)][, metric := "causal rate ÷ neutral rate"],
    accs[, .(rule, arch, v)][, metric := "OFFSET ACCURACY"]))
mets[, metric := factor(metric, levels = c(
    "usable markers (count)", "causal recovery rate", "% of calls that are causal",
    "fold enrichment for causal", "causal rate ÷ neutral rate", "OFFSET ACCURACY"))]
mets[, rule := factor(as.character(rule), levels = names(RULE_COL))]
mets[, best := v == max(v), by = .(metric, arch)]
fwrite(dcast(mets, rule + arch ~ metric, value.var = "v"),
       file.path(OUT, "M_metric_comparison.tsv"), sep = "\t")
print(dcast(mets[best == TRUE], metric ~ arch, value.var = "rule"))

pM <- ggplot(mets, aes(rule, v, fill = rule, alpha = best)) +
    geom_col(width = 0.7) +
    geom_point(data = mets[best == TRUE], aes(y = v), shape = 8, size = 1.9,
               colour = CLINEGO_COL$fg, show.legend = FALSE, inherit.aes = TRUE) +
    coord_flip() +
    facet_grid(arch ~ metric, scales = "free_x") +
    scale_alpha_manual(values = c(`FALSE` = 0.42, `TRUE` = 1), guide = "none") +
    scale_fill_manual(values = RULE_COL, guide = "none") +
    labs(tag = "M", x = NULL, y = NULL) +
    theme_clinego() +
    theme(plot.tag = element_text(face = "bold", size = 13),
          strip.text = element_text(size = 7.5), axis.text.x = element_text(size = 6.5),
          axis.text.y = element_text(size = 7.5), panel.border = element_blank())
clinego_save_both(file.path(OUT, "M_metric_comparison"), pM, w = 14, h = 6)
message("  OK M_metric_comparison")
message("\nwritten to ", OUT)
