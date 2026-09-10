#!/usr/bin/env Rscript
# =============================================================================
# mvp_oracle_figure.R -- the main simulation figure, rebuilt around the oracle.
#
# The claim being made: the pipeline's >=2-of-3 GEA marker panel PREDICTS AS WELL
# AS the true causal loci, in every genetic architecture, while single GEA
# methods do not. So the causal-loci panel is the reference line at ZERO and
# every other panel is its paired difference from it. Reaching the line is the
# success; exceeding it is not the message and carries no p-value here.
#
# Reads only the tables written by benchmarks/mvp_oracle_stats.R (plus the panel
# composition table already on disk). Fits nothing.
#
# Two encodings are built for the two panels whose reading is a judgement call:
#   B1 grouped bars  vs  B2 per-replicate boxplots
#   D1 method lines  vs  D2 method ranking on the recommended panel
#
# Usage:
#   Rscript /pipeline/benchmarks/mvp_oracle_figure.R
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
    library(ggrepel)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_oracle"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

# ---- palette --------------------------------------------------------------
# Same semantics as figures_main so a reader who has seen one figure can read the
# other: green = combining rules, amber = single methods, grey = uninformative
# references, navy = the causal-loci reference, red = the dropped method.
MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
           sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")
shades <- function(base, n) {
    if (n == 1L) return(base)
    c(base,
      grDevices::colorRampPalette(c(base, "white"))(6)[2:(1 + (n - 1) %/% 2 + (n - 1) %% 2)],
      grDevices::colorRampPalette(c(base, "black"))(6)[2:(1 + (n - 1) %/% 2)])[seq_len(n)]
}
green3 <- shades(MINOU[["sage"]],  3)
amber3 <- shades(MINOU[["amber"]], 3)
grey3  <- shades(MINOU[["grey"]],  3)

PANEL_COL <- c("2/3 methods" = green3[1], "3/3 methods" = green3[2],
               "1/3 methods" = green3[3], "RDA" = amber3[1], "EMMAX" = amber3[2],
               "LFMM" = amber3[3], "all loci" = grey3[3],
               "random, size-matched" = grey3[2], "neutral" = MINOU[["grey"]])
METHOD_COL <- c("GFoffset" = MINOU[["teal"]], "LFMM2offset" = MINOU[["amber"]],
                "RDA-uncorrected" = MINOU[["sage"]], "RDA-corrected" = MINOU[["red"]])
CLASS_COL  <- c(causal = shades(MINOU[["sage"]], 3)[3], linked = MINOU[["sage"]],
                neutral = MINOU[["red"]])
ORACLE_COL <- MINOU[["navy"]]

ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
PANEL_ORDER <- c("2/3 methods", "3/3 methods", "1/3 methods", "RDA", "EMMAX",
                 "LFMM", "all loci", "random, size-matched", "neutral")

rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("MISSING: ", f); fread(f) }
save_fig <- function(p, stem, w, h) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h)
    message(sprintf("  OK %s", stem))
}
base_t <- function(p, tag) p + theme_clinego() +
    theme(plot.tag = element_text(face = "bold", size = 13),
          strip.text = element_text(size = 8.5), legend.position = "bottom") +
    labs(tag = tag)

DELTA   <- rd(OUT, "delta_by_arch.tsv")
PERSEED <- rd(OUT, "delta_per_seed.tsv")
WORST   <- rd(OUT, "worst_stratum.tsv")
DMETH   <- rd(OUT, "delta_by_method.tsv")

# =============================================================================
# A -- what each selection rule returns (unchanged from the old figure)
# =============================================================================
# A EXPLAINS B and the caption must say so: under highly polygenic architecture
# `3/3 methods` collapses to ~21 markers and `EMMAX` to ~4, which is exactly why
# those two panels fall off the oracle line in B under that architecture.
message("\n=== A: composition of what each rule returns ===")
PR_KEY <- c(best = "2/3 methods", intersect3 = "3/3 methods", union = "1/3 methods",
            solo_rda = "RDA", solo_emmax = "EMMAX", solo_lfmm = "LFMM",
            truth = "causal loci", all = "all loci", neutral_all = "neutral",
            rand_best1 = "random, size-matched")
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- seeds[arm == "primary"]
PRIM[, arch_lab := factor(arch_level,
                          levels = c("oliogenic", "mod-polygenic", "highly-polygenic"),
                          labels = ARCH_LABELS)]
comp <- rd(EVAL, OFF, "panel_pr_recomputed.tsv")
comp <- merge(comp, PRIM[, .(seed, arch_lab)], by = "seed")
comp[, label := PR_KEY[set]]
comp <- comp[!is.na(label)]
comp[, `:=`(usable = n_causal + n_linked, background = n - n_causal - n_linked)]

OBTAINABLE <- c("2/3 methods", "3/3 methods", "1/3 methods", "RDA", "EMMAX", "LFMM")
a_dt  <- comp[label %in% OBTAINABLE]
a_ord <- a_dt[, .(u = median(usable)), by = label][order(u)]$label
a_dt[, label := factor(label, levels = a_ord)]
a_src <- a_dt[, .(causal = median(n_causal), linked = median(n_linked),
                  background = median(background), total = median(as.numeric(n)),
                  usable = median(usable), n_seeds = .N), by = .(arch_lab, label)]
a_long <- melt(a_src, id.vars = c("arch_lab", "label", "total", "usable", "n_seeds"),
               measure.vars = c("causal", "linked", "background"),
               variable.name = "class", value.name = "markers")
a_long[, class := factor(class, levels = c("causal", "linked", "background"),
                         labels = names(CLASS_COL))]
fwrite(a_long, file.path(OUT, "A_rule_composition.tsv"), sep = "\t")

nc <- comp[label == "causal loci", .(nc = round(median(n_causal))), by = arch_lab]
arch_facet <- setNames(sprintf("%s\n(~%d causal loci)", nc$arch_lab, nc$nc), nc$arch_lab)

pA <- base_t(ggplot(a_long, aes(label, markers, fill = class)) +
    geom_col(width = 0.7) +
    geom_text(data = a_src, aes(label, total, label = as.integer(usable)),
              inherit.aes = FALSE, hjust = -0.25, size = 2.7, colour = CLINEGO_COL$fg) +
    coord_flip() +
    facet_wrap(~ arch_lab, nrow = 1, labeller = labeller(arch_lab = arch_facet)) +
    scale_fill_manual(values = CLASS_COL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "markers selected (median)"), "A") +
    theme(axis.text.y = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 1))
save_fig(pA, "A_rule_composition", w = 11, h = 4.6)

# =============================================================================
# B -- HEADLINE: distance from the causal-loci reference
# =============================================================================
# Zero is the oracle. A panel ON the line has done the job. The two encodings
# differ only in what they admit: B1 shows the medians cleanly, B2 shows the
# 90-replicate spread behind them and is honest about the overlap.
message("\n=== B: distance from the oracle ===")
ORDB <- WORST[order(-worst_delta)]$label                # best worst-case first
b_dt <- DELTA[arch_lab != "pooled"]
b_dt[, `:=`(label = factor(as.character(label), levels = ORDB),
            arch_lab = factor(arch_lab, levels = ARCH_LABELS))]
ps <- copy(PERSEED)
ps[, `:=`(label = factor(as.character(label), levels = ORDB),
          arch_lab = factor(arch_lab, levels = ARCH_LABELS))]

zero_line <- function(p) p +
    geom_hline(yintercept = 0, colour = ORACLE_COL, linewidth = 0.7)

pB1 <- base_t(zero_line(ggplot(b_dt, aes(label, median_delta, fill = arch_lab)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72) +
    geom_linerange(aes(ymin = q25, ymax = q75),
                   position = position_dodge(width = 0.78), linewidth = 0.35,
                   colour = CLINEGO_COL$fg)) +
    scale_fill_manual(values = setNames(c(green3[1], MINOU[["teal"]], MINOU[["navy"]]),
                                        ARCH_LABELS), name = NULL) +
    labs(x = NULL, y = "accuracy relative to the causal loci"), "B") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8)) +
    guides(fill = guide_legend(nrow = 1))
save_fig(pB1, "B1_delta_bars", w = 10.5, h = 5.0)

pB2 <- base_t(zero_line(ggplot(ps, aes(label, delta, fill = label)) +
    geom_boxplot(outlier.size = 0.45, linewidth = 0.3, alpha = 0.92, width = 0.72)) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_manual(values = PANEL_COL, guide = "none") +
    labs(x = NULL, y = "accuracy relative to the causal loci"), "B") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(pB2, "B2_delta_boxplots", w = 12.5, h = 5.0)

# =============================================================================
# C -- robustness: the worst architecture a panel is ever asked to handle
# =============================================================================
message("\n=== C: worst-stratum distance from the oracle ===")
w_dt <- copy(WORST)[, label := factor(as.character(label), levels = rev(ORDB))]
pC <- base_t(ggplot(w_dt, aes(label, worst_delta, fill = label)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, colour = ORACLE_COL, linewidth = 0.7) +
    geom_text(aes(label = worst_arch, y = pmin(worst_delta, 0) - 0.004),
              hjust = 1, size = 2.5, colour = CLINEGO_COL$muted) +
    coord_flip() +
    scale_fill_manual(values = PANEL_COL, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.34, 0.08))) +
    labs(x = NULL, y = "worst-architecture accuracy relative to the causal loci"), "C") +
    theme(axis.text.y = element_text(size = 8))
save_fig(pC, "C_worst_stratum", w = 8.6, h = 4.6)

# =============================================================================
# C2 -- HOW OFTEN a panel falls below the causal loci
# =============================================================================
# The complement of Lind & Lotterhos 2025's Figure 3 statistic ("adaptive
# outperformed 68% of models using neutral markers"). Same number, flipped:
#
#     below_pct = 100 - win_pct = share of replicates where the panel scored
#                                 BELOW the true causal loci
#
# Flipped on purpose. win% puts the recommended panel at 64-70%, which reads as
# "beats the truth two thirds of the time" -- the claim this figure exists to
# retire. Stated as a shortfall rate, low is good, and "better than the causal
# loci" is not representable on the axis at all.
#
# 50% is the reference: a panel indistinguishable from the oracle wins half the
# replicates by chance, so the colour scale is anchored with its midpoint there.
#
# Dropped from BOTH encodings: the pooled column (it averages over the very
# factor the panel is about) and the `random, size-matched` / `neutral` rows
# (known-dead references, and at 96-97% they swallowed the whole red end of the
# scale, flattening the 36-51% contrast that carries the result). `all loci`
# stays -- it is the no-curation default a real user faces, not a control.
message("\n=== C2: how often a panel falls below the causal loci ===")
C2_DROP <- c("random, size-matched", "neutral")
h_dt <- DELTA[arch_lab != "pooled" & !label %in% C2_DROP,
              .(label, arch_lab, n, below_pct = 100 - win_pct)]
h_dt[, arch_lab := factor(arch_lab, levels = ARCH_LABELS)]
hord <- h_dt[, .(m = mean(below_pct)), by = label][order(-m)]$label
h_dt[, label := factor(as.character(label), levels = hord)]
fwrite(h_dt[order(label, arch_lab)], file.path(OUT, "C2_below_oracle.tsv"), sep = "\t")
print(dcast(h_dt, label ~ arch_lab, value.var = "below_pct"))

LEG <- "replicates worse than the causal loci (%)"

pC2 <- base_t(ggplot(h_dt, aes(arch_lab, label, fill = below_pct)) +
    geom_tile(colour = "white", linewidth = 0.9) +
    geom_text(aes(label = sprintf("%.0f%%", below_pct),
                  colour = below_pct < 25 | below_pct > 75),
              size = 3.2, fontface = "bold", show.legend = FALSE) +
    scale_fill_gradientn(
        colours = c("#14614a", MINOU[["sage"]], "grey92", MINOU[["red"]], "#8f2438"),
        values  = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100), labels = c("0", "25", "50", "75", "100"),
        name = LEG) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = CLINEGO_COL$fg)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL), "C") +
    theme(axis.text.x = element_text(size = 8.5),
          axis.text.y = element_text(size = 9),
          panel.border = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          legend.key.width = unit(30, "pt"), legend.key.height = unit(8, "pt"))
save_fig(pC2, "C2_below_oracle_heatmap", w = 7.4, h = 4.6)

# --- same numbers as trends across polygenicity ------------------------------
# The heatmap states each cell; this states the SLOPE. Architecture is ordered
# oligogenic -> highly polygenic, so a flat line means the panel does not care
# how polygenic the trait is, and a rising line means it degrades as the signal
# spreads over more loci. That comparison is the point of the panel.
l_dt <- copy(h_dt)[, label := factor(as.character(label), levels = rev(hord))]
lab_end <- l_dt[arch_lab == ARCH_LABELS[3]]
pC3 <- base_t(ggplot(l_dt, aes(arch_lab, below_pct, colour = label, group = label)) +
    geom_hline(yintercept = 50, linetype = "22", linewidth = 0.4,
               colour = ORACLE_COL) +
    annotate("text", x = 0.62, y = 50, label = "no better or worse\nthan the causal loci",
             hjust = 0, vjust = -0.25, size = 2.7, lineheight = 0.95,
             colour = ORACLE_COL) +
    geom_line(linewidth = 1) + geom_point(size = 2.4) +
    ggrepel::geom_text_repel(data = lab_end, aes(label = label), hjust = 0,
                             nudge_x = 0.12, direction = "y", size = 3,
                             segment.size = 0.3, min.segment.length = 0,
                             show.legend = FALSE) +
    scale_colour_manual(values = PANEL_COL, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.10, 0.42))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(x = NULL, y = LEG), "C") +
    theme(axis.text.x = element_text(size = 9), axis.title.y = element_text(size = 8.5))
save_fig(pC3, "C3_below_oracle_lines", w = 7.6, h = 4.8)

# --- same statistic, method axis kept open -----------------------------------
# C2/C3 take the median across the three working offset methods within each
# replicate before counting, which answers "does the panel work" but cannot say
# whether one offset method is carrying the result. Here every (architecture x
# offset method) cell is its own column, grouped under its architecture.
# RDA-corrected is absent by construction -- see S-RDA.
message("\n=== C4: below-oracle rate, per architecture x offset method ===")
BM <- rd(OUT, "below_by_method.tsv")
BM <- BM[arch_lab != "pooled" & !label %in% C2_DROP]
BM[, `:=`(arch_lab = factor(arch_lab, levels = ARCH_LABELS),
          method_label = factor(method_label,
                                levels = c("GFoffset", "LFMM2offset", "RDA-uncorrected"),
                                labels = c("GF", "LFMM2", "RDA")),
          label = factor(as.character(label), levels = hord))]
fwrite(BM[order(label, arch_lab, method_label)],
       file.path(OUT, "C4_below_oracle_by_method.tsv"), sep = "\t")
print(dcast(BM, label ~ arch_lab + method_label, value.var = "below_pct"))

pC4 <- base_t(ggplot(BM, aes(method_label, label, fill = below_pct)) +
    geom_tile(colour = "white", linewidth = 0.9) +
    geom_text(aes(label = sprintf("%.0f", below_pct),
                  colour = below_pct < 25 | below_pct > 75),
              size = 3.0, fontface = "bold", show.legend = FALSE) +
    facet_grid(. ~ arch_lab, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_gradientn(
        colours = c("#14614a", MINOU[["sage"]], "grey92", MINOU[["red"]], "#8f2438"),
        values  = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100), name = LEG) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = CLINEGO_COL$fg)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL), "C") +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 9),
          panel.border = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          strip.placement = "outside", panel.spacing.x = unit(6, "pt"),
          legend.key.width = unit(30, "pt"), legend.key.height = unit(8, "pt"))
save_fig(pC4, "C4_below_oracle_by_method", w = 9.0, h = 4.6)

# =============================================================================
# D -- the maladaptation step: does the conclusion depend on the offset method?
# =============================================================================
message("\n=== D: offset methods ===")
d_pool <- DMETH[arch_lab == "pooled"]
d_pool[, `:=`(label = factor(as.character(label), levels = ORDB),
              method_label = factor(method_label, levels = names(METHOD_COL)))]

pD1 <- base_t(ggplot(d_pool, aes(label, median_delta,
                                 colour = method_label, group = method_label)) +
    geom_hline(yintercept = 0, colour = ORACLE_COL, linewidth = 0.7) +
    geom_line(linewidth = 0.95) + geom_point(size = 2.2) +
    scale_colour_manual(values = METHOD_COL, name = NULL) +
    labs(x = NULL, y = "accuracy relative to the causal loci"), "D") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8)) +
    guides(colour = guide_legend(nrow = 1))
save_fig(pD1, "D1_method_lines", w = 9.4, h = 4.8)

d_best <- DMETH[label == "2/3 methods" & arch_lab != "pooled"]
d_best[, `:=`(method_label = factor(method_label, levels = names(METHOD_COL)),
              arch_lab = factor(arch_lab, levels = ARCH_LABELS))]
pD2 <- base_t(ggplot(d_best, aes(method_label, median_delta, fill = method_label)) +
    geom_col(width = 0.72) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.35, colour = CLINEGO_COL$fg) +
    geom_hline(yintercept = 0, colour = ORACLE_COL, linewidth = 0.7) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_manual(values = METHOD_COL, guide = "none") +
    labs(x = NULL, y = "accuracy relative to the causal loci"), "D") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))
save_fig(pD2, "D2_method_ranking", w = 9.6, h = 4.6)

message("\nFigures written to ", OUT)
