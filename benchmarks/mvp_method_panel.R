#!/usr/bin/env Rscript
# =============================================================================
# mvp_method_panel.R -- panel D candidates: how the four genomic-offset methods
# compare, and whether the marker-set conclusion depends on which one is used.
#
# WHY THIS PANEL EXISTS. Panels B and C take the median across offset methods
# within each replicate, so method performance is averaged away and a reader
# cannot tell whether the marker-set result is carried by one method. It is not:
# `2/3 methods` ranks first under every working method. That robustness check is
# the main point here; the method ranking itself is the second.
#
# Five encodings of the same numbers, to be chosen on how they read:
#   D1 lines        marker set on x, one line per method   (interaction/parallelism)
#   D2 heatmap      method x marker set matrix             (completeness)
#   D3 forest       paired effect of 2/3 vs each reference (effect size + CI)
#   D4 grouped bars classic position encoding
#   D5 boxplots     per-replicate distribution, faceted    (spread, not just median)
#
# Method colours deliberately avoid the green/amber families used for marker sets
# elsewhere, since methods are a different kind of entity. Red stays reserved for
# the broken method.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
METH_COL <- c("GFoffset"        = MINOU[["teal"]],
              "LFMM2offset"     = MINOU[["navy"]],
              "RDA-uncorrected" = MINOU[["grey"]],
              "RDA-corrected"   = MINOU[["red"]])   # reserved: the dropped method
METHS <- names(METH_COL)

SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")

man <- fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character"))[arm=="primary"]
S   <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S   <- S[seed %in% man$seed & marker_set %in% names(SETLAB)]
S[, `:=`(accuracy = -tau,
         set = factor(SETLAB[marker_set], levels = SETLAB),
         method = factor(method_label, levels = METHS))]

# x order: by accuracy pooled over the three WORKING methods, so the broken one
# does not decide the layout.
ordset <- S[method != "RDA-corrected", .(a = median(accuracy)), by = set][order(-a)]$set
S[, set := factor(as.character(set), levels = as.character(ordset))]
MED <- S[, .(accuracy = median(accuracy), n = .N), by = .(set, method)]
fwrite(MED, file.path(OUT, "D_method_by_set.tsv"), sep = "\t")
print(dcast(MED, set ~ method, value.var = "accuracy"))

base_t <- function(p, tag) p + theme_clinego() +
    theme(plot.tag = element_text(face="bold", size=13), panel.border = element_blank(),
          legend.position = "bottom", strip.text = element_text(size=8.5)) +
    labs(tag = tag)

# ---------------- D1: lines ------------------------------------------------
p1 <- base_t(ggplot(MED, aes(set, accuracy, colour = method, group = method)) +
    geom_line(linewidth = 0.95) + geom_point(size = 2.3) +
    scale_colour_manual(values = METH_COL, name = NULL) +
    labs(x = NULL, y = "Accuracy") +
    guides(colour = guide_legend(nrow = 1)), "D1") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))
clinego_save_both(file.path(OUT, "D1_method_lines"), p1, w = 8.4, h = 4.8)

# ---------------- D2: heatmap ----------------------------------------------
p2 <- base_t(ggplot(MED, aes(set, method, fill = accuracy)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.3f", accuracy)), size = 2.6, colour = "white") +
    scale_fill_gradientn(colours = c(MINOU[["red"]], "grey92", MINOU[["teal"]]),
                         name = "accuracy") +
    scale_y_discrete(limits = rev(METHS)) +
    labs(x = NULL, y = NULL), "D2") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
          legend.position = "right")
clinego_save_both(file.path(OUT, "D2_method_heatmap"), p2, w = 9.4, h = 3.6)

# ---------------- D3: paired-difference forest ------------------------------
W <- dcast(S, seed + method ~ set, value.var = "accuracy")
REF <- "2/3 methods"
fo <- rbindlist(lapply(setdiff(levels(S$set), REF), function(o)
    rbindlist(lapply(METHS, function(m) {
        d <- W[method == m][[o]]; r <- W[method == m][[REF]]
        ok <- is.finite(d) & is.finite(r); if (sum(ok) < 10) return(NULL)
        ht <- suppressWarnings(wilcox.test(r[ok], d[ok], paired = TRUE, conf.int = TRUE))
        data.table(other = o, method = m, est = unname(ht$estimate),
                   lo = ht$conf.int[1], hi = ht$conf.int[2], p = ht$p.value)
    }))))
fo[, `:=`(other = factor(other, levels = rev(setdiff(levels(S$set), REF))),
          method = factor(method, levels = METHS))]
fwrite(fo, file.path(OUT, "D_paired_effects.tsv"), sep = "\t")
p3 <- base_t(ggplot(fo, aes(est, other, colour = method)) +
    geom_vline(xintercept = 0, colour = CLINEGO_COL$fg, linewidth = 0.35) +
    geom_linerange(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = 0.7), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.7), size = 1.9) +
    scale_colour_manual(values = METH_COL, name = NULL) +
    labs(x = "accuracy gain of 2/3 methods over the set on the left", y = NULL) +
    guides(colour = guide_legend(nrow = 1)), "D3")
clinego_save_both(file.path(OUT, "D3_method_forest"), p3, w = 8.6, h = 4.8)

# ---------------- D4: grouped bars ------------------------------------------
p4 <- base_t(ggplot(MED, aes(set, accuracy, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.74) +
    scale_fill_manual(values = METH_COL, name = NULL) +
    coord_cartesian(ylim = c(0.45, NA)) +
    labs(x = NULL, y = "Accuracy") +
    guides(fill = guide_legend(nrow = 1)), "D4") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))
clinego_save_both(file.path(OUT, "D4_method_bars"), p4, w = 9.4, h = 4.8)

# ---------------- D5: per-replicate distributions ---------------------------
p5 <- base_t(ggplot(S, aes(set, accuracy, fill = method)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3, alpha = 0.9) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_manual(values = METH_COL, guide = "none") +
    labs(x = NULL, y = "Accuracy"), "D5") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5))
clinego_save_both(file.path(OUT, "D5_method_boxplots"), p5, w = 12.5, h = 4.4)

message("  OK D1-D5")
