#!/usr/bin/env Rscript
# =============================================================================
# mvp_method_figure.R -- how the four genomic-offset engines compare ON THE
# RECOMMENDED >=2-of-3 MARKER PANEL, in ABSOLUTE accuracy.
#
# WHY THIS EXISTS. The main figure's panel B scores every engine against ITS OWN
# causal-loci run, so its columns are not comparable to each other. Read as if
# they were, it appears to say gradientForest does badly under high polygenicity
# (53% vs LFMM2's 77%). It does not: GF's causal-loci accuracy there is 0.759,
# the highest oracle value anywhere in that figure, so GF is measured against the
# hardest denominator. This figure removes the denominator entirely and plots the
# accuracies themselves.
#
# Three encodings of the same numbers, to be chosen on how they read:
#   V1 grouped bars      medians + IQR, classic position encoding
#   V2 boxplots          the per-replicate distribution behind those medians
#   V3 paired forest     each engine minus GF, paired within replicate, with CI
#                        -- the direct answer to "is GF worse?"
#
# Usage:  Rscript /pipeline/benchmarks/mvp_method_figure.R
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(ggplot2)})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "method_figure"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
           sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")
# Red is reserved for the engine the pipeline drops (see S-RDA): 2.65% of its
# gardens predict BACKWARDS, and on 5 replicates it returns byte-identical
# accuracy for every marker panel.
METH_COL <- c("gradientForest\noffset" = MINOU[["teal"]],
              "LFMM2\noffset"          = MINOU[["navy"]],
              "RDA\noffset"            = MINOU[["sage"]],
              "RDA\noffset (corrected)"= MINOU[["red"]])
METH_LEV <- names(METH_COL)

ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")

rd <- function(...) { f <- file.path(...); if (!file.exists(f)) stop("MISSING: ", f); fread(f) }
save_fig <- function(p, stem, w, h) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h); message("  OK ", stem)
}
base_t <- function(p, tag) p + theme_clinego() +
    theme(plot.tag = element_text(face = "bold", size = 13),
          strip.text = element_text(size = 8.5), legend.position = "bottom") +
    labs(tag = tag)

seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- mvp_prim(seeds)
PRIM[, arch_lab := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]

S <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
S <- S[marker_set == "gea_best"]
S[, accuracy := -tau]                       # NEVER abs(): see header
S <- merge(S, PRIM[, .(seed, arch_lab)], by = "seed")
S[, method := factor(c(GFoffset = METH_LEV[1], LFMM2offset = METH_LEV[2],
                       `RDA-uncorrected` = METH_LEV[3],
                       `RDA-corrected` = METH_LEV[4])[method_label], levels = METH_LEV)]
stopifnot(!any(is.na(S$method)))
message(sprintf("%d rows, %d seeds, %d engines", nrow(S), uniqueN(S$seed), uniqueN(S$method)))

MED <- S[, .(n = .N, accuracy = median(accuracy),
             q25 = quantile(accuracy, .25), q75 = quantile(accuracy, .75)),
         by = .(arch_lab, method)]
setorder(MED, arch_lab, -accuracy)
fwrite(MED, file.path(OUT, "method_accuracy_by_arch.tsv"), sep = "\t")
print(dcast(MED, method ~ arch_lab, value.var = "accuracy"))

YLAB <- "prediction accuracy (-Kendall's tau)"

# ---------------- V1: grouped bars -----------------------------------------
p1 <- base_t(ggplot(MED, aes(method, accuracy, fill = method)) +
    geom_col(width = 0.72) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.4, colour = CLINEGO_COL$fg) +
    geom_text(aes(label = sprintf("%.3f", accuracy)), vjust = -0.6, size = 2.7,
              colour = CLINEGO_COL$fg) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_manual(values = METH_COL, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = NULL, y = YLAB), "V1") +
    theme(axis.text.x = element_text(size = 7.5, lineheight = 0.95),
          axis.title.y = element_text(size = 8.5))
save_fig(p1, "V1_method_bars", w = 9.6, h = 4.4)

# ---------------- V2: per-replicate distributions ---------------------------
p2 <- base_t(ggplot(S, aes(method, accuracy, fill = method)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3, alpha = 0.92, width = 0.7) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_manual(values = METH_COL, guide = "none") +
    labs(x = NULL, y = YLAB), "V2") +
    theme(axis.text.x = element_text(size = 7.5, lineheight = 0.95),
          axis.title.y = element_text(size = 8.5))
save_fig(p2, "V2_method_boxplots", w = 9.6, h = 4.4)

# ---------------- V3: paired difference vs gradientForest -------------------
# Paired within replicate: the same seed, the same gardens, the same marker
# panel -- only the offset engine differs. Anything peculiar to a replicate
# cancels, which an unpaired comparison of medians cannot claim.
W <- dcast(S, seed + arch_lab ~ method, value.var = "accuracy")
REF <- METH_LEV[1]
fo <- rbindlist(lapply(setdiff(METH_LEV, REF), function(m)
    rbindlist(lapply(ARCH_LABELS, function(a) {
        d <- W[arch_lab == a]; x <- d[[m]]; r <- d[[REF]]
        ok <- is.finite(x) & is.finite(r); if (sum(ok) < 10) return(NULL)
        ht <- suppressWarnings(wilcox.test(x[ok], r[ok], paired = TRUE, conf.int = TRUE))
        data.table(method = m, arch_lab = a, n = sum(ok),
                   est = unname(ht$estimate), lo = ht$conf.int[1], hi = ht$conf.int[2],
                   p = ht$p.value)
    }))))
fo[, `:=`(method = factor(method, levels = rev(METH_LEV)),
          arch_lab = factor(arch_lab, levels = ARCH_LABELS))]
fwrite(fo, file.path(OUT, "method_paired_vs_GF.tsv"), sep = "\t")
print(fo[, .(method, arch_lab, est = round(est, 3), lo = round(lo, 3),
             hi = round(hi, 3), p = signif(p, 2))])

p3 <- base_t(ggplot(fo, aes(est, method, colour = method)) +
    geom_vline(xintercept = 0, colour = MINOU[["teal"]], linewidth = 0.8) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.55) +
    geom_point(size = 2.2) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_colour_manual(values = METH_COL, guide = "none") +
    labs(x = "accuracy relative to gradientForest offset (paired, 95% CI)",
         y = NULL), "V3") +
    theme(axis.text.y = element_text(size = 7.5, lineheight = 0.95),
          axis.title.x = element_text(size = 8.5))
save_fig(p3, "V3_paired_vs_GF", w = 10.0, h = 3.8)

message("\nWrote ", OUT)
