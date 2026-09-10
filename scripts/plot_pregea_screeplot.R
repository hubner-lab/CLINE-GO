#!/usr/bin/env Rscript
# plot_pregea_screeplot.R — LD-pruned genotype PCA screeplot for preGEA.
#
# Reuses the LEA .eigenvalues file already produced by prestructure mode's
# `pca()` call (W['pca_eigenvalues']) — no PCA is re-run here. Adds a
# broken-stick null model (the standard "how many axes carry real signal"
# heuristic) and marks the swept LFMM-K band + the EMMAX #PC ceiling so the
# preGEA Shiny tab can show, at a glance, whether the ladders it swept are
# even in the plausible structure-axis range.
#
# The count of above-broken-stick axes is an UPPER BOUND on the number of
# real structure axes, not a target for #PCs/#K — see docs/rda_research.md
# C2 (EMMAX block) and the preGEA Shiny note (pregea_tw_upper_bound):
# "#PCs = K_best from sNMF" is not supported by any source. That caveat is
# rendered in Shiny (Code Structure Rule 8), not baked into this plot.

library(ggplot2)
library(data.table)
library(qs)

source("/pipeline/scripts/R/utils/theme_clinego.R")

args <- commandArgs(trailingOnly = TRUE)
################################################################################
EIGENVALUES <- args[1]                       # LEA .eigenvalues, one value/line, no header
K_BEST      <- as.numeric(args[2])
K_RANGE_STR <- args[3]                       # comma-separated swept LFMM K values, e.g. "2,3,4,5"
N_PCS_MAX   <- as.numeric(args[4])           # shared EMMAX/RDA #PC ladder ceiling
OUT_PNG     <- args[5]
OUT_TSV     <- args[6]
INTER_DIR   <- args[7]
################################################################################

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)
dir.create(INTER_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

k_range <- as.integer(strsplit(K_RANGE_STR, ",")[[1]])

message("INFO: Reading eigenvalues from ", EIGENVALUES)
eig <- as.numeric(readLines(EIGENVALUES))
n <- length(eig)

# Broken-stick null model: expected proportion of variance for rank k out of
# n axes = (1/n) * sum_{i=k}^{n} (1/i). An axis is "above broken-stick" when
# its observed proportion exceeds this null expectation — the standard
# significance heuristic for a PCA screeplot (Legendre & Legendre).
broken_stick <- function(n) {
    vapply(seq_len(n), function(k) sum(1 / (k:n)) / n, numeric(1))
}

var_explained <- eig / sum(eig)
bs <- broken_stick(n)

scree_dt <- data.table(
    pc            = seq_len(n),
    eigenvalue    = eig,
    var_explained = var_explained,
    cum_var       = cumsum(var_explained),
    broken_stick  = bs,
    above_bs      = var_explained > bs
)
n_above_bs <- sum(scree_dt$above_bs)
message("INFO: ", n_above_bs, " axis/axes above the broken-stick null (upper bound, not a target)")

fwrite(scree_dt, OUT_TSV, sep = "\t", quote = FALSE)

# Plot: only structural elements (Code Structure Rule 8) — the K range and
# #PC ceiling are geom_vline/geom_hline reference lines (data-derived,
# allowed); the "upper bound not target" caveat text lives in the Shiny note.
plot_dt <- scree_dt[seq_len(min(n, max(N_PCS_MAX, max(k_range), 20) + 5))]
g <- ggplot(plot_dt, aes(x = pc, y = var_explained)) +
    geom_col(aes(fill = above_bs)) +
    geom_line(aes(y = broken_stick), color = CLINEGO_THRESHOLD, linewidth = 0.6) +
    geom_point(aes(y = broken_stick), color = CLINEGO_THRESHOLD, size = 1.2) +
    geom_vline(xintercept = range(k_range), linetype = "dashed", color = CLINEGO_COL$secondary) +
    geom_hline(yintercept = 0, color = "transparent") +
    # CLINEGO_NEUTRAL is an alias of CLINEGO_RETAINED (both #0B775E) — using it
    # here made above/below-broken-stick bars render identically. CLINEGO_THRESHOLD
    # (dark plum) actually distinguishes "below the null" from "above it".
    scale_fill_manual(values = c(`TRUE` = CLINEGO_RETAINED, `FALSE` = CLINEGO_THRESHOLD), guide = "none") +
    labs(x = "PC", y = "Proportion of variance explained",
        title = "LD-pruned genotype PCA screeplot",
        subtitle = sprintf("%d axis/axes above broken-stick null (n=%d total) | K_best=%s | swept K range dashed",
                           n_above_bs, n, K_BEST)) +
    theme_clinego()

ggsave(OUT_PNG, g, width = 8, height = 5, dpi = 300)
ggsave(OUT_SVG, g, width = 8, height = 5, device = svglite::svglite, bg = "transparent")
qsave(g, OUT_QS)

message("INFO: preGEA screeplot complete: ", OUT_PNG)
