#!/usr/bin/env Rscript
# Custom Manhattan plot for association analysis
# Generates clean, customizable Manhattan plots using ggplot2
# Optionally highlights significant regions with method-specific shapes
#
# DERIVED COPY of scripts/plot_manhattan.R (pipeline, 2026-10-16 state) for the Trifolium paper figures (plan
# ~/.claude/plans/fancy-scribbling-garden.md, Phase 5: LFMM supplement at K 5). Lines marked `# DERIVED` differ; the
# change is entry D10 of manuscript_fig3/pipeline_deviations.md, awaiting the user's approval:
#   D10 several traits (arg 5, draw order) stacked in one figure instead of one PNG per trait; minou colours
#       (chromosomes FIG3_CHR, significant SNPs in the predictor's FIG3_PRED_COL, threshold line MS_INK) instead of
#       get_chr_colors() / get_trait_colors() / red; theme_ms() + ms_save() at a fixed mm size (arg 6 = output stem,
#       arg 8 = WxH mm) instead of theme_manhattan() + ggsave(10 x 4 in); the scattermore raster sized to the panel
#       (add_scatter_layer() hard-codes 3000 x 1200 px); no subtitle (the count goes to the caption), plain y title
#       (plotmath needs the Symbol font, as D7); sig points smaller (2.5 -> 0.9) for the mm size; no QQ, background
#       PNG or coords JSON (Shiny-only outputs); per trait, n_sig is asserted equal to the pipeline sig table (arg 7).
# Run inside the image: Rscript plot_manhattan_ms.R <pvalues.tsv> <adjust> <K> <method> <trait1,trait2,...>
#   <out_stem> <sig_snps.tsv> <W>x<H>

library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(qvalue)
library(scales)
library(scattermore)
library(patchwork)  # DERIVED (D10)

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/manhattan_utils.R")
source("/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R")  # DERIVED (D10)

args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_TABLE    = args[1]   # SNPID chr pos TRAITS...
ADJUST         = args[2]   # e.g., "bonf_0.05" or "qval_0.05"
Kbest          = args[3] %>% as.numeric
METHOD         = args[4]   # EMMAX or LFMM (current plot method)
TRAITS         = str_split(args[5], ',')[[1]]  # DERIVED (D10): several traits, in draw order
PLOT_STEM      = args[6]   # DERIVED (D10): output stem (.svg + .png)
SIG_TABLE      = args[7]   # DERIVED (D10): pipeline sig table, for the n_sig assertion
SIZE_MM        = as.numeric(str_split(args[8], 'x')[[1]])  # DERIVED (D10): width x height in mm
################################

message(paste0('INFO: Manhattan plots for ', METHOD, ' - ', paste(TRAITS, collapse = ', ')))  # DERIVED (D10)
message(paste0('INFO: K = ', Kbest))
message(paste0('INFO: Adjustment: ', ADJUST))

# DERIVED (D10): add_scatter_layer() with the raster sized to the panel (pipeline: pixels = c(3000, 1200) fixed)
add_scatter_layer_ms <- function(data, chr_colors, alpha = 0.5, pointsize, pixels) {
    data <- data[is.finite(data$log10p) & is.finite(data$pos_cum), ]
    if (nrow(data) == 0) {
        return(ggplot2::geom_blank())
    }
    data$point_color <- chr_colors[as.character(data$chr_f)]
    scattermore::geom_scattermore(data = data,
                                  ggplot2::aes(x = pos_cum, y = log10p, color = point_color),
                                  alpha = alpha, pointsize = pointsize,
                                  pixels = pixels, interpolate = FALSE)
}
PANEL_H <- SIZE_MM[2] / length(TRAITS)                                        # DERIVED (D10)
PIXELS  <- round(c(SIZE_MM[1] - 10, PANEL_H - 10) / 25.4 * 300)               # DERIVED (D10): data area (minus ~10 mm of title/axes) at 300 px/in, as the pipeline
POINTSZ <- 6                                                                  # DERIVED (D10): ~0.5 mm at 300 px/in

################################ Main

# Load association data
message('INFO: Loading association data')
snps_assoc <- fread(ASSOC_TABLE, sep = '\t', header = T)
snps_assoc$chr <- as.character(snps_assoc$chr)
n_rows <- nrow(snps_assoc)
message(paste0('INFO: Loaded ', n_rows, ' rows'))
sig_table <- fread(SIG_TABLE, select = c("SNPID", "trait"))                  # DERIVED (D10)

# Parse adjustment parameters
adjustment         <- ADJUST %>% str_split('_') %>% unlist %>% .[1]
original_threshold <- ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric

message(paste0('INFO: Adjustment method: ', adjustment))
message(paste0('INFO: Threshold: ', original_threshold))

panel_counts <- list()                                                        # DERIVED (D10)
panels <- lapply(seq_along(TRAITS), function(i) {                             # DERIVED (D10): one panel per trait
TRAIT <- TRAITS[i]

# Check if trait exists
if (!(TRAIT %in% colnames(snps_assoc))) {
    stop(paste0("Trait '", TRAIT, "' not found in association table"))
}
sig_color <- unname(FIG3_PRED_COL[TRAIT])                                     # DERIVED (D10)
stopifnot(!is.na(sig_color))                                                  # DERIVED (D10)

# Compute threshold — compute_pval_threshold filters NA internally (handles WZA
# windows with no valid SNPs for the current trait).
thresh_result  <- compute_pval_threshold(snps_assoc[[TRAIT]], adjustment, original_threshold)
threshold_ok   <- thresh_result$status == "ok"
pval_threshold <- thresh_result$threshold   # NA when !threshold_ok

message(paste0('INFO: ', TRAIT, ' threshold status=', thresh_result$status,
               ' threshold=', pval_threshold,
               ' n_tested=', thresh_result$n_tested,
               ' n_na_dropped=', thresh_result$n_na_dropped))

threshold_log10 <- if (threshold_ok) -log10(pval_threshold) else NA_real_

# Prepare data for Manhattan plot
plot_data <- snps_assoc %>%
    dplyr::select(SNPID, chr, pos, !!sym(TRAIT)) %>%
    dplyr::rename(pvalue = !!sym(TRAIT))

manhattan_prep <- prepare_manhattan_data(plot_data, "chr", "pos", "pvalue")
plot_df <- manhattan_prep$data
chr_info <- manhattan_prep$chr_info

# Get chromosome colors
n_chr <- length(unique(plot_df$chr_f))
chr_colors <- rep(FIG3_CHR, length.out = n_chr)                              # DERIVED (D10)
names(chr_colors) <- levels(plot_df$chr_f)

# Mark significant SNPs (none if threshold not available)
plot_df <- plot_df %>%
    dplyr::mutate(is_significant = if (threshold_ok) log10p >= threshold_log10 else FALSE)

n_sig <- sum(plot_df$is_significant, na.rm = TRUE)
message(paste0('INFO: Found ', n_sig, ' significant SNPs above threshold'))
stopifnot(n_sig == sig_table[trait == TRAIT, .N])                             # DERIVED (D10)
panel_counts[[TRAIT]] <<- data.table(trait = TRAIT, n_sig = n_sig, threshold = pval_threshold)  # DERIVED (D10)

# y-axis ceiling
y_max <- max(plot_df$log10p, na.rm = TRUE) * 1.1
if (threshold_ok) y_max <- max(y_max, threshold_log10 * 1.2)

################################ Plot 1: Simple Manhattan (static export — threshold line baked in)

p_simple <- ggplot() +
    add_scatter_layer_ms(data = plot_df %>% filter(!is_significant),          # DERIVED (D10)
                      chr_colors = chr_colors,
                      alpha = 0.5, pointsize = POINTSZ, pixels = PIXELS) +    # DERIVED (D10)
    scale_color_identity()

p_simple <- p_simple +
    geom_point(data = plot_df %>% filter(is_significant),
               aes(x = pos_cum, y = log10p),
               color = sig_color, alpha = 0.9, size = 0.9) +                  # DERIVED (D10)
    scale_x_continuous(
        labels = chr_info$chr_f,
        breaks = chr_info$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.05)),
        limits = c(0, y_max)
    ) +
    labs(
        title    = paste0(METHOD, " - ", TRAIT),
        x        = if (i == length(TRAITS)) "Chromosome" else NULL,          # DERIVED (D10)
        y        = "-log10 P"                                                 # DERIVED (D10)
    ) +
    theme_ms() + theme(plot.title = element_text(size = MS_BASE, colour = MS_INK))  # DERIVED (D10)

# Add threshold line to static export only when threshold is valid
if (threshold_ok) {
    p_simple <- p_simple +
        geom_hline(yintercept = threshold_log10, linetype = "dashed",
                   color = MS_INK, linewidth = 0.3)                           # DERIVED (D10)
}
p_simple
})                                                                            # DERIVED (D10)

fig <- wrap_plots(panels, ncol = 1)                                           # DERIVED (D10)
ms_save(PLOT_STEM, fig, SIZE_MM[1], SIZE_MM[2])                               # DERIVED (D10)
fwrite(rbindlist(panel_counts), paste0(PLOT_STEM, "_counts.tsv"), sep = "\t")  # DERIVED (D10)
message(paste0('INFO: Saved ', PLOT_STEM, ' (.svg, .png, _counts.tsv)'))      # DERIVED (D10)

message('INFO: Manhattan plot complete')
