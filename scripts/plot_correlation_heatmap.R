#!/usr/bin/env Rscript
# Correlogram of a wide numeric factor table, optionally joined with a second
# block of factors.
#
# Generalized (Phase 3) from the old climate-only signature, which took a
# climate table as arg 1 and positionally cbind()ed the trait columns of a
# metadata table onto it. Two callers now need different combinations:
#
#   Environmental (mode=climate)  A = climate_present_site.tsv, B = NULL
#   Phenotypic    (mode=traits)   A = metadata.tsv,             B = NULL
#   Phenotypic + climate          A = metadata_climate_valid.tsv, B = climate_present_site.tsv
#
# The positional cbind() it replaces was only correct because the climate rule
# was careful to pass the climate-valid metadata subset (same row count AND
# same row order). That invariant lived in a Snakefile comment, not in code.
# The two blocks are joined on `sample` whenever both carry it, and a row-count
# mismatch without a shared key is a hard error instead of a silent misalignment
# (factor_correlation(), scripts/R/lib/factor_correlation.R, since 2026-10-30).

library(ggplot2)
library(data.table)
library(qs)
library(ggcorrplot)

source("/pipeline/scripts/R/lib/factor_correlation.R")

args = commandArgs(trailingOnly=TRUE)
#################################
TABLE_A   = args[1]    # required wide table (numeric factor columns + metadata columns)
TABLE_B   = args[2]    # optional second block, "NULL" to skip
OUT_PNG   = args[3]
INTER_DIR = args[4]
TITLE     = if (length(args) >= 5 && nzchar(args[5])) args[5] else "Correlogram of factors"
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- file.path(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Plot theme
plot_theme <-
  theme_classic(base_size = 8, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = 'bold', size = 10, angle = -90),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 12))

#################################

# The matrix: factor_correlation() (scripts/R/lib/factor_correlation.R) joins the blocks,
# drops constant factors and correlates site-level factors over one row per site.
# correlation_table.R writes the same matrix as a table (bffcfc).
cor_matrix <- factor_correlation(TABLE_A, TABLE_B)$r

# Build the heatmap using ggcorrplot
gHM <- ggcorrplot(cor_matrix,
                  type = "lower",
                  lab = TRUE,
                  digits = 1,
                  lab_size = 4) +
  ggtitle(TITLE) +
  plot_theme

# Determine size based on number of factors
scale_factor <- 0.7
width <- max(4, ncol(cor_matrix) * scale_factor)
height <- max(4, ncol(cor_matrix) * scale_factor)

# Save
ggsave(OUT_PNG, gHM, width = width, height = height, units = "in")
ggsave(OUT_SVG, gHM, width = width, height = height, units = "in",
       device = svglite::svglite, bg = 'transparent')
qsave(gHM, OUT_QS)
