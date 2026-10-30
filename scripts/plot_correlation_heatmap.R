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
source("/pipeline/scripts/R/utils/factor_order.R")
source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
TABLE_A   = args[1]    # required wide table (numeric factor columns + metadata columns)
TABLE_B   = args[2]    # optional second block, "NULL" to skip
OUT_PNG   = args[3]
INTER_DIR = args[4]
TITLE     = if (length(args) >= 5 && nzchar(args[5])) args[5] else "Correlogram of factors"
FACTORS   = if (length(args) >= 6 && !args[6] %in% c("", "all", "NULL")) args[6] else ""
            # comma-separated factors to draw (Climate.predictors): cuts each block that names
            # any of them, factors_to_draw(); "all" or absent draws every factor (c9edc9)
#################################
if (length(args) < 4 || length(args) > 6)
  stop("plot_correlation_heatmap.R takes 4 to 6 arguments, got ", length(args))

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- file.path(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Plot theme
plot_theme <-
  theme_clinego(base_size = 10) +
  theme(plot.title.position = "plot",
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10))

#################################

# The matrix: factor_correlation() (scripts/R/lib/factor_correlation.R) joins the blocks,
# drops constant factors and correlates site-level factors over one row per site, on EVERY
# factor; correlation_table.R writes it as a table (bffcfc). The figure then draws the factors
# the rule names, each table's block in place (traits before the predictors) and numbered names
# in numeric order within it (bio_1, bio_2, ..., bio_10; 9e5dd0's rule), so each pair it shows
# has the table's r.
cor_matrix <- factor_correlation(TABLE_A, TABLE_B)$r
keep <- if (nzchar(FACTORS)) trimws(strsplit(FACTORS, ",")[[1]]) else character(0)
blocks <- lapply(Filter(function(t) !identical(t, "NULL") && nzchar(t), list(TABLE_A, TABLE_B)),
                 function(t) factor_cols(read_factor_block(t)))
draw <- factors_to_draw(colnames(cor_matrix), blocks, keep)
missing <- setdiff(keep, colnames(cor_matrix))
if (length(missing) > 0)
  message("INFO: listed factor(s) not in the matrix (constant, all-NA or absent): ",
          paste(missing, collapse = ", "))
draw <- unlist(lapply(blocks, function(b) {
  x <- intersect(b, draw)
  x[numbered_order(x)]
}), use.names = FALSE)
message("INFO: drawing ", length(draw), " of ", ncol(cor_matrix), " factors: ",
        paste(draw, collapse = ", "))

if (length(draw) >= 2) {
  gHM <- ggcorrplot(cor_matrix[draw, draw],
                    type = "lower",
                    lab = TRUE,
                    digits = 1,
                    lab_size = 4) +
    ggtitle(TITLE) +
    plot_theme
} else {
  # Empty state: one factor has no pair to correlate (CLAUDE.md rule 8's placeholder exception).
  gHM <- ggplot() +
    annotate("text", x = 0, y = 0, size = 4.5, colour = CLINEGO_COL$muted,
             label = if (length(draw) == 1) paste0("One factor (", draw, "): no pair to correlate")
                     else "No factor to correlate") +
    ggtitle(TITLE) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 12, colour = CLINEGO_COL$fg))
}

# Determine size based on number of factors
scale_factor <- 0.7
width <- max(4, length(draw) * scale_factor)
height <- max(4, length(draw) * scale_factor)

# Save
ggsave(OUT_PNG, gHM, width = width, height = height, units = "in")
ggsave(OUT_SVG, gHM, width = width, height = height, units = "in",
       device = svglite::svglite, bg = 'transparent')
qsave(gHM, OUT_QS)
