#!/usr/bin/env Rscript
# correlation_table.R — the correlation matrix a correlogram draws, as a long pair table (bffcfc):
#   climate/tables/present/climate_correlation.tsv   beside correlation_heatmap.png
#   Traits/tables/trait_correlation.tsv              beside correlation_heatmap_traits.png
#   Traits/tables/trait_climate_correlation.tsv      beside correlation_heatmap_traits_climate.png
# factor_correlation() computes it, the function plot_correlation_heatmap.R draws from, on the
# same tables, so the numbers are the figure's (the figure rounds them to one digit).
#
# One row per pair of non-constant factors; constant ones are left out, as in the figure:
#   factor_1, factor_2  the pair, in the matrix's order
#   r                   Pearson correlation, pairwise-complete
#   basis               sites: both factors are constant within every site, so the pair was
#                       correlated over one row per site; samples: over the sample rows
#   n_obs               rows (sites or samples) with both values present
#
#   Rscript correlation_table.R <table_a.tsv> <table_b.tsv|NULL> <out.tsv>

suppressPackageStartupMessages(library(data.table))

source("/pipeline/scripts/R/lib/factor_correlation.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
    stop("usage: correlation_table.R <table_a.tsv> <table_b.tsv|NULL> <out.tsv> (got ",
         length(args), " arguments)")
################################################################################
TABLE_A <- args[1]
TABLE_B <- args[2]   # "NULL" for a single block
OUT_TSV <- args[3]
################################################################################

pairs <- correlation_pairs(factor_correlation(TABLE_A, TABLE_B))
fwrite(pairs, OUT_TSV, sep = "\t", quote = FALSE)

message("INFO: wrote ", nrow(pairs), " factor pair(s); |r| >= 0.7 in ",
        sum(abs(pairs$r) >= 0.7, na.rm = TRUE))
