#!/usr/bin/env Rscript
# tracy_widom_table.R — the Tracy-Widom test of every PCA axis as a table,
# PreStructure/tables/tracy_widom.tsv (6bf499). plot_pca.R runs the same test for its plot but
# wrote no numbers; the table has its own rule because pca_plot's projections feed EMMAX, so a
# new output there would rerun every scan on an existing tree.
#
# LEA's own routine runs on a copy of the eigenvalues file pca_plot wrote: R_tracyWidom, the C
# code tracy.widom() calls for a pcaProject. The values are LEA's, under clearer names:
#   axis               N             the PCA axis (1 = PC1)
#   eigenvalue         eigenvalues
#   tw_statistic       twstats
#   tw_pvalue          pvalues
#   effective_n        effectn       the effective number of markers (Patterson et al. 2006)
#   variance_fraction  percentage    the axis's share of the summed eigenvalues: a fraction
#
#   Rscript tracy_widom_table.R <file.eigenvalues> <out.tsv>

suppressPackageStartupMessages({
    library(LEA)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
    stop("usage: tracy_widom_table.R <file.eigenvalues> <out.tsv> (got ", length(args),
         " arguments)")
################################################################################
EIGENVALUES <- args[1]
OUT_TSV     <- args[2]
################################################################################

# A private copy, so LEA's .tracywidom side file lands in a temp dir and not beside pca_plot's
# outputs. LEA requires the .eigenvalues extension.
work <- tempfile("tracy_widom_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
ev_copy <- file.path(work, "pca.eigenvalues")
tw_file <- file.path(work, "pca.tracywidom")
if (!file.copy(EIGENVALUES, ev_copy))
    stop("tracy_widom_table.R: cannot copy ", EIGENVALUES)

.C("R_tracyWidom", as.character(ev_copy), as.character(tw_file), PACKAGE = "LEA")
tw <- read.table(tw_file, header = TRUE)

lea_cols <- c("N", "eigenvalues", "twstats", "pvalues", "effectn", "percentage")
if (!identical(colnames(tw), lea_cols))
    stop("tracy_widom_table.R: LEA wrote columns ", paste(colnames(tw), collapse = ", "),
         ", expected ", paste(lea_cols, collapse = ", "))

out <- data.table(axis              = tw$N,
                  eigenvalue        = tw$eigenvalues,
                  tw_statistic      = tw$twstats,
                  tw_pvalue         = tw$pvalues,
                  effective_n       = tw$effectn,
                  variance_fraction = tw$percentage)
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE)

message("INFO: Tracy-Widom: ", nrow(out), " axes, ", sum(out$tw_pvalue < 0.05),
        " with p < 0.05 (", paste0("PC", out$axis[out$tw_pvalue < 0.05], collapse = ", "), ")")
