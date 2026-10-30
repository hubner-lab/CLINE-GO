#!/usr/bin/env Rscript
# cross_entropy_table.R — sNMF's cross-entropy for every K and run as a table,
# PreStructure/tables/cross_entropy_K{start}-{end}.tsv (f59f46). plot_cross_entropy.R draws
# each K's best run and wrote no numbers, so choosing sNMF.k_best from pipeline output meant
# reading the .snmfProject by hand.
#
# One row per (K, run): LEA's cross.entropy(). best_run marks the run with the lowest
# cross-entropy at that K: the run extract_clusters.R takes the K's Q-matrix from, and the
# point plot_cross_entropy.R draws.
#
#   Rscript cross_entropy_table.R <project.snmfProject> <k_start> <k_end> <out.tsv>

suppressPackageStartupMessages({
    library(LEA)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L)
    stop("usage: cross_entropy_table.R <project.snmfProject> <k_start> <k_end> <out.tsv> (got ",
         length(args), " arguments)")
################################################################################
SNMF_PROJECT <- args[1]
K_START      <- as.integer(args[2])
K_END        <- as.integer(args[3])
OUT_TSV      <- args[4]
################################################################################

snmf_res <- load.snmfProject(SNMF_PROJECT)

out <- rbindlist(lapply(K_START:K_END, function(K) {
    ce <- as.numeric(cross.entropy(snmf_res, K = K))
    data.table(K = K, run = seq_along(ce), cross_entropy = ce,
               best_run = seq_along(ce) == which.min(ce))
}))
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE)

best <- out[best_run == TRUE]
message("INFO: cross-entropy of the best run per K: ",
        paste0("K", best$K, " = ", signif(best$cross_entropy, 5), collapse = ", "))
