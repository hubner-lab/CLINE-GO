#!/usr/bin/env Rscript
# assoc_lambda.R — the genomic inflation factor (lambda_GC) of every trait column of one
# GEA/GWAS method's p-value table, written beside it as {method}_lambda_K{k}.tsv (00d7a8).
#
# lambda_GC = the median 1-df chi-square statistic of the trait's non-NA p-values over
# qchisq(0.5, 1) — lambda_gc() in pregea_calibration.R, the formula preGEA uses. `correction`
# says whether the method already rescaled its p-values by an inflation factor, which holds
# lambda near 1 by construction, so a lambda of 1.00 there is not evidence of calibration:
#   genomic_control  LFMM with genomic_control: true (lfmm2.test divides by its GIF)
#   rdadapt_gif      RDA (rdadapt() divides the Mahalanobis distances by their GIF)
#   none             EMMAX, the GAPIT models, LFMM with genomic_control: false
# It is resolved from the method's configs entry by method_pvalue_correction() in
# common.smk and passed in, never inferred here.
#
#   Rscript assoc_lambda.R <pvalues.tsv> <method> <correction> <out.tsv>

suppressPackageStartupMessages(library(data.table))

source("/pipeline/scripts/R/lib/pregea_calibration.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L)
    stop("usage: assoc_lambda.R <pvalues.tsv> <method> <correction> <out.tsv> (got ",
         length(args), " arguments)")
################################################################################
PVALUES    <- args[1]
METHOD     <- args[2]
CORRECTION <- args[3]   # genomic_control | rdadapt_gif | none
OUT_TSV    <- args[4]
################################################################################
if (!CORRECTION %in% c("genomic_control", "rdadapt_gif", "none"))
    stop("assoc_lambda.R: correction must be genomic_control, rdadapt_gif or none, got '",
         CORRECTION, "'")

pv <- fread(PVALUES, colClasses = c(chr = "character"))
traits <- setdiff(colnames(pv), c("SNPID", "chr", "pos"))
if (length(traits) == 0L)
    stop("assoc_lambda.R: ", PVALUES, " has no trait column")

out <- rbindlist(lapply(traits, function(tr) {
    p <- pv[[tr]]
    data.table(method = METHOD, trait = tr, n_tests = sum(!is.na(p)),
               lambda_gc = lambda_gc(p), correction = CORRECTION)
}))
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE)

message("INFO: ", METHOD, " lambda_GC (correction: ", CORRECTION, "): ",
        paste0(out$trait, " = ", signif(out$lambda_gc, 4), collapse = ", "))
