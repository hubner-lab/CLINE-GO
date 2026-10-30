#!/usr/bin/env Rscript
# sig_threshold_table.R — the per-trait significance threshold behind one GEA/GWAS method's
# sig-SNP table, written beside its p-value table as {method}_thresholds_K{k}.tsv (bb27e1).
#
# One row per trait column of the p-value table:
#   rule, value     the method's configs entry (adjust and threshold)
#   status          compute_pval_threshold()'s verdict; not "ok" means no threshold and no call
#   threshold       the p-value cutoff the rule gave (NA unless status is ok)
#   n_tested        the non-NA p-values it was derived from: Bonferroni's family size, which is
#                   counted per trait column
#   n_na_dropped    the NA p-values left out
#   n_significant   p <= threshold: the trait's rows in the method's sig-SNP table
# Every column comes from find_significant_snps_per_trait(), the function find_sig_snps.R
# writes the sig-SNP table with, so the two tables cannot disagree. The file name carries no
# rule: a changed adjust rewrites this file instead of leaving a stale variant beside it.
#
#   Rscript sig_threshold_table.R <pvalues.tsv> <adjust> <method> <out.tsv>
#     adjust: <rule>_<value>, as in the configs entry (e.g. bonf_0.05)

suppressPackageStartupMessages({
    library(data.table)
    library(qvalue)
    library(parallel)
})

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/lib/sig_snps.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L)
    stop("usage: sig_threshold_table.R <pvalues.tsv> <adjust> <method> <out.tsv> (got ",
         length(args), " arguments)")
################################################################################
PVALUES <- args[1]
ADJUST  <- args[2]   # e.g. bonf_0.05
METHOD  <- args[3]
OUT_TSV <- args[4]
################################################################################
# Split at the first underscore, as find_sig_snps.R does: a rule name never holds one.
RULE  <- sub("_.*$", "", ADJUST)
VALUE <- as.numeric(sub("^[^_]*_", "", ADJUST))
if (!grepl("_", ADJUST, fixed = TRUE) || is.na(VALUE))
    stop("sig_threshold_table.R: adjust must be <rule>_<value>, got '", ADJUST, "'")

pv  <- fread(PVALUES, colClasses = c(chr = "character"))
res <- find_significant_snps_per_trait(pv, RULE, VALUE, 1)

n_sig <- res$sig_snps[, .N, by = trait]
out <- res$diagnostics[, .(method = METHOD, trait, rule = RULE, value = VALUE, status,
                           threshold, n_tested, n_na_dropped)]
out[, n_significant := n_sig$N[match(trait, n_sig$trait)]]
out[is.na(n_significant), n_significant := 0L]
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE)

message("INFO: ", METHOD, " ", ADJUST, " thresholds written for ", nrow(out), " trait(s): ",
        paste0(out$trait, " (", out$status, ", ", out$n_significant, " significant)",
               collapse = ", "))
