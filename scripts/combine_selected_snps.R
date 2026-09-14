#!/usr/bin/env Rscript
# Combine significant SNPs from different methods into Selected_SNPs table

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
    library(magrittr)
})

source("/pipeline/scripts/R/lib/combine_sigsnps.R")

args = commandArgs(trailingOnly = TRUE)
options(scipen = 99999)
################
SIGSNPS_FILES       = args[1]  # space-separated list of files
METHOD              = args[2]  # Union | Cross-method | Cross-method per-trait (or legacy aliases)
CLUMPING_DISTANCE   = as.numeric(args[3])
PREDICTORS_SELECTED = str_split(args[4], ',')[[1]]
OUTPUT              = args[5]
OUTPUT_PER_TRAIT    = if (length(args) >= 6) args[6] else "NULL"  # selected_snps_per_trait.tsv
################
if (length(args) < 5 || any(is.na(args[1:5])))
    stop("combine_selected_snps.R needs 5 positional args (+ optional per-trait output); got ", length(args))

message('INFO: Combining significant SNPs')
message(paste0('INFO: Strategy: ', METHOD))
message(paste0('INFO: Clumping distance: ', CLUMPING_DISTANCE))

sigSNPs_vec <- str_split(SIGSNPS_FILES, ' ')[[1]]
sigSNPs_vec <- sigSNPs_vec[sigSNPs_vec != ""]

# Extract method names from file paths
methods_vec <- sapply(str_split(sigSNPs_vec, '/'), function(x) x[length(x) - 1])
sigSNPs_vec <- setNames(sigSNPs_vec, methods_vec)

message(paste0('INFO: Input files: ', paste(sigSNPs_vec, collapse = ', ')))

# Load tables
sigSNPs_lst <- lapply(sigSNPs_vec, function(x) {
    dt <- fread(x, colClasses = c(chr = "character"))
    if (nrow(dt) == 0) {
        return(data.table(SNPID = character(), chr = character(),
                          pos = integer(), trait = character(),
                          method = character(), pvalue = numeric()))
    }
    dt
}) %>% setNames(methods_vec)

combined <- combine_sigsnps_with_traits(sigSNPs_lst, METHOD, CLUMPING_DISTANCE, PREDICTORS_SELECTED)

combined$snps %>% fwrite(OUTPUT, sep = '\t')
message(paste0('INFO: Saved to ', OUTPUT))

# Per-(SNP, trait) minimum p — the trait-scoped companion of selected_snps.tsv's
# trait-agnostic min_pvalue, read by create_regions.R for the per-trait region table
# (audit 2026-09-13 SC3).
if (OUTPUT_PER_TRAIT != "NULL") {
    combined$trait_pvalues %>% fwrite(OUTPUT_PER_TRAIT, sep = '\t')
    message(paste0('INFO: Saved per-trait p-values to ', OUTPUT_PER_TRAIT))
}
message('INFO: Complete')
