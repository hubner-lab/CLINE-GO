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
# 5 or 6 arguments, the FIRST being ONE argv slot holding the whole space-separated
# file list. An unquoted list splits into N entries and shifts every later
# positional; with three or more files that lands here, with two it is caught by
# the input-clobber check below. Same failure as combine_pheno_pvalues.R.
if (!(length(args) %in% c(5L, 6L)) || any(is.na(args[1:5])))
    stop("Usage: combine_selected_snps.R '<file1 file2 ...>' <method> <clumping_distance> ",
         "<predictors_csv> <output> [<per_trait_output>]\n  got ", length(args),
         " argument(s): ", paste(args, collapse = " | "),
         "\n  The file list must be ONE quoted argument.")

message('INFO: Combining significant SNPs')
message(paste0('INFO: Strategy: ', METHOD))
message(paste0('INFO: Clumping distance: ', CLUMPING_DISTANCE))

sigSNPs_vec <- str_split(SIGSNPS_FILES, ' ')[[1]]
sigSNPs_vec <- sigSNPs_vec[sigSNPs_vec != ""]

# Second line of defence: never write on top of an input, whatever the argv shape.
clash <- intersect(normalizePath(sigSNPs_vec, mustWork = FALSE),
                   normalizePath(c(OUTPUT, OUTPUT_PER_TRAIT), mustWork = FALSE))
if (length(clash) > 0) {
    stop("Refusing to overwrite an input file: ", paste(clash, collapse = ", "))
}

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
