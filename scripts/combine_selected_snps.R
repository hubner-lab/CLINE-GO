#!/usr/bin/env Rscript
# Combine significant SNPs from different methods into Selected_SNPs table

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
    library(magrittr)
})

source("/pipeline/scripts/R/lib/regions.R")          # distance_for_chr()
source("/pipeline/scripts/R/lib/region_distance.R")  # resolve_clumping_distance_strict()
source("/pipeline/scripts/R/lib/combine_sigsnps.R")

args = commandArgs(trailingOnly = TRUE)
options(scipen = 99999)
################
SIGSNPS_FILES       = args[1]  # space-separated list of files
METHOD              = args[2]  # Union | Cross-method | Cross-method per-trait (or legacy aliases)
CLUMPING_SPEC       = args[3]  # agreement distance: bp, or auto_genome_wide | auto_per_chromosome
PREDICTORS_SELECTED = str_split(args[4], ',')[[1]]
OUTPUT              = args[5]
OUTPUT_PER_TRAIT    = if (length(args) >= 6) args[6] else "NULL"  # selected_snps_per_trait.tsv
LD_DECAY_PATH       = if (length(args) >= 7) args[7] else "NULL"  # resolves an auto_* spec
LD_RULE             = if (length(args) >= 8) args[8] else "critical"  # clumping_ld_rule
LD_DECAY_GROUP      = if (length(args) >= 9) args[9] else "All"
################
# 5, 6 or 9 arguments, the FIRST being ONE argv slot holding the whole space-separated
# file list. An unquoted list splits into N entries and shifts every later
# positional; with three or more files that lands here, with two it is caught by
# the input-clobber check below. Same failure as combine_pheno_pvalues.R.
if (!(length(args) %in% c(5L, 6L, 9L)) || any(is.na(args[1:5])))
    stop("Usage: combine_selected_snps.R '<file1 file2 ...>' <method> <clumping_distance> ",
         "<predictors_csv> <output> [<per_trait_output> [<ld_decay_table> <ld_rule> ",
         "<ld_decay_group>]]\n  got ", length(args),
         " argument(s): ", paste(args, collapse = " | "),
         "\n  The file list must be ONE quoted argument.")

# The same distance as the region clumping: an auto_* spec is resolved from the
# LD-decay table exactly as create_regions.R resolves it.
CLUMPING_DISTANCE <- resolve_clumping_distance_strict(CLUMPING_SPEC, LD_DECAY_PATH,
                                                      LD_RULE, LD_DECAY_GROUP)

message('INFO: Combining significant SNPs')
message(paste0('INFO: Strategy: ', METHOD))
message(paste0('INFO: Clumping distance: ', format_clumping_distance(CLUMPING_DISTANCE)))

sigSNPs_vec <- str_split(SIGSNPS_FILES, ' ')[[1]]
sigSNPs_vec <- sigSNPs_vec[sigSNPs_vec != ""]

# Second line of defence: never write on top of an input, whatever the argv shape.
clash <- intersect(normalizePath(sigSNPs_vec, mustWork = FALSE),
                   normalizePath(c(OUTPUT, OUTPUT_PER_TRAIT), mustWork = FALSE))
if (length(clash) > 0) {
    stop("Refusing to overwrite an input file: ", paste(clash, collapse = ", "))
}

message(paste0('INFO: Input files: ', paste(sigSNPs_vec, collapse = ', ')))

# Load tables
sigSNPs_lst <- lapply(sigSNPs_vec, function(x) {
    dt <- fread(x, colClasses = c(chr = "character"))
    if (nrow(dt) == 0) {
        return(data.table(SNPID = character(), chr = character(),
                          pos = integer(), trait = character(),
                          method = character(), pvalue = numeric()))
    }
    if (!'method' %in% names(dt)) stop("ERROR: ", x, " has no `method` column")
    dt
})

# Each table's method is the one its own `method` column names (find_sig_snps.R writes it),
# not the folder it sits in (caf4df). A table without rows names none: only then is the
# folder used, tables/methods/<method>/ as find_sig_snps writes it. Argv order = column order.
methods_vec <- vapply(seq_along(sigSNPs_lst), function(i) {
    m <- unique(sigSNPs_lst[[i]]$method)
    if (length(m) > 1) {
        stop("ERROR: ", sigSNPs_vec[i], " holds more than one method: ", paste(m, collapse = ", "))
    }
    if (length(m) == 1) m else basename(dirname(sigSNPs_vec[i]))
}, character(1))
if (anyDuplicated(methods_vec)) {
    dup <- methods_vec[duplicated(methods_vec)]
    stop("ERROR: two input files are the same method (", paste(unique(dup), collapse = ", "), "): ",
         paste(sigSNPs_vec[methods_vec %in% dup], collapse = ", "))
}
names(sigSNPs_lst) <- methods_vec

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
