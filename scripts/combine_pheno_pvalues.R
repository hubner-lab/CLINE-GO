#!/usr/bin/env Rscript
# Combine per-trait p-value files into a single table
# Used in DROP mode where each trait's EMMAX was run on different sample sets

library(data.table)
library(dplyr)
library(stringr)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
# EXACTLY three arguments, the FIRST being ONE argv slot holding the whole
# space-separated file list. Without this an unquoted "{params.files_str}" in a
# rule splits into N entries, every positional shifts left, and the script writes
# its combined table OVER one of its own inputs at exit 0 — invisible to Snakemake,
# since the damaged file is an input, not an output. Reproduced before the guard:
# "Saved combined p-values to .../p_flowering_time.tsv", which was an input.
if (length(args) != 3L) {
    stop("Usage: combine_pheno_pvalues.R '<file1 file2 ...>' <pvalue_out> <qvalue_out>\n",
         "  got ", length(args), " argument(s): ", paste(args, collapse = " | "),
         "\n  The file list must be ONE quoted argument.")
}
################
FILES_STR = args[1]   # space-separated list of per-trait p-value TSV paths
OUTPUT = args[2]       # output file path for combined p-values
QVAL_OUTPUT = args[3]  # output file path for combined q-values
################

# Parse input files
files <- FILES_STR %>% str_split(' ') %>% unlist
files <- files[files != ""]

# Second line of defence: never write on top of an input, whatever the argv shape.
clash <- intersect(normalizePath(files, mustWork = FALSE),
                   normalizePath(c(OUTPUT, QVAL_OUTPUT), mustWork = FALSE))
if (length(clash) > 0) {
    stop("Refusing to overwrite an input file: ", paste(clash, collapse = ", "))
}

message(paste0('INFO: Combining ', length(files), ' per-trait p-value files'))
message(paste0('INFO: Files: ', paste(files, collapse = ', ')))

# Read and merge all files by full outer join on SNPID, chr, pos
pval_dt <- lapply(files, fread) %>%
    Reduce(function(x, y) merge(x, y, by = c("SNPID", "chr", "pos"), all = TRUE), .)

traits <- pval_dt %>% dplyr::select(-SNPID, -chr, -pos) %>% colnames
message(paste0('INFO: Traits: ', paste(traits, collapse = ', ')))
message(paste0('INFO: Combined table: ', nrow(pval_dt), ' SNPs x ', length(traits), ' traits'))

# Compute q-values per trait (fallback to BH if qvalue fails)
message("INFO: Computing q-values")
qval_dt <- lapply(pval_dt %>% dplyr::select(-SNPID, -chr, -pos), function(pvals) {
    tryCatch(
        qvalue(pvals)$qvalues,
        error = function(e) {
            message(paste0("WARNING: qvalue failed (", e$message, "), using BH adjustment"))
            p.adjust(pvals, method = 'BH')
        }
    )
}) %>%
    do.call(cbind, .) %>%
    as.data.table() %>%
    cbind(pval_dt %>% dplyr::select(SNPID, chr, pos), .)

# Write output
fwrite(pval_dt, OUTPUT, sep = '\t')
message(paste0('INFO: Saved combined p-values to ', OUTPUT))

fwrite(qval_dt, QVAL_OUTPUT, sep = '\t')
message(paste0('INFO: Saved combined q-values to ', QVAL_OUTPUT))
