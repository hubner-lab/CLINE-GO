#!/usr/bin/env Rscript
# Filter samples to those with valid (non-NA) coordinates for climate/GEA analysis.
# Unlike phenotypes, climate values cannot be meaningfully imputed, so missing
# coordinates are always dropped -- no MEAN/MEDIAN option, no config knob.
# Samples with missing coordinates are NOT removed from the main metadata.tsv;
# they remain available to GWAS/phenotype and other non-spatial analyses.

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
################
METADATA = args[1]             # TSV: site, sample, latitude, longitude, [trait1, trait2, ...]
OUTPUT_SAMPLES_LIST = args[2]  # plink --keep format (FID IID, no header) of coord-valid samples
OUTPUT_METADATA = args[3]      # metadata.tsv filtered to coord-valid samples, VCF order preserved
OUTPUT_SUMMARY = args[4]       # coord_missing_summary.tsv
################

message('INFO: Filtering samples with valid coordinates')
message(paste0('INFO: Metadata: ', METADATA))

meta <- fread(METADATA, sep = '\t', colClasses = c("site" = "character", "sample" = "character"))
message(paste0('INFO: Loaded ', nrow(meta), ' samples'))

if (!all(c('latitude', 'longitude') %in% colnames(meta))) {
    stop('ERROR: Metadata missing required "latitude"/"longitude" columns')
}

keep <- !is.na(meta$latitude) & !is.na(meta$longitude)
n_total <- nrow(meta)
n_available <- sum(keep)
n_missing <- n_total - n_available
missing_samples <- if (n_missing > 0) paste(meta$sample[!keep], collapse = ',') else ""

message(paste0('INFO: Climate coordinates: ', n_available, '/', n_total, ' samples have valid coordinates'))
if (n_missing > 0) {
    message(paste0('WARNING: Dropping ', n_missing, ' sample(s) with missing coordinates from climate/GEA analysis: ', missing_samples))
    message('WARNING: MAF/missingness filtering was computed on the full processing-module ',
            'cohort BEFORE this drop. Removing samples here can leave those filters inconsistent with ',
            'the actual GEA sample set -- a SNP that passed MAF/missingness genome-wide can become ',
            'monomorphic or zero-variance once these samples are excluded, and is NOT re-filtered. This ',
            'can produce degenerate model fits (e.g. exact p=0/-log10(p)=Inf) for a small number of SNPs, ',
            'especially for methods sensitive to zero-variance genotypes (e.g. EMMAX). Check GEA results ',
            'if they look anomalous for a subset of SNPs -- consider filtering suspect SNPs upstream on ',
            'the actual analysis cohort if this occurs.')
}

meta_valid <- meta[keep, ]

# Sample list: plink --keep format (FID IID), use --const-fid 0 in plink
samples_list <- data.table(FID = 0, IID = meta_valid$sample)
fwrite(samples_list, OUTPUT_SAMPLES_LIST, sep = ' ', col.names = FALSE, quote = FALSE)

# Filtered metadata (VCF order preserved -- input meta was already ordered by align_metadata)
fwrite(meta_valid, OUTPUT_METADATA, sep = '\t', quote = FALSE)

# Summary (single-row table, mirrors prepare_phenotypes.R's per-trait summary)
summary_dt <- data.table(
    n_total = n_total,
    n_available = n_available,
    n_missing = n_missing,
    missing_samples = missing_samples
)
fwrite(summary_dt, OUTPUT_SUMMARY, sep = '\t', quote = FALSE)

message(paste0('INFO: Wrote ', nrow(meta_valid), ' coord-valid samples to ', OUTPUT_SAMPLES_LIST))
message(paste0('INFO: Wrote filtered metadata to ', OUTPUT_METADATA))
message(paste0('INFO: Wrote coordinate summary to ', OUTPUT_SUMMARY))
message('INFO: Complete')
