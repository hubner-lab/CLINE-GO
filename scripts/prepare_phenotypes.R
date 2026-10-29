#!/usr/bin/env Rscript
# Prepare per-trait phenotype files from metadata for downstream association analysis
# Handles missing values via imputation (MEAN/MEDIAN) or per-trait sample subsetting (DROP)

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
################
METADATA = args[1]       # TSV: site, sample, latitude, longitude, trait1, trait2, ...
MISSING_STRATEGY = args[2]  # MEAN, MEDIAN, or DROP
OUTPUT_DIR = args[3]     # Directory for per-trait outputs
OUTPUT_SUMMARY = args[4] # Path for missing value summary TSV
################

message('INFO: Preparing phenotype files')
message(paste0('INFO: Metadata: ', METADATA))
message(paste0('INFO: Missing strategy: ', MISSING_STRATEGY))
message(paste0('INFO: Output directory: ', OUTPUT_DIR))

# Read metadata
meta <- fread(METADATA, sep = '\t', colClasses = c("site" = "character", "sample" = "character"))
message(paste0('INFO: Loaded ', nrow(meta), ' samples, ', ncol(meta), ' columns'))

# Trait columns: everything after site, sample, latitude, longitude, checked by name (1d4be2;
# common.smk refuses any other header at parse time). [-(1:4)] is empty for a 4-column file,
# where 5:ncol(meta) was c(5, 4) and named longitude a trait.
meta_cols <- c('site', 'sample', 'latitude', 'longitude')
if (!identical(colnames(meta)[1:4], meta_cols)) {
    stop('ERROR: metadata columns 1-4 must be ', paste(meta_cols, collapse = ', '),
         '; got ', paste(colnames(meta)[1:4], collapse = ', '))
}
trait_cols <- colnames(meta)[-(1:4)]
if (length(trait_cols) == 0) {
    stop('ERROR: No phenotype columns found (expected columns 5+ after site, sample, latitude, longitude)')
}
message(paste0('INFO: Traits: ', paste(trait_cols, collapse = ', ')))

# Warn for non-numeric values
for (trait in trait_cols) {
    if (!is.numeric(meta[[trait]])) {
        message(paste0('WARNING: Non-numeric values in trait "', trait, '", coercing to numeric'))
        meta[[trait]] <- suppressWarnings(as.numeric(meta[[trait]]))
    }
}

# Validate strategy
if (!MISSING_STRATEGY %in% c('MEAN', 'MEDIAN', 'DROP')) {
    stop(paste0('ERROR: Unknown MISSING_STRATEGY "', MISSING_STRATEGY, '". Use MEAN, MEDIAN, or DROP'))
}

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Build summary table
summary_dt <- lapply(trait_cols, function(trait) {
    vals <- meta[[trait]]
    n_total <- length(vals)
    n_missing <- sum(is.na(vals))
    n_available <- n_total - n_missing
    missing_samples <- if (n_missing > 0) paste(meta$sample[is.na(vals)], collapse = ',') else ""
    imputed_value <- if (n_missing > 0 && MISSING_STRATEGY != 'DROP') {
        if (MISSING_STRATEGY == 'MEAN') round(mean(vals, na.rm = TRUE), 6)
        else round(median(vals, na.rm = TRUE), 6)
    } else NA_real_

    data.frame(
        trait = trait,
        n_total = n_total,
        n_missing = n_missing,
        n_available = n_available,
        missing_samples = missing_samples,
        strategy = MISSING_STRATEGY,
        imputed_value = imputed_value,
        stringsAsFactors = FALSE
    )
}) %>% do.call(rbind, .)

# Output per strategy
if (MISSING_STRATEGY %in% c('MEAN', 'MEDIAN')) {
    # Impute missing values and output single combined file
    for (trait in trait_cols) {
        n_miss <- sum(is.na(meta[[trait]]))
        if (n_miss > 0) {
            fill_val <- if (MISSING_STRATEGY == 'MEAN') mean(meta[[trait]], na.rm = TRUE)
                        else median(meta[[trait]], na.rm = TRUE)
            meta[[trait]][is.na(meta[[trait]])] <- fill_val
            message(paste0('INFO: Imputed ', n_miss, ' missing values in "', trait, '" with ', MISSING_STRATEGY, ' = ', round(fill_val, 6)))
        }
    }
    out_pheno <- meta %>% dplyr::select(sample, all_of(trait_cols))
    out_path <- file.path(OUTPUT_DIR, 'all_phenotypes.tsv')
    fwrite(out_pheno, out_path, sep = '\t', quote = FALSE)
    message(paste0('INFO: Wrote combined phenotypes (', nrow(out_pheno), ' samples, ', length(trait_cols), ' traits) to ', out_path))

    # Site-level means for piemap (all samples, after imputation)
    for (trait in trait_cols) {
        site_means <- meta[, .(trait_mean = mean(get(trait), na.rm = TRUE)), by = site]
        setnames(site_means, 'trait_mean', trait)
        fwrite(site_means, file.path(OUTPUT_DIR, paste0(trait, '_site_means.tsv')), sep = '\t', quote = FALSE)
        message(paste0('INFO: Wrote site means for "', trait, '" (', nrow(site_means), ' sites)'))
    }

} else {
    # DROP: output per-trait files with only non-missing samples
    for (trait in trait_cols) {
        keep <- !is.na(meta[[trait]])
        meta_trait <- meta[keep, ]
        message(paste0('INFO: Trait "', trait, '": ', sum(keep), '/', nrow(meta), ' samples available'))
        if (sum(keep) < nrow(meta)) {
            message('WARNING: MAF/missingness filtering was computed on the full ',
                    'processing-module cohort BEFORE this per-trait sample drop. Dropping samples for ',
                    'a missing phenotype here can leave those filters inconsistent with the actual GWAS ',
                    'sample set for this trait -- a SNP that passed MAF/missingness genome-wide can ',
                    'become monomorphic or zero-variance once these samples are excluded, and is NOT ',
                    're-filtered. This can produce degenerate model fits (e.g. exact p=0/-log10(p)=Inf) ',
                    'for a small number of SNPs, especially for methods sensitive to zero-variance ',
                    'genotypes (e.g. EMMAX). Check GWAS results for "', trait, '" if they look anomalous ',
                    'for a subset of SNPs.')
        }

        # Phenotype file: sample + trait value
        pheno_dt <- meta_trait %>% dplyr::select(sample, all_of(trait))
        fwrite(pheno_dt, file.path(OUTPUT_DIR, paste0(trait, '_phenotype.tsv')), sep = '\t', quote = FALSE)

        # Sample list: plink --keep format (FID IID), use --const-fid 0 in plink
        samples_list <- data.table(FID = 0, IID = meta_trait$sample)
        fwrite(samples_list, file.path(OUTPUT_DIR, paste0(trait, '_samples.list')), sep = ' ', col.names = FALSE, quote = FALSE)

        # Subsetted metadata: site, sample, latitude, longitude
        meta_sub <- meta_trait %>% dplyr::select(site, sample, latitude, longitude)
        fwrite(meta_sub, file.path(OUTPUT_DIR, paste0(trait, '_metadata.tsv')), sep = '\t', quote = FALSE)

        # Site-level means for piemap (non-missing samples only)
        site_means <- meta_trait[, .(trait_mean = mean(get(trait), na.rm = TRUE)), by = site]
        setnames(site_means, 'trait_mean', trait)
        fwrite(site_means, file.path(OUTPUT_DIR, paste0(trait, '_site_means.tsv')), sep = '\t', quote = FALSE)
        message(paste0('INFO: Wrote site means for "', trait, '" (', nrow(site_means), ' sites)'))
    }
}

# Write summary
fwrite(summary_dt, OUTPUT_SUMMARY, sep = '\t', quote = FALSE)
message(paste0('INFO: Wrote phenotype summary to ', OUTPUT_SUMMARY))
message('INFO: Complete')
