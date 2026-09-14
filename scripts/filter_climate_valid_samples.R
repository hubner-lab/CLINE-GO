#!/usr/bin/env Rscript
# Further narrow the coord-valid sample set to exclude samples whose raster climate
# extraction returned NA (e.g. a coordinate landing on an ocean/NoData pixel), reported by
# download_climate_present.R. No-op (output == input) when climate_na_excluded.tsv is empty --
# the common case. Samples excluded here are NOT removed from metadata.tsv/metadata_climate.tsv;
# they remain available to GWAS/phenotype/structure/PCA/sNMF and coordinate-only plotting
# (piemaps), which keep using the wider coord-valid set.

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
################
COORD_SAMPLES = args[1]        # coord_valid_samples.list: plink --keep format (FID IID, no header)
METADATA = args[2]             # metadata_climate.tsv: coord-valid samples, VCF order preserved
EXCLUDED = args[3]             # climate_na_excluded.tsv: sample, site, latitude, longitude, reason, distance_km
OUTPUT_SAMPLES_LIST = args[4]  # plink --keep format of climate-valid samples
OUTPUT_METADATA = args[5]      # metadata_climate.tsv further filtered to climate-valid samples
################

message('INFO: Narrowing coord-valid samples to climate-valid samples')

meta <- fread(METADATA, sep = '\t', colClasses = c("site" = "character", "sample" = "character"))
message(paste0('INFO: Loaded ', nrow(meta), ' coord-valid samples'))

excluded <- fread(EXCLUDED, sep = '\t', colClasses = list(character = "sample"))
n_excluded <- nrow(excluded)

if (n_excluded == 0) {
    message('INFO: No climate-NA exclusions -- climate-valid set equals coord-valid set')
    meta_valid <- meta
} else {
    message(paste0('WARNING: Excluding ', n_excluded, ' sample(s) with NA climate values from climate-dependent steps: ',
                   paste(excluded$sample, collapse = ',')))
    meta_valid <- meta[!(meta$sample %in% excluded$sample), ]
}

# Sample list: plink --keep format (FID IID), use --const-fid 0 in plink
samples_list <- data.table(FID = 0, IID = meta_valid$sample)
fwrite(samples_list, OUTPUT_SAMPLES_LIST, sep = ' ', col.names = FALSE, quote = FALSE)

# Filtered metadata (order preserved from input)
fwrite(meta_valid, OUTPUT_METADATA, sep = '\t', quote = FALSE)

message(paste0('INFO: Wrote ', nrow(meta_valid), ' climate-valid samples to ', OUTPUT_SAMPLES_LIST))
message(paste0('INFO: Wrote filtered metadata to ', OUTPUT_METADATA))
message('INFO: Complete')
