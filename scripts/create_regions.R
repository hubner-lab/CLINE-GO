#!/usr/bin/env Rscript
# Create per-trait and combined regions from significant SNPs using single-linkage clustering
#
# snp_clumping_distance: "auto_per_chromosome" | "auto_genome_wide" | "auto" | integer bp

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
})

source("/pipeline/scripts/R/lib/region_distance.R")
source("/pipeline/scripts/R/lib/regions.R")

args = commandArgs(trailingOnly = TRUE)
options(scipen = 99999)
################
SELECTED_SNPS          = args[1]
CLUMPING_DISTANCE_SPEC = args[2]
OUTPUT_PER_TRAIT       = args[3]
OUTPUT_COMBINED        = args[4]
LD_DECAY_PATH          = if (length(args) >= 5) args[5] else "NULL"
R2_THRESHOLD           = if (length(args) >= 6) as.numeric(args[6]) else 0.2
LD_DECAY_GROUP         = if (length(args) >= 7) args[7] else "All"
TRAIT_PVALUES          = if (length(args) >= 8) args[8] else "NULL"   # selected_snps_per_trait.tsv
################

message('INFO: Creating regions from selected SNPs using single-linkage clustering')

dist_spec <- resolve_clumping_distance(CLUMPING_DISTANCE_SPEC, LD_DECAY_PATH,
                                        R2_THRESHOLD, LD_DECAY_GROUP)

snps <- fread(SELECTED_SNPS, colClasses = c(chr = "character"))

if (nrow(snps) == 0) {
    message('WARNING: No SNPs to process')
    empty_dt <- data.table(
        region_id = character(), trait = character(), chr = character(),
        start = integer(), end = integer(), length = integer(),
        snp_count = integer(), snp_ids = character(),
        methods = character(), min_pvalue = numeric()
    )
    fwrite(empty_dt, OUTPUT_PER_TRAIT, sep = '\t')
    fwrite(empty_dt, OUTPUT_COMBINED,  sep = '\t')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Processing ', nrow(snps), ' SNPs'))

# Per-trait regions. The per-(SNP, trait) p-value table is what lets a trait's region
# report the trait's OWN best p rather than the SNP-level minimum over every trait.
trait_pvalues <- NULL
if (TRAIT_PVALUES != "NULL") {
    if (!file.exists(TRAIT_PVALUES)) stop("per-trait p-value table not found: ", TRAIT_PVALUES)
    trait_pvalues <- fread(TRAIT_PVALUES, colClasses = c(chr = "character"))
}
message('INFO: Creating per-trait regions...')
per_trait <- build_per_trait_regions(snps, dist_spec, trait_pvalues)

if (is.null(per_trait) || nrow(per_trait) == 0) {
    message('WARNING: No per-trait regions created')
    per_trait <- data.table(
        region_id = character(), trait = character(), chr = character(),
        start = integer(), end = integer(), length = integer(),
        snp_count = integer(), snp_ids = character(),
        methods = character(), min_pvalue = numeric(),
        other_traits = character(), other_snp_count = integer()
    )
}

fwrite(per_trait, OUTPUT_PER_TRAIT, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(per_trait), ' per-trait regions to ', OUTPUT_PER_TRAIT))

# Combined regions
message('INFO: Creating combined regions...')
combined <- build_combined_regions(snps, dist_spec)

if (is.null(combined) || nrow(combined) == 0) {
    message('WARNING: No combined regions created')
    combined <- data.table(
        region_id = character(), chr = character(),
        start = integer(), end = integer(), length = integer(),
        snp_count = integer(), snp_ids = character(),
        traits = character(), methods = character(), min_pvalue = numeric()
    )
}

fwrite(combined, OUTPUT_COMBINED, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(combined), ' combined regions to ', OUTPUT_COMBINED))

message('INFO: Complete')
