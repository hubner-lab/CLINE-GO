#!/usr/bin/env Rscript
# Run per-region GO enrichment analysis using clusterProfiler::enricher()

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(clusterProfiler)
    library(qs)
})

source("/pipeline/scripts/R/lib/gff_parsing.R")
source("/pipeline/scripts/R/lib/enrichment.R")

args = commandArgs(trailingOnly = TRUE)
################
GENES_PER_REGION = args[1]
GFF              = args[2]
GO_FIELD         = args[3]
GFF_FEATURE      = args[4]
TABLES_DIR       = args[5]
INTERMEDIATE_DIR = args[6]
################

message('INFO: Running per-region GO enrichment analysis')
message(paste0('INFO: GO field: ', GO_FIELD))

if (is.null(GO_FIELD) || GO_FIELD == 'NULL' || GO_FIELD == '') {
    message('WARNING: GO_FIELD not specified, skipping enrichment')
    quit(save = "no", status = 0)
}

genes <- fread(GENES_PER_REGION)
if (nrow(genes) == 0) {
    message('WARNING: No genes to process')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Loaded ', nrow(genes), ' region-gene pairs'))

# Extract trait from region_id (format: chr_start-end_trait; trait_from_region_id(), enrichment.R)
genes[, trait := trait_from_region_id(region_id)]

if (anyNA(genes$trait)) {
    message('ERROR: Could not extract a trait from region_id ',
            genes$region_id[is.na(genes$trait)][1], ' (expected <chr>_<start>-<end>_<trait>)')
    quit(save = "no", status = 1)
}

message(paste0('INFO: Traits: ', paste(unique(genes$trait), collapse = ', ')))

# Build term2gene from GFF
bg <- build_term2gene_from_gff(GFF, GFF_FEATURE, GO_FIELD)
if (is.null(bg)) quit(save = "no", status = 0)

term2gene         <- bg$term2gene
term2name         <- bg$term2name
all_genes_with_go <- bg$all_genes_with_go

# Create output directories per trait
for (trait in unique(genes$trait)) {
    dir.create(file.path(TABLES_DIR, trait),       recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(INTERMEDIATE_DIR, trait), recursive = TRUE, showWarnings = FALSE)
}

# Process each region
genes_by_region <- genes[, .(gene_ids = list(unique(gene_id))), by = c("region_id", "trait")]
all_results <- list()
result_counter <- 0L

for (i in seq_len(nrow(genes_by_region))) {
    result <- run_enrichment_for_region(
        region_id       = genes_by_region$region_id[i],
        trait           = genes_by_region$trait[i],
        genes_in_region = genes_by_region$gene_ids[[i]],
        term2gene       = term2gene,
        term2name       = term2name,
        all_genes_with_go = all_genes_with_go
    )

    if (!is.null(result)) {
        result_counter <- result_counter + 1L
        result_dt <- save_enrichment_result(result, TABLES_DIR, INTERMEDIATE_DIR)
        all_results[[result_counter]] <- result_dt
    }
}

# Per-trait summary files
if (length(all_results) > 0) {
    combined <- rbindlist(all_results)
    for (trait in unique(combined$trait)) {
        trait_rows <- combined[trait == trait]
        if (nrow(trait_rows) > 0) {
            summary_file <- file.path(TABLES_DIR, trait,
                                       paste0('Enrichment_', trait, '_summary.tsv'))
            fwrite(trait_rows, summary_file, sep = '\t')
            message(paste0('INFO: Trait ', trait, ': ',
                           uniqueN(trait_rows$region_id), ' regions, ',
                           nrow(trait_rows), ' terms'))
        }
    }
    message(paste0('INFO: Complete: ', result_counter, ' regions with enrichment'))
} else {
    message('WARNING: No regions had significant enrichment')
}

message('INFO: Complete')
