#!/usr/bin/env Rscript
# EMMAX association analysis for phenotype traits
# Receives pre-computed TPED, kinship, and PCA projections from separate rules.
# Formats inputs for EMMAX, runs per-trait, parses output, computes q-values.

library(data.table)
library(dplyr)
library(magrittr)
library(purrr)
library(stringr)
library(qvalue)

source("/pipeline/scripts/R/utils/emmax_core.R")
source("/pipeline/scripts/R/utils/pval_threshold.R")

args = commandArgs(trailingOnly=TRUE)
########################
VCF              = args[1]                # Filtered VCF (possibly per-trait subset for DROP mode)
TPED_PREFIX      = args[2]                # Path prefix for pre-computed .tped/.tfam
KINSHIP_FILE     = args[3]                # Pre-computed .aIBS.kinf
PCA_PROJECTIONS  = args[4]                # LEA PCA projections (space-sep, no header)
K_BEST           = args[5] %>% as.numeric # Number of PCs as covariates
PHENOTYPE_FILE   = args[6]                # TSV: sample, trait1, [trait2, ...]
TABLES_DIR       = args[7]                # Output directory
SAMPLES_ORDER    = if (length(args) >= 8) args[8] else NULL  # Optional: full sample order file
OUTPUT_PVALUES   = if (length(args) >= 9) args[9] else NULL  # Optional: explicit output p-values path
TRAIT_SUBSET     = if (length(args) >= 10) args[10] else NULL # Optional: run only this trait (per-factor caching)
########################

message("INFO: Starting EMMAX phenotype analysis")
message(paste0("INFO: VCF: ", VCF))
message(paste0("INFO: TPED prefix: ", TPED_PREFIX))
message(paste0("INFO: Kinship: ", KINSHIP_FILE))
message(paste0("INFO: PCA projections: ", PCA_PROJECTIONS))
message(paste0("INFO: K = ", K_BEST))
message(paste0("INFO: Phenotype file: ", PHENOTYPE_FILE))
message(paste0("INFO: Samples order: ", ifelse(is.null(SAMPLES_ORDER), "not provided", SAMPLES_ORDER)))

# Create output directory
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Extract SNP IDs and sample names from VCF ---
message("INFO: Extracting SNP IDs from VCF")
snpid <- read_vcf_snpids(VCF)

message("INFO: Extracting sample names from VCF")
SampleName <- read_vcf_samples(VCF)
message(paste0("INFO: ", length(SampleName), " samples in VCF"))

# --- Load PCA projections (autodetect format, optional DROP-mode subsetting) ---
message("INFO: Loading PCA projections")
covariates <- load_pca_covariates(PCA_PROJECTIONS, K_BEST, SAMPLES_ORDER, SampleName)
message(paste0("INFO: Covariates matrix: ", nrow(covariates), " rows x ", ncol(covariates), " cols"))

# --- Read phenotype file ---
message("INFO: Reading phenotype file")
pheno <- fread(PHENOTYPE_FILE, sep = '\t', header = TRUE, colClasses = c("sample" = "character"))
trait_cols <- colnames(pheno)[colnames(pheno) != 'sample']
# Per-factor caching: subset to a single trait when TRAIT_SUBSET is provided
if (!is.null(TRAIT_SUBSET)) {
    if (!(TRAIT_SUBSET %in% trait_cols))
        stop(paste0("TRAIT_SUBSET '", TRAIT_SUBSET, "' not found in phenotype file. Available: ",
                    paste(trait_cols, collapse = ", ")))
    trait_cols <- TRAIT_SUBSET
    message(paste0("INFO: Single-trait mode — processing only: ", TRAIT_SUBSET))
}
message(paste0("INFO: Traits: ", paste(trait_cols, collapse = ', ')))
message(paste0("INFO: Phenotype samples: ", nrow(pheno)))

# Match phenotype sample order to VCF sample order
pheno_ordered <- data.table(sample = SampleName) %>%
    left_join(pheno, by = 'sample')

if (any(is.na(pheno_ordered$sample))) {
    stop("ERROR: Some VCF samples not found in phenotype file")
}

# --- Run EMMAX per trait ---
work_dir <- dirname(TPED_PREFIX)
results <- lapply(trait_cols, function(trait) {
    message(paste0("INFO: Processing trait: ", trait))

    trait_vals <- pheno_ordered[[trait]]
    if (all(is.na(trait_vals))) {
        message(paste0("WARNING: All values NA for trait ", trait, ", skipping"))
        return(NULL)
    }

    # EMMAX phenotype file: FID IID value (no header)
    phen_file <- file.path(work_dir, paste0('EMMAX_phenotype_', trait, '.tsv'))
    data.table(FID = SampleName, IID = SampleName, value = trait_vals) %>%
        fwrite(phen_file, sep = '\t', col.names = FALSE)

    # EMMAX covariate file: FID IID 1 PC1 PC2 ... (no header)
    covar_file <- file.path(work_dir, paste0('EMMAX_covariates_', trait, '.tsv'))
    data.table(FID = SampleName, IID = SampleName, intercept = 1) %>%
        cbind(covariates) %>%
        fwrite(covar_file, sep = '\t', col.names = FALSE)

    # Run EMMAX
    out_prefix <- file.path(work_dir, paste0('EMMAX_OUT_', trait))
    run_emmax(paste('-v -d 10',
                    '-t', TPED_PREFIX,
                    '-p', phen_file,
                    '-k', KINSHIP_FILE,
                    '-c', covar_file,
                    '-o', out_prefix))

    # Parse .ps output
    ps_file <- paste0(out_prefix, '.ps')
    if (!file.exists(ps_file)) {
        message(paste0("WARNING: EMMAX output not found for trait ", trait))
        return(NULL)
    }

    df <- fread(ps_file) %>%
        setNames(c('SNPID', 'beta', 'SE.beta', trait)) %>%
        dplyr::mutate(
            SNPID = !!snpid,
            chr = SNPID %>% str_extract("(.*):", group = 1),
            pos = SNPID %>% str_extract(":(.*)", group = 1) %>% as.numeric
        ) %>%
        dplyr::select(SNPID, chr, pos, all_of(trait))

    message(paste0("INFO: ", trait, " — ", nrow(df), " SNPs"))
    return(df)
})

# Remove NULLs and merge
results <- results[!sapply(results, is.null)]

if (length(results) == 0) {
    stop("ERROR: No traits produced results")
}

pval_dt <- results %>%
    reduce(function(x, y) left_join(x, y, by = c("SNPID", "chr", "pos")))

message(paste0("INFO: Combined p-value table: ", nrow(pval_dt), " SNPs x ", length(trait_cols), " traits"))

# --- Compute q-values (with safe fallback to BH if qvalue fails) ---
message("INFO: Computing q-values")
qval_dt <- lapply(pval_dt %>% dplyr::select(-SNPID, -chr, -pos), function(pvals) {
    compute_qvalues_safe(pvals)
}) %>%
    do.call(cbind, .) %>%
    as.data.table() %>%
    cbind(pval_dt %>% dplyr::select(SNPID, chr, pos), .)

# --- Write output ---
if (!is.null(OUTPUT_PVALUES)) {
    pval_path <- OUTPUT_PVALUES
    qval_path <- sub('pvalues', 'qvalues', OUTPUT_PVALUES)
} else {
    emmax_dir <- file.path(TABLES_DIR, "EMMAX")
    dir.create(emmax_dir, recursive = TRUE, showWarnings = FALSE)
    pval_path <- file.path(emmax_dir, paste0("EMMAX_pvalues_K", K_BEST, ".tsv"))
    qval_path <- file.path(emmax_dir, paste0("EMMAX_qvalues_K", K_BEST, ".tsv"))
}

fwrite(pval_dt, pval_path, sep = '\t')
fwrite(qval_dt, qval_path, sep = '\t')

message(paste0("INFO: Wrote p-values to ", pval_path))
message(paste0("INFO: Wrote q-values to ", qval_path))
message("INFO: EMMAX phenotype analysis complete")
