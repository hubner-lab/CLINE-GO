#!/usr/bin/env Rscript
# EMMAX association analysis — single-trait mode
# Runs EMMAX for ONE bioclimatic trait; called once per trait (per-factor caching).
# Uses pre-computed TPED/TFAM and kinship from separate Snakemake rules.

library(dplyr)
library(data.table)
library(magrittr)
library(purrr)
library(stringr)

source("/pipeline/scripts/R/utils/emmax_core.R")

args = commandArgs(trailingOnly=TRUE)
########################
VCF = args[1]
Kbest = args[2] %>% as.numeric
TRAIT_FILE = args[3]  # tsv file with climate/trait values
COVARIATES_FILE = args[4]  # eigenvec file from PCA
TRAIT = args[5]        # Single trait column name (e.g. "bio_1")
INTER_DIR = args[6]
SAMPLES_FILE = args[7]
OUT_FILE = args[8]     # Exact output path for per-trait pvalue TSV
TPED_PREFIX = args[9]  # Pre-computed TPED/TFAM prefix
KINSHIP_FILE = args[10] # Pre-computed .aBN.kinf
SAMPLES_ORDER = if (length(args) >= 11) args[11] else NULL  # Optional: full sample order file
                # (VCF is a coord-valid subset when samples with missing lat/lon were dropped
                # for climate/GEA analysis — enables PCA row subsetting, see load_pca_covariates)
########################

message("INFO: Starting EMMAX analysis (single-trait mode)")
message(paste0("INFO: K = ", Kbest))
message(paste0("INFO: Trait: ", TRAIT))
message(paste0("INFO: TPED prefix: ", TPED_PREFIX))
message(paste0("INFO: Kinship file: ", KINSHIP_FILE))
message(paste0("INFO: Output: ", OUT_FILE))

# Load trait data
trait <- fread(TRAIT_FILE, sep = '\t', header = T)

# Load covariates (PCs from LEA projections or PLINK eigenvec).
# SAMPLES_ORDER + the current VCF's sample list let load_pca_covariates() subset
# PCA rows automatically when the VCF is a coord-valid subset of the full cohort.
covariates <- load_pca_covariates(COVARIATES_FILE, Kbest, SAMPLES_ORDER, read_vcf_samples(VCF))

################################ Functions

FUN_emmax <- function(VCF, trait, covariates, OUT, TPED_PREFIX, KINSHIP_FILE) {

    tfam <- paste0(TPED_PREFIX, '.tfam')

    message("INFO: Read VCF and extract SNPID")
    snpid <- read_vcf_snpids(VCF)

    message("INFO: Read VCF and extract SampleNames")
    SampleName <- read_vcf_samples(VCF)

    # FILES produced
    traitname <- colnames(trait)[1]
    PHEN <- paste0(OUT, 'EMMAX_phenotype_', traitname, '.tsv')
    COVAR <- paste0(OUT, 'EMMAX_covariates_', traitname, '.tsv')
    EMMAXOUT <- paste0(OUT, 'EMMAX_OUT_', traitname)

    # Prepare phenotype file — use FID/IID from TFAM to ensure exact match
    tfam_ids <- fread(tfam, header = FALSE, select = list(character = 1:2),
                      col.names = c("FAMID", "INDID"))
    emmax.phen <- cbind(tfam_ids, trait) %T>%
        write.table(PHEN, sep = '\t', col.names = F, row.names = F, quote = F)

    if (!is.null(covariates)) {
        message('RUN emmax with PCA and kinship corrections')

        # Prepare CV file
        emmax.phen %>%
            dplyr::select(FAMID, INDID) %>%
            dplyr::mutate(smth = 1) %>%  # intercept
            cbind(covariates) %>%
            write.table(COVAR, sep = '\t', col.names = F, row.names = F, quote = F)

        # Run EMMAX
        run_emmax(paste('-v -d 10 -t', TPED_PREFIX,
                        '-p', PHEN,
                        '-k', KINSHIP_FILE,
                        '-c', COVAR,
                        '-o', EMMAXOUT))
    } else {
        message('RUN emmax without PCA correction')
        run_emmax(paste('-v -d 10 -t', TPED_PREFIX,
                        '-p', PHEN,
                        '-k', KINSHIP_FILE,
                        '-o', EMMAXOUT))
    }

    EMMAXOUT.ps <- paste0(EMMAXOUT, '.ps')

    # Load results
    df <- fread(EMMAXOUT.ps) %>%
        setNames(c('SNPID', 'beta', 'SE.beta', traitname)) %>%
        dplyr::mutate(
            SNPID = !!snpid,
            chr = SNPID %>% str_extract("(.*):", group = 1),
            pos = SNPID %>% str_extract(":(.*)", group = 1) %>% as.numeric
        ) %>%
        dplyr::select(SNPID, chr, pos, everything()) %>%
        dplyr::select(-beta, -SE.beta)

    return(df)
}

################################ Main

# Select only the single requested trait
trait <- trait %>% dplyr::select(!!TRAIT)

# Run EMMAX for the single trait
pval_dt <- FUN_emmax(VCF, trait, covariates, INTER_DIR, TPED_PREFIX, KINSHIP_FILE)

message(pval_dt %>% str)

# Save per-trait result to exact output path
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
pval_dt %>% fwrite(OUT_FILE, sep = '\t')

message("INFO: EMMAX analysis complete — wrote ", OUT_FILE)
