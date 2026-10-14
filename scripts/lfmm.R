#!/usr/bin/env Rscript
# LFMM association analysis — single-trait mode
# Trains LFMM model on LD-pruned data and tests on full dataset for ONE trait.
# Called once per trait (per-factor caching).

library(LEA)
library(dplyr)
library(data.table)
library(WGCNA)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
#################################
LFMM_LD_IMP = args[1]   # LD-pruned imputed LFMM (for training)
LFMM_IMP = args[2]       # Full imputed LFMM (for testing)
CLIMATE = args[3]        # Climate data (scaled)
Kbest = args[4] %>% as.numeric
TRAIT = args[5]          # Single trait column name (e.g. "bio_1")
VCFSNP = args[6]         # SNP positions
OUT_FILE = args[7]       # Exact output path for per-trait pvalue TSV
# Optional arg 8 (preGEA only) — "TRUE"|"FALSE", default TRUE keeps every
# existing GEA call byte-identical. preGEA's K ladder passes FALSE so lambda
# stays informative and the p-value histogram isn't pre-flattened by LEA's
# own genomic-control recalibration (docs/rda_research.md C.0) — with
# genomic.control=TRUE, lambda is forced toward 1 by construction and a
# lambda-vs-K sweep would be flat regardless of whether K is correct.
GENOMIC_CONTROL = if (length(args) >= 8) toupper(args[8]) == "TRUE" else TRUE
#################################

message('INFO: Starting LFMM analysis (single-trait mode)')
message(paste0('INFO: K = ', Kbest))
message(paste0('INFO: Trait: ', TRAIT))
message(paste0('INFO: Output: ', OUT_FILE))

# Load LD-pruned imputed LFMM (for training)
message('INFO: Read LFMM_LD (training data)')
lfmm_ld_imp <- fread(LFMM_LD_IMP, sep = ' ', header = F)
message(lfmm_ld_imp %>% str)

# Load the single requested predictor
message('INFO: Read predictor')
predictors <- fread(CLIMATE, sep = '\t', header = T) %>%
    dplyr::select(!!TRAIT)
message(predictors %>% str)

# Scale predictor
message('INFO: Scale predictor')
predictors[, names(predictors) := lapply(.SD, scale)]
message(predictors %>% str)

# Load full imputed LFMM (for testing)
message('INFO: Read LFMM full (testing data)')
lfmm_imp <- fread(LFMM_IMP, sep = ' ', header = F)

# Load SNP positions
message('INFO: Read VCFSNP')
vcfsnp <- fread(VCFSNP, sep = ' ', header = F) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr),
                  SNPID = paste0(chr, ':', pos)) %>%
    dplyr::select(SNPID, chr, pos)

# Column i of the imputed matrix is row i of the vcfsnp. A count mismatch must stop here:
# data.table recycles a shorter p-value vector into pval_dt when the lengths divide evenly.
if (ncol(lfmm_imp) != nrow(vcfsnp)) {
    stop(sprintf('lfmm.R: %s has %d SNP columns but %s has %d rows -- the SNP labels would not line up',
                 LFMM_IMP, ncol(lfmm_imp), VCFSNP, nrow(vcfsnp)))
}

# Run LFMM2 for the single predictor
message(paste0('INFO: Train LFMM model on LD dataset with ', TRAIT, ' predictor'))
lfmm.model <- lfmm2(lfmm_ld_imp,
                    env = predictors[, ..TRAIT],
                    K = Kbest)
message(lfmm.model %>% str)

message('INFO: Extract p-values from the model (genomic.control=', GENOMIC_CONTROL, ')')
lfmm.res <- lfmm2.test(lfmm.model,
                       input = lfmm_imp,
                       env = predictors[, ..TRAIT],
                       genomic.control = GENOMIC_CONTROL)
message(lfmm.res %>% str)

# Build per-trait result table
pval_dt <- data.table(vcfsnp)
pval_dt[[TRAIT]] <- as.numeric(lfmm.res$pvalues)

# Save per-trait result to exact output path
message('INFO: Save p-values')
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
pval_dt %>%
    fwrite(OUT_FILE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

message(paste0('INFO: LFMM analysis complete — wrote ', OUT_FILE))
