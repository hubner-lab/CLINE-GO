# emmax_core.R — shared helpers for EMMAX association scripts
# Usage: source("/pipeline/scripts/R/utils/emmax_core.R")
# Requires: data.table, dplyr loaded by the sourcing script

EMMAX_BIN <- "/pipeline/scripts/emmax-intel64"
EMMAX_RUN <- "/pipeline/scripts/emmax_run.sh"

# Run EMMAX through the preflight wrapper (scripts/emmax_run.sh) and STOP on a
# non-zero exit.
#
# Two things this fixes. (1) The wrapper refuses to launch the x86-64 EMMAX
# binary on a host that cannot run it natively — under x86-on-arm64 emulation
# EMMAX's Intel OpenMP/MKL runtime hangs forever instead of failing, which
# stalls Snakemake and orphans the container. (2) Every call site used bare
# system() and ignored its status, so an EMMAX failure surfaced later as a
# missing .ps file: emmax.R then errored inside fread() with an unrelated
# message, and emmax_phenotypes.R downgraded it to a warning and returned NULL,
# silently dropping that trait from the association results.
run_emmax <- function(args) {
    cmd <- paste(shQuote(EMMAX_RUN), shQuote(EMMAX_BIN), args)
    message(paste0("INFO: Running: ", cmd))
    status <- system(cmd)
    if (status != 0) {
        stop(paste0("ERROR: EMMAX exited with status ", status,
                    " — see the messages above for the cause. Command: ", cmd))
    }
    invisible(status)
}

# Extract CHROM:POS SNP IDs from a VCF (reads only non-header lines, first 2 columns).
read_vcf_snpids <- function(vcf_path) {
    data.table::fread(cmd = paste('grep -v "##"', vcf_path, '| cut -f1-2')) %>%
        setNames(c('CHROM', 'POS')) %>%
        dplyr::mutate(SNPID = paste0(CHROM, ':', POS)) %>%
        .$SNPID
}

# Extract sample names from the VCF #CHROM header line (columns 10+).
read_vcf_samples <- function(vcf_path) {
    data.table::fread(cmd = paste("head -1000", vcf_path, "| grep '#CHROM' | cut -f10-"),
                      header = FALSE) %>%
        as.character()
}

# Load PCA projections from a LEA or PLINK eigenvec file.
# k:                  number of PCs to use as covariates
# samples_order_path: path to a single-column file listing full sample order (DROP mode)
# vcf_samples:        character vector of samples currently in the VCF (DROP mode)
# When samples_order_path + vcf_samples are provided and the PCA has more rows than the VCF,
# rows are subset to the VCF samples (DROP mode for phenotype association).
load_pca_covariates <- function(pca_path, k, samples_order_path = NULL, vcf_samples = NULL) {
    # k <= 0 means "no PCA covariates" (e.g. EMMAX kinship-only reference panel in
    # preGEA's #PC ladder). Must be checked BEFORE format autodetection below:
    # the LEA branch (`1:k` -> `1:0` -> c(1,0)) errors on index 0, but the PLINK
    # branch (`3:(2+k)` -> `3:2` -> c(3,2)) does NOT error — it silently selects
    # the wrong two columns in reversed order, which is a far worse failure mode
    # (plausible-looking wrong output instead of a crash).
    if (k <= 0) {
        message("INFO: k=0 — no PCA covariates (kinship-only model)")
        return(NULL)
    }

    cov_raw <- data.table::fread(pca_path, sep = ' ', header = FALSE)

    # Autodetect PLINK eigenvec format (FID + IID prefix columns)
    if (ncol(cov_raw) > 2 &&
        length(unique(cov_raw[[1]])) == 1 &&
        length(unique(cov_raw[[2]])) == nrow(cov_raw)) {
        message("INFO: Detected PLINK eigenvec format")
        pca_cols <- 3:(2 + k)
    } else {
        message("INFO: Detected LEA projections format")
        pca_cols <- 1:k
    }

    # DROP mode: subset rows when PCA covers the full cohort but VCF is per-trait subset
    if (!is.null(vcf_samples) && !is.null(samples_order_path) &&
        nrow(cov_raw) > length(vcf_samples)) {
        message("INFO: PCA has more rows than VCF samples — subsetting (DROP mode)")
        all_samples <- data.table::fread(samples_order_path, header = FALSE,
                                         colClasses = "character")$V1
        message(paste0("INFO: Full sample order: ", length(all_samples), " samples"))
        keep_idx <- which(all_samples %in% vcf_samples)
        cov_raw  <- cov_raw[keep_idx, ]
        message(paste0("INFO: Subsetted PCA to ", nrow(cov_raw), " rows"))
    }

    if (!is.null(vcf_samples) && nrow(cov_raw) != length(vcf_samples)) {
        stop(paste0("ERROR: PCA row count (", nrow(cov_raw),
                    ") does not match VCF sample count (", length(vcf_samples),
                    "). Provide samples_order_path for DROP mode."))
    }

    cov_raw[, pca_cols, with = FALSE] %>% setNames(paste0('PC', 1:k))
}
