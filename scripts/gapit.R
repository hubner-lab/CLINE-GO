#!/usr/bin/env Rscript
# GAPIT3 association analysis wrapper
# Runs GAPIT models (GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK)
# Uses pre-computed GD/GM, kinship (BN from emmax-kin), and LEA PCA as covariates.
# Outputs standardized pvalue TSVs matching EMMAX/LFMM format.

library(data.table)
library(dplyr)
library(stringr)
library(GAPIT)

args <- commandArgs(trailingOnly = TRUE)
########################
GD_FILE         <- args[1]   # GAPIT numeric GD (Taxa + SNP columns)
GM_FILE         <- args[2]   # GAPIT numeric GM (Name, Chromosome, Position)
TRAIT_FILE      <- args[3]   # Climate site scaled or phenotype file
PCA_FILE        <- args[4]   # LEA PCA projections (space-sep, no header)
KINSHIP_FILE    <- args[5]   # .aBN.kinf from emmax-kin (n x n, no header/rownames)
K_BEST          <- as.integer(args[6])
MODELS          <- strsplit(args[7], ",")[[1]]  # Comma-separated model names
WORKDIR         <- args[8]   # Temp working directory for GAPIT native output
TABLES_DIR      <- args[9]   # Output dir for standardized pvalue TSVs
PREDICTORS      <- strsplit(args[10], ",")[[1]] # Trait column names
SAMPLES_FILE    <- args[11]  # Sample names file (for Taxa matching)
NATIVE_OUTDIR   <- args[12]  # Directory for organized GAPIT native output
SAMPLES_SUBSET  <- if (length(args) >= 13 && args[13] != "NULL") args[13] else NULL
OUTPUT_PREFIX   <- if (length(args) >= 14 && args[14] != "NULL") args[14] else NULL  # For DROP mode: trait name prefix
# Number of PCA covariates to enter as fixed effects (GEA.configs params.n_pcs). This is a
# SEPARATE quantity from K_BEST: K_BEST is sNMF k_best and is only a filename tag here
# (the "_K{k}" in the output name, which the calling Snakemake rule declares as its output
# path), whereas N_PCS controls how many PC columns are used. They were the same value until
# params.n_pcs became honoured, and conflating them made the script write
# "{trait}_pvalues_K0.tsv" while the rule expected "_K5.tsv". Every rule (GEA and GWAS)
# passes it; the K_BEST fallback only serves hand-run calls. A 14-argument rule call
# silently discards the configured n_pcs — tests/python/test_method_param_wiring.py pins it.
N_PCS           <- if (length(args) >= 15 && args[15] != "NULL") as.integer(args[15]) else K_BEST
########################

message("INFO: Starting GAPIT analysis")
message("INFO: Models: ", paste(MODELS, collapse = ", "))
message("INFO: Traits: ", paste(PREDICTORS, collapse = ", "))
message("INFO: K = ", K_BEST, " (filename tag)   n_pcs = ", N_PCS, " (PCA covariates)")

# --- Load GD and GM ---
message("INFO: Loading GD and GM")
myGD <- fread(GD_FILE, header = TRUE)
myGM <- fread(GM_FILE, header = TRUE)
message("INFO: GD: ", nrow(myGD), " samples x ", ncol(myGD) - 1, " SNPs")
message("INFO: GM: ", nrow(myGM), " SNPs")

# --- Load traits ---
message("INFO: Loading trait data")
traits <- fread(TRAIT_FILE, sep = '\t', header = TRUE)

# --- Load PCA covariates ---
message("INFO: Loading PCA projections")
pca_raw <- fread(PCA_FILE, sep = ' ', header = FALSE)
# N_PCS == 0 is a legal setting (the registry declares GAPIT n_pcs with min=0):
# run with kinship only, no PC covariates. It must be special-cased because `1:0`
# is c(1, 0), not an empty range — pca_raw[, 1:0] silently returns one column while
# paste0("PC", 1:0) produces two names, which is where this used to abort.
pca_cols <- if (N_PCS > 0) {
    x <- pca_raw[, 1:N_PCS]
    setnames(x, paste0("PC", 1:N_PCS))
    x
} else {
    message("INFO: n_pcs = 0 — running without PC covariates (kinship only)")
    NULL
}

# --- Load kinship ---
message("INFO: Loading kinship matrix")
kin_raw <- fread(KINSHIP_FILE, header = FALSE)
message("INFO: Kinship dimensions: ", nrow(kin_raw), " x ", ncol(kin_raw))

# --- Load sample names for ordering ---
sample_names <- myGD$Taxa

# --- Subset if DROP mode ---
if (!is.null(SAMPLES_SUBSET)) {
    message("INFO: DROP mode — subsetting to sample list: ", SAMPLES_SUBSET)
    keep_samples <- fread(SAMPLES_SUBSET, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")$IID

    # Subset GD
    keep_idx <- which(sample_names %in% keep_samples)
    myGD <- myGD[keep_idx, ]

    # Subset PCA (NULL when n_pcs = 0 — nothing to subset)
    if (!is.null(pca_cols)) pca_cols <- pca_cols[keep_idx, ]

    # Subset kinship (rows and columns)
    kin_raw <- kin_raw[keep_idx, ..keep_idx]

    sample_names <- myGD$Taxa
    message("INFO: Subsetted to ", length(sample_names), " samples")
}

# --- Build Y (phenotype matrix) ---
# GAPIT expects: Taxa, trait1, trait2, ...
message("INFO: Building phenotype matrix")
myY <- data.table(Taxa = sample_names)

# Match trait data to sample order
# Trait file might have a sample/Taxa column or be in the same order as metadata
if ("sample" %in% colnames(traits)) {
    traits_ordered <- data.table(sample = sample_names) %>%
        merge(traits, by = "sample", sort = FALSE)
    for (pred in PREDICTORS) {
        myY[[pred]] <- traits_ordered[[pred]]
    }
} else {
    # Traits are in same order as samples (from metadata)
    # Need to load samples file to match order
    all_samples <- fread(SAMPLES_FILE, header = TRUE, colClasses = c("site" = "character", "sample" = "character"))
    if ("sample" %in% colnames(all_samples)) {
        traits_with_sample <- cbind(data.table(sample = all_samples$sample), traits)
        traits_ordered <- data.table(sample = sample_names) %>%
            merge(traits_with_sample, by = "sample", sort = FALSE)
        for (pred in PREDICTORS) {
            myY[[pred]] <- traits_ordered[[pred]]
        }
    } else {
        # Assume same order
        for (pred in PREDICTORS) {
            myY[[pred]] <- traits[[pred]]
        }
    }
}

# --- Build CV (covariates) ---
# GAPIT expects: Taxa, PC1, PC2, ...
message("INFO: Building covariates matrix")
myCV <- if (is.null(pca_cols)) NULL else cbind(data.table(Taxa = sample_names), pca_cols)

# --- Build KI (kinship) ---
# GAPIT expects: Taxa column + n x n matrix
message("INFO: Building kinship matrix for GAPIT")
myKI <- data.table(Taxa = sample_names)
myKI <- cbind(myKI, kin_raw)

# --- Create work directory ---
dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
setwd(WORKDIR)

# --- Run GAPIT (one model at a time to avoid GAPIT internal state corruption) ---
message("INFO: Running GAPIT with models: ", paste(MODELS, collapse = ", "))
failed_models <- c()
for (model in MODELS) {
    message("INFO: Running model: ", model)
    tryCatch({
        gapit_result <- GAPIT(
            Y = as.data.frame(myY),
            GD = as.data.frame(myGD),
            GM = as.data.frame(myGM),
            # NULL when n_pcs = 0; with PCA.total = 0 below that leaves GAPIT
            # fitting kinship only, which is the intended "no structure correction" rung.
            CV = if (is.null(myCV)) NULL else as.data.frame(myCV),
            KI = as.data.frame(myKI),
            model = model,
            PCA.total = 0,  # We provide our own PCA via CV
            # GAPIT.Genotype.View builds its LD-decay histogram bins with an unguarded
            # seq(): when n.select > 500 and max(ns) > 1000 it evaluates
            #     seq(max(ns)/4 + 200, max(ns)/3, 200)
            #     seq(max(ns)/3 + 500, max(ns)/2, 500)
            # which need max(ns) >= 2400 and >= 3000 respectively, where ns is the vector
            # of pairwise marker distances inside the LD window. Any dataset whose max
            # pairwise distance falls in 1000 < max(ns) < 3000 therefore hits
            # seq(from, to, by) with to < from and aborts with "wrong sign in 'by'
            # argument", taking the whole model down with it. Measured on the MVP
            # benchmark (50 kb linkage groups): 4 of the first 6 replicates crashed.
            # This view is genotype QC plotting only — marker density, heterozygosity,
            # LD decay. Association results, p-value tables and per-model native output
            # are unaffected by disabling it.
            Geno.View.output = FALSE,
            # GAPIT's default is Random.model = TRUE, which runs GAPIT.RandomModel AFTER
            # the association scan: it refits an lme4 mixed model for every SNP the model
            # called significant, to estimate that SNP's variance component, and writes
            # BLUP/Pred tables and per-model PDFs. Nothing downstream reads any of it —
            # this script's only product is the p-value table assembled below from
            # GAPIT.Association.GWAS_Results.{model}.{trait}.csv, which is written by the
            # scan itself, before RandomModel runs.
            #
            # Its cost scales with the NUMBER OF SIGNIFICANT SNPs, not with the marker
            # count, so it is invisible on a null trait and dominant on a real one.
            # Measured on MVP1232548 cell c3, 1000 individuals x 10 325 SNPs, all eight
            # models in one call: the orthogonal control axis bio_2 (few significant SNPs)
            # completed five models in 20 min, while the confounded axis bio_1 spent
            # 19+ min inside RandomModel for GLM alone and had completed one.
            #
            # Turning it off cannot move a p-value: the scan is already finished when
            # RandomModel is invoked. Verified rather than assumed — the
            # GAPIT.Association.GWAS_Results CSVs for GLM/MLM/CMLM/ECMLM/SUPER on bio_2
            # are byte-identical before and after this change (baseline kept in
            # benchmarks/mvp_eval/probe_p11/randommodel_baseline/).
            Random.model = FALSE,
            file.output = TRUE
        )
    }, error = function(e) {
        message("WARNING: GAPIT model ", model, " failed: ", e$message)
        failed_models <<- c(failed_models, model)
    })
}

setwd(old_wd)

if (length(failed_models) == length(MODELS)) {
    stop("All GAPIT models failed: ", paste(failed_models, collapse = ", "))
}
if (length(failed_models) > 0) {
    message("WARNING: Failed models: ", paste(failed_models, collapse = ", "))
}

# --- Parse GAPIT output and create standardized pvalue tables ---
message("INFO: Parsing GAPIT results")

for (model in MODELS) {
    message("INFO: Processing model: ", model)

    # Collect p-values across traits
    pval_list <- list()

    for (trait in PREDICTORS) {
        # GAPIT output file pattern
        result_file <- file.path(WORKDIR, paste0("GAPIT.Association.GWAS_Results.", model, ".", trait, ".csv"))

        if (!file.exists(result_file)) {
            message("WARNING: Result file not found for ", model, "/", trait, ": ", result_file)
            next
        }

        gwas <- fread(result_file)

        # GAPIT result columns: SNP, Chr, Pos, P.value, MAF, nobs, ...
        pvals <- data.table(
            SNPID = gwas$SNP,
            chr = as.character(gwas$Chr),
            pos = as.integer(gwas$Pos)
        )
        pvals[[trait]] <- gwas$P.value

        pval_list[[trait]] <- pvals
    }

    if (length(pval_list) == 0) {
        message("WARNING: No results for model ", model, ", skipping")
        next
    }

    # Merge all traits into single table
    pval_dt <- pval_list[[1]]
    if (length(pval_list) > 1) {
        for (i in 2:length(pval_list)) {
            pval_dt <- merge(pval_dt, pval_list[[i]], by = c("SNPID", "chr", "pos"), all = TRUE)
        }
    }

    # Sort by chr, pos
    pval_dt <- pval_dt[order(chr, pos)]

    # Write standardized output
    model_dir <- file.path(TABLES_DIR, model)
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

    file_prefix <- if (!is.null(OUTPUT_PREFIX)) OUTPUT_PREFIX else model
    out_file <- file.path(model_dir, paste0(file_prefix, "_pvalues_K", K_BEST, ".tsv"))
    fwrite(pval_dt, out_file, sep = "\t")
    message("INFO: Wrote ", out_file, " (", nrow(pval_dt), " SNPs x ", length(PREDICTORS), " traits)")

    # Copy native GAPIT output files to organized directory
    native_model_dir <- file.path(NATIVE_OUTDIR, model)
    dir.create(native_model_dir, recursive = TRUE, showWarnings = FALSE)

    native_files <- list.files(WORKDIR, pattern = paste0("GAPIT.*", model), full.names = TRUE)
    if (length(native_files) > 0) {
        file.copy(native_files, native_model_dir, overwrite = TRUE)
        message("INFO: Copied ", length(native_files), " native GAPIT files to ", native_model_dir)
    }
}

message("INFO: GAPIT analysis complete")
