#!/usr/bin/env Rscript
# rda_param_sweep.R -- isolate the three parameters that separate CLINE-GO's RDA GEA scan
# from the one Lotterhos 2023 ran on the same seeds:
#
#   1. scale=        vegan::rda response scaling. Ours TRUE (correlation-based); BOTH
#                    Lotterhos and the Capblancq & Forester 2021 tutorial leave it at FALSE.
#   2. condition_pcs Ours resolves through the @k_best sentinel (5 here); Lotterhos runs two
#                    arms, 0 PCs and Condition(PC1+PC2) = 2 PCs.
#   3. pmax          Ours combines a partial and an unconstrained fit via
#                    pmax(p_partial, p_unconstrained) -- the continuous-p analogue of the
#                    tutorial's Reduce(intersect, ...). Lotterhos reports each arm separately.
#
# This is a BENCHMARK PROBE, not pipeline code. It replicates scripts/rda.R's data loading
# and re-uses the pipeline's own rdadapt() and load_pca_covariates() verbatim, so the only
# things that vary across cells are the three parameters above. It writes nothing into any
# {PROJECT}_results/ tree.
#
# Usage:
#   Rscript benchmarks/rda_param_sweep.R --seed=1232548 [--project=MVP1232548] [--k=2]

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(dplyr)
  library(robust)
  library(qvalue)
})

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || identical(x, "")) y else x
A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
  a <- sub("^--", "", a)
  A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
SEED_ID <- A$seed
PROJECT <- A$project %||% paste0("MVP", SEED_ID)
K_FIXED <- as.integer(A$k %||% "2")   # Lotterhos and the tutorial both hard-code rdadapt(fit, 2)
if (is.null(SEED_ID)) stop("Usage: rda_param_sweep.R --seed=N")

ROOT <- "/pipeline"
source(file.path(ROOT, "scripts/R/lib/rdadapt.R"))
source(file.path(ROOT, "scripts/R/utils/emmax_core.R"))

RES  <- file.path(ROOT, paste0(PROJECT, "_results"))
INT  <- file.path(RES, "_intermediate")
WORK <- Sys.glob(file.path(RES, "_work", "maf*", "ld*"))[1]
WORKP<- dirname(WORK)

LFMM_FULL     <- file.path(INT, "climate_subset/lfmm_imp_full_climate.lfmm")
CLIMATE       <- file.path(RES, "climate/tables/present/climate_present_site_scaled.tsv")
VCFSNP_FULL   <- Sys.glob(file.path(WORKP, "*.vcfsnp"))[1]
PCA           <- Sys.glob(file.path(WORK, "*.pca", "*.projections"))[1]
SAMPLES_ORDER <- file.path(INT, "samples/samples_order.list")
CLIMATE_VALID <- file.path(INT, "samples/coord_valid_samples.list")
OUTDIR        <- file.path(ROOT, "benchmarks/mvp_eval/rda_sweep", PROJECT)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(LFMM_FULL, CLIMATE, VCFSNP_FULL, PCA, SAMPLES_ORDER, CLIMATE_VALID)) {
  if (is.na(f) || !file.exists(f)) stop("missing input: ", f)
}

PREDICTORS <- c("bio_1", "bio_2")
message("INFO: project ", PROJECT, "  K(rdadapt) = ", K_FIXED)

# ---- data, replicating scripts/rda.R sections 1-5 --------------------------------------
climate_dt <- fread(CLIMATE, sep = "\t", header = TRUE)
Env <- climate_dt[, ..PREDICTORS]
n_samples <- nrow(Env)

Y_full <- as.matrix(fread(LFMM_FULL, sep = " ", header = FALSE))
stopifnot(nrow(Y_full) == n_samples, max(Y_full, na.rm = TRUE) <= 2)
col_sds  <- apply(Y_full, 2, sd, na.rm = TRUE)
kept_idx <- setdiff(seq_len(ncol(Y_full)), which(col_sds == 0 | is.na(col_sds)))
Y <- Y_full[, kept_idx, drop = FALSE]
m_fitted <- ncol(Y_full)
message("INFO: ", n_samples, " samples x ", ncol(Y), " markers (",
        m_fitted - ncol(Y), " zero-variance dropped)")

vcfsnp <- fread(VCFSNP_FULL, header = FALSE)
stopifnot(nrow(vcfsnp) == m_fitted)
snp_dt <- data.table(chr = as.character(vcfsnp$V1), pos = as.integer(vcfsnp$V2))
snp_dt[, SNPID := paste0(chr, ":", pos)]

climate_valid_ids <- fread(CLIMATE_VALID, header = FALSE, colClasses = "character")$V2

# ---- one cell -------------------------------------------------------------------------
fit_cell <- function(scale_flag, n_cond, use_pmax) {
  Variables <- as.data.frame(Env)
  rhs <- PREDICTORS
  if (n_cond > 0) {
    cov <- load_pca_covariates(PCA, n_cond, SAMPLES_ORDER, climate_valid_ids)
    stopifnot(!is.null(cov), nrow(cov) == n_samples)
    Variables <- cbind(Variables, as.data.frame(cov))
    rhs <- c(rhs, paste0("Condition(", paste(colnames(cov), collapse = "+"), ")"))
  }
  f  <- as.formula(paste("Y ~", paste(rhs, collapse = " + ")))
  fit <- rda(f, data = Variables, scale = scale_flag)
  kmax <- fit$CCA$rank
  k    <- min(K_FIXED, kmax)
  if (k < 2) stop("rank < 2 for cell scale=", scale_flag, " cond=", n_cond)
  ra <- rdadapt(fit, k)
  p  <- ra$p.values
  gif <- ra$gif_lambda

  gif_unc <- NA_real_
  if (use_pmax && n_cond > 0) {
    fu <- rda(as.formula(paste("Y ~", paste(PREDICTORS, collapse = " + "))),
              data = Variables, scale = scale_flag)
    ku <- min(k, fu$CCA$rank)          # K pinned across fits — see knowledge note 50
    rau <- rdadapt(fu, ku)
    p <- pmax(p, rau$p.values)
    gif_unc <- rau$gif_lambda
  }
  list(p = p, k = k, gif = gif, gif_unc = gif_unc)
}

CELLS <- list(
  list(id = "A_theirs_uncorrected", scale = FALSE, cond = 0, pmax = FALSE,
       note = "reproduces Lotterhos's uncorrected arm"),
  list(id = "B_theirs_corrected",   scale = FALSE, cond = 2, pmax = FALSE,
       note = "reproduces Lotterhos's Condition(PC1+PC2) arm"),
  list(id = "C_scale_only",         scale = TRUE,  cond = 2, pmax = FALSE,
       note = "B + our scale=TRUE — isolates the scaling effect"),
  list(id = "D_ours_current",       scale = TRUE,  cond = 5, pmax = TRUE,
       note = "current CLINE-GO default (5 PCs + pmax)"),
  list(id = "E_ours_no_pmax",       scale = TRUE,  cond = 5, pmax = FALSE,
       note = "D without the intersection rule — isolates pmax"),
  list(id = "F_scaleFALSE_5pc",     scale = FALSE, cond = 5, pmax = FALSE,
       note = "5 PCs without scaling — isolates conditioning depth"),
  list(id = "G_scaleTRUE_2pc_pmax",  scale = TRUE,  cond = 2, pmax = TRUE,
       note = "pmax at the tutorial's conditioning depth — is 5 PCs needed?"),
  list(id = "H_scaleFALSE_2pc_pmax", scale = FALSE, cond = 2, pmax = TRUE,
       note = "pmax without scaling — does pmax's gain depend on scale=TRUE?")
)

rows <- list()
for (cl in CELLS) {
  message("\n=== ", cl$id, " (scale=", cl$scale, ", cond=", cl$cond, ", pmax=", cl$pmax, ") ===")
  r <- fit_cell(cl$scale, cl$cond, cl$pmax)
  full_p <- rep(NA_real_, m_fitted)
  full_p[kept_idx] <- r$p
  out <- copy(snp_dt)[, climate_multivariate := full_p]
  setorder(out, chr, pos)
  fwrite(out, file.path(OUTDIR, paste0(cl$id, "_pvalues.tsv")), sep = "\t")
  rows[[length(rows) + 1L]] <- data.table(
    cell = cl$id, scale = cl$scale, condition_pcs = cl$cond, pmax = cl$pmax,
    K = r$k, gif_lambda = round(r$gif, 4), gif_lambda_unconstrained = round(r$gif_unc, 4),
    median_p = round(median(full_p, na.rm = TRUE), 4),
    min_p = signif(min(full_p, na.rm = TRUE), 3),
    n_q05 = sum(p.adjust(full_p, "BH") < 0.05, na.rm = TRUE),
    n_bonf05 = sum(full_p < 0.05 / sum(!is.na(full_p)), na.rm = TRUE),
    note = cl$note)
  print(rows[[length(rows)]][, .(cell, gif_lambda, median_p, min_p, n_q05, n_bonf05)])
}
summ <- rbindlist(rows)
fwrite(summ, file.path(OUTDIR, "sweep_summary.tsv"), sep = "\t")
message("\nINFO: wrote ", OUTDIR, "/sweep_summary.tsv")
print(summ[, .(cell, scale, condition_pcs, pmax, gif_lambda, median_p, n_q05, n_bonf05)])
