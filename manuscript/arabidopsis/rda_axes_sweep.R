#!/usr/bin/env Rscript
# rda_axes_sweep.R: rdadapt candidates by number of retained axes, on the RDA that pregea fits.
#
# rdadapt's count depends on the number of axes K, and the pipeline's ladder holds one count per rung
# at K = the significant axes, plus one plain reference at K = 2. On the whole genome every axis sits
# at the permutation floor, so that K is always the predictor count and no table shows what K does.
#
# Reads the same inputs as pregea_rda_setup.R and fits the same model through the image's rdadapt.R
# and load_pca_covariates(). The fit needs no permutations. Stops unless K = rank reproduces the
# ladder's n_candidates_partial on every swept rung and K = 2 reproduces n_candidates_plain.
#
# Writes under ANALYSIS/Arabidopsis1135/pregea/:
#   rda_axes_sweep.tsv  condition_pcs K n_markers gif_lambda n_bonf n_q05 n_q10 n_bonf_shared_with_K2
#   rda_axis_eig.tsv    condition_pcs axis eigenvalue share broken_stick above_broken_stick
#
# With the argument `full` it fits production's model instead: every imputed SNP at condition_pcs 0.
# That matrix only exists after mode=gea and the ladder has no row for it, so the check is production
# itself: K = 2 must reproduce every p-value of RDA_pvalues_K6.tsv. The pruned tables above are the
# evidence behind the axes choice and stay untouched; this writes rda_axes_sweep_full.tsv,
# rda_axis_eig_full.tsv and rda_axes_full_hits.tsv (K SNPID chr pos pvalue). The full sweep table adds
# min_p, n_top_shared_with_K2 (of each K's top-ranked SNPs, as many as K=2 calls, how many are K=2
# calls) and n_q05_shared_with_K2.
#
# broken_stick is the share each axis would hold if the constrained variance split at random, the
# null the pipeline already draws on the genomic PCA screeplot.
#
# Run from the repo root, in the image:
#   docker run --rm --name arab-rda-axes --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp \
#     -e OPENBLAS_NUM_THREADS=1 -v $PWD:/data -w /data cline-go:6cdedcc \
#     Rscript manuscript/arabidopsis/rda_axes_sweep.R [full]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)       # load_pca_covariates() uses %>%
  library(here)
  library(yaml)
  library(vegan)
  library(robust)
  library(qvalue)
})
source("/pipeline/scripts/R/utils/emmax_core.R")  # load_pca_covariates()
source("/pipeline/scripts/R/lib/rdadapt.R")

project <- "Arabidopsis1135"
full_scope <- identical(commandArgs(trailingOnly = TRUE), "full")
rungs <- if (full_scope) 0L else c(0L, 3L, 6L)

results <- here(paste0(project, "_results"))
config <- read_yaml(here(paste0("config_", project, ".yaml")))
predictors <- strsplit(config$PreGEA$predictors, ",", fixed = TRUE)[[1]]
out_dir <- here("ANALYSIS", project, "pregea")

#-------------------------------------------------------------------------- inputs
pca <- Sys.glob(file.path(results, "_work", "*", "ld*", "*.pca", "*.projections"))
stopifnot(length(pca) == 1L)
samples_file <- file.path(results, "_intermediate", "samples", "samples_order.list")
samples <- readLines(samples_file)

climate <- fread(file.path(results, "climate", "tables", "present", "climate_present_site_scaled.tsv"),
                 colClasses = c(sample = "character"))
stopifnot(identical(climate$sample, samples), all(predictors %in% names(climate)))
Env <- as.data.frame(climate[, ..predictors])

matrix_name <- if (full_scope) "lfmm_imp_full_climate.lfmm" else "lfmm_imp_climate.lfmm"
Y <- as.matrix(fread(file.path(results, "_intermediate", "climate_subset", matrix_name),
                     sep = " ", header = FALSE))
stopifnot(nrow(Y) == nrow(Env))
varies <- apply(Y, 2, sd) > 0
Y <- Y[, varies, drop = FALSE]
m <- ncol(Y)

if (full_scope) {
  # Production's SNP order: the imputed filtered set, one line per matrix column.
  vcfsnp <- Sys.glob(file.path(results, "_work", "*", "*_K6imp_maf*.vcfsnp"))
  stopifnot(length(vcfsnp) == 1L)
  snps <- fread(vcfsnp, header = FALSE, select = 1:2, col.names = c("chr", "pos"),
                colClasses = list(character = 1))
  stopifnot(nrow(snps) == length(varies))
  snps <- snps[varies][, SNPID := paste0(chr, ":", pos)]
  production <- fread(file.path(results, "GEA", "tables", "methods", "RDA", "RDA_pvalues_K6.tsv"),
                      colClasses = c(chr = "character", SNPID = "character"))
  production <- production[match(snps$SNPID, SNPID)]
  stopifnot(!anyNA(production$climate_multivariate))
} else {
  ladder <- fread(file.path(results, "PreGEA", "tables", "rda", "rda_condition_ladder.tsv"))
}

#-------------------------------------------------------------------------- sweep
n_bonf_of <- function(p) sum(p < 0.05 / m, na.rm = TRUE)

sweep_rung <- function(cond_pcs) {
  covariates <- load_pca_covariates(pca, cond_pcs, samples_file, samples)
  variables <- Env
  rhs <- predictors
  if (!is.null(covariates)) {
    variables <- cbind(variables, as.data.frame(covariates))
    rhs <- c(rhs, paste0("Condition(", paste(colnames(covariates), collapse = "+"), ")"))
  }
  fit <- rda(as.formula(paste("Y ~", paste(rhs, collapse = " + "))), data = variables, scale = TRUE)

  K_max <- fit$CCA$rank
  eig <- as.numeric(fit$CCA$eig)
  stick <- rev(cumsum(1 / rev(seq_len(K_max)))) / K_max
  stopifnot(isTRUE(all.equal(sum(stick), 1)))
  eig_table <- data.table(condition_pcs = cond_pcs, axis = seq_len(K_max), eigenvalue = eig,
                          share = eig / sum(eig), broken_stick = stick)

  calls <- lapply(2:K_max, function(K) rdadapt(fit, K))
  names(calls) <- 2:K_max
  hits_K2 <- which(calls[["2"]]$p.values < 0.05 / m)
  sweep_table <- rbindlist(lapply(2:K_max, function(K) {
    call <- calls[[as.character(K)]]
    data.table(condition_pcs = cond_pcs, K = K, n_markers = m, gif_lambda = call$gif_lambda,
               n_bonf = n_bonf_of(call$p.values),
               n_q05 = sum(call$q.values < 0.05, na.rm = TRUE),
               n_q10 = sum(call$q.values < 0.10, na.rm = TRUE),
               n_bonf_shared_with_K2 = sum(which(call$p.values < 0.05 / m) %in% hits_K2))
  }))

  if (full_scope) {
    difference <- all.equal(calls[["2"]]$p.values, production$climate_multivariate, tolerance = 1e-6)
    if (!isTRUE(difference)) stop("K=2 does not reproduce RDA_pvalues_K6.tsv: ", difference)
    message("INFO: K=2 reproduces every production p-value, gif ", signif(calls[["2"]]$gif_lambda, 4))
    # A K with no call still has a ranking: does it rank K=2's calls on top, or other SNPs?
    q05_K2 <- calls[["2"]]$q.values < 0.05
    sweep_table[, `:=`(
      min_p = vapply(K, function(k) min(calls[[as.character(k)]]$p.values), numeric(1)),
      n_top_shared_with_K2 = vapply(K, function(k) {
        sum(order(calls[[as.character(k)]]$p.values)[seq_along(hits_K2)] %in% hits_K2)
      }, integer(1)),
      n_q05_shared_with_K2 = vapply(K, function(k) {
        sum(calls[[as.character(k)]]$q.values < 0.05 & q05_K2)
      }, integer(1))
    )]
    hits <- rbindlist(lapply(2:K_max, function(K) {
      called <- which(calls[[as.character(K)]]$p.values < 0.05 / m)
      snps[called, .(K = K, SNPID, chr, pos, pvalue = calls[[as.character(K)]]$p.values[called])]
    }))
    return(list(sweep = sweep_table, eig = eig_table, hits = hits))
  }

  expected <- ladder[condition_pcs == cond_pcs, n_candidates_partial]
  observed <- sweep_table[K == K_max, n_bonf]
  if (observed != expected) stop("rung ", cond_pcs, ": K=", K_max, " gives ", observed,
                                 ", the ladder has ", expected)
  if (cond_pcs == 0L) {
    expected <- ladder[condition_pcs == 0L, n_candidates_plain]
    observed <- sweep_table[K == 2L, n_bonf]
    if (observed != expected) stop("plain reference: K=2 gives ", observed, ", the ladder has ", expected)
  }
  message("INFO: rung ", cond_pcs, " reproduces the ladder")
  list(sweep = sweep_table, eig = eig_table)
}

parts <- lapply(rungs, sweep_rung)
sweep <- rbindlist(lapply(parts, `[[`, "sweep"))
eig <- rbindlist(lapply(parts, `[[`, "eig"))
eig[, above_broken_stick := share > broken_stick]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- if (full_scope) "_full" else ""
fwrite(sweep[, gif_lambda := signif(gif_lambda, 4)],
       file.path(out_dir, paste0("rda_axes_sweep", suffix, ".tsv")), sep = "\t")
numeric_cols <- c("eigenvalue", "share", "broken_stick")
fwrite(eig[, (numeric_cols) := lapply(.SD, signif, 4), .SDcols = numeric_cols],
       file.path(out_dir, paste0("rda_axis_eig", suffix, ".tsv")), sep = "\t")
if (full_scope) {
  fwrite(rbindlist(lapply(parts, `[[`, "hits")), file.path(out_dir, "rda_axes_full_hits.tsv"), sep = "\t")
}
message("Wrote ", file.path(out_dir, paste0("rda_axes_sweep", suffix, ".tsv")), " and its siblings")
