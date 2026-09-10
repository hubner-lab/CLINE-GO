#!/usr/bin/env Rscript
# gf_sensitivity_check.R — imputation-sensitivity check for Gradient Forest.
#
# Compares two fits of the SAME adaptive SNP set:
#   A) the imputed, individual-level model  (adaptive_model.qs)
#   B) the observed-calls-only, site-allele-frequency model (frequency_model.qs)
#
# The question it answers is not "which model is better" — it is "did the imputation decide
# the answer". On a panel with heavy missingness, agreement between a fit that invented ~44%
# of its input and a fit that invented none of it is the strongest available evidence that
# the offset is a property of the data rather than of the imputer.
#
# Writes a long key/value TSV (columns: quantity, value), matching the shape of
# geometric_offset_diagnostics.tsv so the Shiny loader convention applies unchanged.
#
# Args (positional):
#   1  MODEL_IMPUTED    — adaptive_model.qs   (individual-level, imputed)
#   2  MODEL_FREQUENCY  — frequency_model.qs  (site-level, observed only)
#   3  CLIM_PRESENT     — climate/tables/present/climate_present_site.tsv
#   4  CLIM_FUTURE      — climate/tables/future/climate_future_*_site.tsv
#   5  PREDICTORS       — comma-separated predictor names
#   6  SAMPLES          — metadata with sample + site columns
#   7  OUTPUT           — imputation_sensitivity.tsv

suppressPackageStartupMessages({
    library(data.table)
    library(gradientForest)
    library(qs)
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
MODEL_IMPUTED   <- args[1]
MODEL_FREQUENCY <- args[2]
CLIM_PRESENT    <- args[3]
CLIM_FUTURE     <- args[4]
PREDICTORS      <- trimws(str_split(args[5], ',')[[1]])
SAMPLES         <- args[6]
OUTPUT          <- args[7]

# ── Thresholds ──────────────────────────────────────────────────────────────────────────
# RHO_PASS: conventional "strong agreement" bar for a rank correlation.
# RHO_WARN: below this the agreement is not even nominally significant at the site counts
#   this pipeline sees — the two-tailed p<0.05 critical value for Spearman is ~0.62 at n=11
#   and ~0.65 at n=10, so 0.60 is the point where "they agree" stops being defensible.
#   The realised p-value is reported alongside so the reader is not stuck with a fixed bar.
# MIN_SITES / MIN_SNPS: below these the check itself is uninformative and reports not_run
#   rather than a status it cannot support.
RHO_PASS  <- 0.80
RHO_WARN  <- 0.60
TOP_K     <- 3L      # how many leading predictors must agree
MIN_SITES <- 8L
MIN_SNPS  <- 30L

out <- list()
put <- function(k, v) out[[length(out) + 1L]] <<- data.table(quantity = k, value = as.character(v))

fail_out <- function(status, reason) {
    put('status', status); put('reason', reason)
    fwrite(rbindlist(out), OUTPUT, sep = '\t')
    message('INFO: sensitivity check -> ', status, ' (', reason, ')')
    quit(save = 'no', status = 0)
}

gi <- tryCatch(qread(MODEL_IMPUTED),   error = function(e) NULL)
gf <- tryCatch(qread(MODEL_FREQUENCY), error = function(e) NULL)
if (is.null(gi) || is.null(gf)) fail_out('not_run', 'one or both Gradient Forest models could not be read')

# The frequency model writes a sentinel instead of a fitted object when too few markers
# survive observed-only aggregation (see gradient_forest_model.R). Recognise it and report
# not_run with the reason, rather than erroring on a missing $result.
if (identical(gf$status, 'insufficient')) {
    put('n_snps_frequency', gf$n_snps %||% 0L)
    fail_out('not_run', gf$reason %||% 'too few SNPs survived observed-only aggregation')
}
if (is.null(gi$result) || is.null(gf$result)) {
    fail_out('not_run', 'a model object carries no per-SNP R2 vector')
}

n_snps_imp  <- length(gi$result)
n_snps_freq <- length(gf$result)
put('n_snps_imputed', n_snps_imp)
put('n_snps_frequency', n_snps_freq)
put('mean_r2_imputed', signif(mean(gi$result), 5))
put('mean_r2_frequency', signif(mean(gf$result), 5))

if (n_snps_freq < MIN_SNPS) {
    put('min_snps_required', MIN_SNPS)
    fail_out('not_run', sprintf('only %d SNPs survived observed-only aggregation (need >=%d)',
                                n_snps_freq, MIN_SNPS))
}

# ── Site-level offsets from both models ─────────────────────────────────────────────────
samples <- fread(SAMPLES, colClasses = c(site = 'character', sample = 'character'))
ep <- fread(CLIM_PRESENT); ef <- fread(CLIM_FUTURE)
if (!'sample' %in% names(ep)) fail_out('not_run', 'present climate table has no sample column')
ep$site <- samples$site[match(ep$sample, samples$sample)]
ef$site <- samples$site[match(ef$sample, samples$sample)]

usites <- unique(samples$site)
EP <- unique(ep[, c('site', PREDICTORS), with = FALSE], by = 'site')
EF <- unique(ef[, c('site', PREDICTORS), with = FALSE], by = 'site')
EP <- EP[match(usites, EP$site)]; EF <- EF[match(usites, EF$site)]
n_sites <- nrow(EP)
put('n_sites', n_sites)

if (n_sites < MIN_SITES) {
    put('min_sites_required', MIN_SITES)
    fail_out('not_run', sprintf(
        'only %d sites (need >=%d) — a rank correlation over this few carries no information',
        n_sites, MIN_SITES))
}

offset_of <- function(model) {
    tp <- predict(model, as.data.frame(EP[, ..PREDICTORS]))
    tf <- predict(model, as.data.frame(EF[, ..PREDICTORS]))
    sqrt(rowSums((tf - tp)^2))
}
off_imp  <- tryCatch(offset_of(gi), error = function(e) NULL)
off_freq <- tryCatch(offset_of(gf), error = function(e) NULL)
if (is.null(off_imp) || is.null(off_freq)) fail_out('not_run', 'offset projection failed for one model')

ct  <- suppressWarnings(cor.test(off_imp, off_freq, method = 'spearman'))
rho <- unname(ct$estimate); pv <- ct$p.value
put('offset_spearman', signif(rho, 4))
# Exact Spearman p-values underflow to a literal 0 at these correlations; reporting "0"
# reads as a formatting bug rather than as "smaller than double precision can hold".
put('offset_spearman_p', if (!is.na(pv) && pv < 1e-10) '<1e-10' else signif(pv, 4))
put('offset_pearson', signif(suppressWarnings(cor(off_imp, off_freq)), 4))

# ── Predictor-ranking agreement ─────────────────────────────────────────────────────────
imp_i <- sort(importance(gi), decreasing = TRUE)
imp_f <- sort(importance(gf), decreasing = TRUE)
top_i <- names(imp_i)[seq_len(min(TOP_K, length(imp_i)))]
top_f <- names(imp_f)[seq_len(min(TOP_K, length(imp_f)))]
n_shared <- length(intersect(top_i, top_f))
put('top_predictors_imputed', paste(top_i, collapse = ','))
put('top_predictors_frequency', paste(top_f, collapse = ','))
put('top_predictors_shared', n_shared)
put('top_predictors_k', TOP_K)
put('top_predictor_agrees', identical(top_i[1], top_f[1]))

# ── Status ──────────────────────────────────────────────────────────────────────────────
# Both axes must hold for a pass: the same drivers AND the same site ranking. A model can
# agree on which climate variable matters while ranking the sites differently, and it is the
# site ranking that becomes the map.
status <- if (rho >= RHO_PASS && n_shared >= 2L) {
    'pass'
} else if (rho >= RHO_WARN && n_shared >= 1L) {
    'warning'
} else {
    'fail'
}
put('status', status)
put('rho_pass_threshold', RHO_PASS)
put('rho_warn_threshold', RHO_WARN)

reason <- sprintf('offset rank agreement rho=%.3f (p=%.3g), %d of top-%d predictors shared',
                  rho, pv, n_shared, TOP_K)
put('reason', reason)

fwrite(rbindlist(out), OUTPUT, sep = '\t')
message('INFO: sensitivity check -> ', status, ' | ', reason)
