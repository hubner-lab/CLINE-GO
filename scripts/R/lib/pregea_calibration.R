# pregea_calibration.R — lambda_GC, and the calibration band preGEA picks and checks rungs with
# Usage: source("/pipeline/scripts/R/lib/pregea_calibration.R")
# Requires: base R only
#
# The band is [deflation_floor, 1 + lambda_tol]: PREGEA_DEFLATION_FLOOR and PREGEA_LAMBDA_TOL in
# common.smk, passed to both scripts that use it. pregea_recommend.R takes the smallest EMMAX
# n_pcs whose LD-pruned lambda falls inside it; pregea_transfer_guard.R measures lambda again on
# the full marker set at the recommended values and says whether the calibration carried over.

# lambda_GC — the genomic inflation factor: the median 1-df chi-square statistic of the
# p-values over the null median, qchisq(0.5, 1); NA below 3 tested p-values. The same formula
# pregea_ladder_stats.R computes inline (it also needs the statistics for its genomic-control
# step). pregea_transfer_guard.R and assoc_lambda.R (the GEA/GWAS lambda tables) call this one.
lambda_gc <- function(pvals) {
    pv <- pvals[!is.na(pvals)]
    if (length(pv) < 3) return(NA_real_)
    chisq1 <- suppressWarnings(qchisq(pv, df = 1, lower.tail = FALSE))
    stats::median(chisq1, na.rm = TRUE) / qchisq(0.5, df = 1)
}

# TRUE where lambda lies inside the band; a missing lambda is outside it.
lambda_in_band <- function(lambda, deflation_floor, lambda_tol) {
    !is.na(lambda) & lambda >= deflation_floor & lambda <= 1 + lambda_tol
}

# Transfer-guard verdict for one method. "transfers" only when the full-set lambda is within
# max_delta of the pruned one AND inside the band: a small shift that leaves the band (SIMDATA
# EMMAX, 0.95 pruned -> 0.75 full) used to pass on |delta| < 0.3 alone (961f35). "unavailable"
# when either lambda is missing.
transfer_verdict <- function(lambda_pruned, lambda_full, deflation_floor, lambda_tol,
                             max_delta = 0.3) {
    if (is.na(lambda_pruned) || is.na(lambda_full)) return("unavailable")
    if (abs(lambda_full - lambda_pruned) < max_delta &&
        lambda_in_band(lambda_full, deflation_floor, lambda_tol)) "transfers" else "diverges"
}
