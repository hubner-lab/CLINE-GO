#!/usr/bin/env Rscript
# pregea_ladder_stats.R — per-rung calibration diagnostics for ANY preGEA
# ladder engine (LFMM K sweep, EMMAX #PC sweep). One script, two engines —
# the wide p-value table shape (SNPID/chr/pos/<trait...>) is identical, so
# the diagnostics computation is engine-agnostic; only the genomic-control
# handling differs (see APPLY_GC below).
#
# Nothing else in the pipeline computes lambda_GC today
# (docs/methods_review_suboptimal.md S2) — this is the first place it exists.

suppressPackageStartupMessages({
    library(data.table)
    library(qvalue)
})

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/io_pvalues.R")

args <- commandArgs(trailingOnly = TRUE)
################################################################################
PVALUES_TSV  <- args[1]
ENGINE       <- args[2]                      # "lfmm" | "emmax"
RUNG_PARAM   <- args[3]                      # "K" | "n_pcs"
RUNG_VALUE   <- args[4]
FDR          <- as.numeric(args[5])
BONF_ALPHA   <- as.numeric(args[6])
# APPLY_GC: LFMM's rung table is fit with genomic.control=FALSE (preGEA arg 8
# to lfmm.R) precisely so the RAW p-values here are informative for
# lambda/histogram diagnosis. To make the reported hit counts reflect what the
# production genomic.control=TRUE GEA run would actually produce, GC
# recalibration is applied HERE, downstream of the raw diagnostics, exactly
# mirroring LEA's own lfmm2.test(genomic.control=TRUE) formula. EMMAX passes
# FALSE — EMMAX has no equivalent recalibration step in the production run.
APPLY_GC     <- toupper(args[7]) == "TRUE"
OUT_STATS_TSV<- args[8]
################################################################################

dir.create(dirname(OUT_STATS_TSV), recursive = TRUE, showWarnings = FALSE)

message("INFO: preGEA ladder stats — engine=", ENGINE, " ", RUNG_PARAM, "=", RUNG_VALUE,
       " apply_gc=", APPLY_GC)

pvals_dt <- read_pvalues_tsv(PVALUES_TSV)
fixed_cols <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
trait_cols <- setdiff(names(pvals_dt), fixed_cols)

# hist_shape classification thresholds — a pragmatic diagnostic heuristic
# (not a formal test): compares 20-bin histogram density in the first bin
# (near p=0), last bin (near p=1), and the mean of the interior bins against
# the uniform-null expectation (ratio 1.0 = matches null exactly).
classify_shape <- function(first_ratio, last_ratio, mid_ratio) {
    if (first_ratio > 1.3 && last_ratio < 1.3 && mid_ratio > 0.7 && mid_ratio < 1.3) {
        "flat_spike"     # good: excess near 0 (signal), flat elsewhere (well-calibrated null)
    } else if (first_ratio > 1.3 && last_ratio > 1.3) {
        "u_shape"        # excess at both tails: anti-conservative / miscalibrated
    } else if (mid_ratio > 1.3 && first_ratio < 1.3) {
        "hump"           # excess in the middle: under-dispersion
    } else if (first_ratio < 0.7) {
        "depleted"       # no spike near 0: over-correction (K/PCs too high)
    } else {
        "flat_spike"
    }
}

stats_one_trait <- function(pvec, trait_name) {
    n_total  <- length(pvec)
    pvec_nona <- pvec[!is.na(pvec)]
    n_tests  <- length(pvec_nona)

    if (n_tests < 3) {
        return(data.table(trait = trait_name,
                          metric = c("n_tests", "hist_shape"),
                          value  = c(as.character(n_tests), "degenerate")))
    }

    # lambda_GC — the standard genomic-inflation-factor formula: median chi-sq
    # (1 df) statistic over the median of the null chi-sq(1) distribution.
    chisq1 <- suppressWarnings(qchisq(pvec_nona, df = 1, lower.tail = FALSE))
    lambda <- stats::median(chisq1, na.rm = TRUE) / qchisq(0.5, df = 1)

    # LEA's lfmm2.test(genomic.control = TRUE) divides the chi-sq statistics by
    # the GIF WHATEVER it is, with no floor at 1 — and so does this project's own
    # reimplementation of the same quantity in benchmarks/lib_detection.R:77,
    # which is what scores the pipeline against ground truth. A max(lambda, 1)
    # floor used to sit here; it silently turned APPLY_GC into a no-op on every
    # DEFLATED rung (lambda < 1), which is exactly the regime this project lives
    # in (EMMAX deflation under high selfing; over-corrected high-K LFMM rungs,
    # the ones classify_shape() calls "depleted"). Since production LFMM runs
    # with genomic.control = TRUE, the floored hit counts were NOT "what the
    # production run would produce" — the stated purpose of applying GC here.
    # Measured on a synthetic deflated null at lambda = 0.90, Bonferroni 0.05,
    # 200k tests: 79 hits floored vs 100 unfloored, a 21% understatement of the
    # metric the ladder plots and pregea_recommend.R read.
    if (APPLY_GC) {
        lambda_use <- if (is.finite(lambda) && lambda > 0) lambda else 1
        p_for_hits <- pchisq(chisq1 / lambda_use, df = 1, lower.tail = FALSE)
    } else {
        lambda_use <- NA_real_
        p_for_hits <- pvec_nona
    }

    bonf_res <- tryCatch(compute_pval_threshold(p_for_hits, "bonf", BONF_ALPHA),
                         error = function(e) list(status = "error"))
    qval_res <- tryCatch(compute_pval_threshold(p_for_hits, "qval", FDR),
                         error = function(e) list(status = "error"))
    # Inclusive '<=' — compute_pval_threshold() returns a cutoff that is itself a
    # member of the call set under 'qval' (see its header). A strict '<' made
    # hits_qval one short (plus ties) on every rung, which biased the rung
    # comparison the preGEA recommender reads.
    hits_bonf <- if (identical(bonf_res$status, "ok")) sum(p_for_hits <= bonf_res$threshold, na.rm = TRUE) else NA_integer_
    hits_qval <- if (identical(qval_res$status, "ok")) sum(p_for_hits <= qval_res$threshold, na.rm = TRUE) else NA_integer_

    # Histogram shape — computed on RAW p (not GC-corrected): the point is
    # diagnosing the rung's null behavior, which GC-correction would mask.
    ks <- suppressWarnings(stats::ks.test(pvec_nona, "punif"))
    hist_flatness_ks <- unname(ks$statistic)

    n_bins <- 20
    breaks <- seq(0, 1, length.out = n_bins + 1)
    counts <- hist(pvec_nona, breaks = breaks, plot = FALSE)$counts
    expected_per_bin <- n_tests / n_bins
    first_ratio <- counts[1] / expected_per_bin
    last_ratio  <- counts[n_bins] / expected_per_bin
    mid_ratio   <- mean(counts[2:(n_bins - 1)]) / expected_per_bin
    frac_p_gt_half <- mean(pvec_nona > 0.5)
    hist_shape <- classify_shape(first_ratio, last_ratio, mid_ratio)

    data.table(
        trait  = trait_name,
        metric = c("n_tests", "lambda_gc", "lambda_gc_applied",
                  "hits_bonf", "hits_qval",
                  "hist_flatness_ks", "hist_spike0", "frac_p_gt_half", "hist_shape"),
        # lambda_gc      = the rung's raw inflation factor (diagnostic)
        # lambda_gc_applied = the divisor actually used for hits_* (NA when
        #                   APPLY_GC is FALSE) — without it the hit columns
        #                   cannot be audited against the production run.
        value  = c(as.character(n_tests), as.character(lambda),
                  as.character(lambda_use),
                  as.character(hits_bonf), as.character(hits_qval),
                  as.character(hist_flatness_ks), as.character(first_ratio > 1.3),
                  as.character(frac_p_gt_half), hist_shape)
    )
}

per_trait <- rbindlist(lapply(trait_cols, function(tr) stats_one_trait(pvals_dt[[tr]], tr)))

# Pooled row — all traits' p-values stacked into one vector. Useful for a
# single ladder-summary number per rung when there are many traits.
pooled_pvec <- unlist(lapply(trait_cols, function(tr) pvals_dt[[tr]]))
pooled <- stats_one_trait(pooled_pvec, "__pooled__")

out <- rbind(per_trait, pooled)
out[, `:=`(engine = ENGINE, rung_param = RUNG_PARAM, rung_value = RUNG_VALUE)]
setcolorder(out, c("engine", "rung_param", "rung_value", "trait", "metric", "value"))

fwrite(out, OUT_STATS_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote ladder rung stats: ", OUT_STATS_TSV, " (", nrow(out), " rows)")
