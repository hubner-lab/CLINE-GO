# scripts/R/utils/pval_threshold.R
#
# The contract that matters most here: THE RETURNED CUTOFF IS INCLUSIVE. For 'top'
# and 'qval' the returned value is itself a member of the intended call set, so a
# caller using `p < threshold` drops the boundary observation. sig_snps.R:67 uses
# `<=`; these tests pin the property that makes that correct.

# --- helpers ---------------------------------------------------------------

# 90 uniform + 10 near-zero: enough signal for Storey's pi0 estimator to converge.
# Plain uniform p-values often make qvalue() fail, which the code turns into
# status = "too_few_tests" — a flaky test rather than a meaningful one.
p_with_signal <- function(n_null = 90L, n_signal = 10L, seed = 42L) {
    set.seed(seed)
    c(runif(n_null), runif(n_signal, min = 0, max = 1e-6))
}

# --- bonf ------------------------------------------------------------------

test_that("bonf threshold is exactly alpha / n_tested", {
    res <- compute_pval_threshold(seq(0.01, 0.20, by = 0.01), "bonf", 0.05)
    expect_identical(res$status, "ok")
    expect_identical(res$n_tested, 20L)
    expect_equal(res$threshold, 0.05 / 20)
})

test_that("NA p-values are dropped and counted, and change the bonf denominator", {
    res <- compute_pval_threshold(c(0.1, 0.2, NA, 0.3, NA), "bonf", 0.05)
    expect_identical(res$n_tested, 3L)
    expect_identical(res$n_na_dropped, 2L)
    expect_equal(res$threshold, 0.05 / 3)   # not 0.05 / 5
})

test_that("every returned threshold lies in [0, 1]", {
    ps <- p_with_signal()
    for (args in list(list("bonf", 0.05), list("top", 10), list("custom", 0.001),
                      list("qval", 0.1))) {
        res <- compute_pval_threshold(ps, args[[1]], args[[2]])
        if (identical(res$status, "ok")) {
            expect_gte(res$threshold, 0)
            expect_lte(res$threshold, 1)
        }
    }
})

# --- qval ------------------------------------------------------------------

test_that("qval refuses below 10 tests and does NOT silently fall back to BH", {
    # Deliberate: pval_threshold.R:39-40 says no automatic fallback is applied, so
    # the caller can surface an actionable warning instead of a quietly different
    # multiple-testing correction. compute_qvalues_safe() is the one that falls back.
    res <- expect_message(
        compute_pval_threshold(c(0.001, 0.002, 0.5, 0.6, 0.7), "qval", 0.1),
        "requires >=10 tests"
    )
    expect_identical(res$status, "too_few_tests")
    expect_true(is.na(res$threshold))
    expect_identical(res$n_tested, 5L)
})

test_that("a successful qval threshold is itself one of the input p-values", {
    ps  <- p_with_signal()
    res <- compute_pval_threshold(ps, "qval", 0.1)
    skip_if_not(identical(res$status, "ok"), "qvalue() did not converge on this fixture")
    # max_pvalue_fdr() returns the LARGEST p whose q < fdr, so the cutoff is a
    # member of the call set — which is why selection must be `p <= threshold`.
    expect_true(any(ps == res$threshold))
    expect_gte(sum(ps <= res$threshold), 1L)
})

test_that("max_pvalue_fdr returns NA when nothing passes and when all input is NA", {
    expect_true(is.na(max_pvalue_fdr(rep(NA_real_, 5), 0.05)))
    expect_true(is.na(max_pvalue_fdr(numeric(0), 0.05)))
})

# --- top -------------------------------------------------------------------

test_that("top threshold is the Nth smallest p-value", {
    ps  <- c(0.5, 0.01, 0.3, 0.02, 0.4)
    res <- compute_pval_threshold(ps, "top", 3)
    expect_identical(res$status, "ok")
    expect_equal(res$threshold, 0.3)          # sorted: .01 .02 .3 | .4 .5
})

test_that("top selection with `p <= threshold` returns exactly N when there are no ties", {
    ps  <- c(0.5, 0.01, 0.3, 0.02, 0.4)
    res <- compute_pval_threshold(ps, "top", 3)
    expect_identical(sum(ps <= res$threshold), 3L)
})

test_that("top with ties at the Nth p returns MORE than N — documented as correct", {
    ps  <- c(0.01, 0.02, 0.30, 0.30, 0.40)
    res <- compute_pval_threshold(ps, "top", 3)
    expect_equal(res$threshold, 0.30)
    expect_identical(sum(ps <= res$threshold), 4L)
})

test_that("top over-request returns a status, while bare max_pvalue_top() errors", {
    res <- expect_message(compute_pval_threshold(c(0.1, 0.2), "top", 10), "top N=10")
    expect_identical(res$status, "too_few_tests")
    expect_true(is.na(res$threshold))

    expect_error(max_pvalue_top(c(0.1, 0.2), 10), "larger than the number of non-NA")
})

test_that("max_pvalue_top ignores NA when counting against topN", {
    expect_equal(max_pvalue_top(c(0.1, NA, 0.2, NA, 0.3), 2), 0.2)
})

# --- custom ----------------------------------------------------------------

test_that("custom passes a positive number through unchanged", {
    res <- compute_pval_threshold(c(0.1, 0.2, 0.3), "custom", 1e-4)
    expect_identical(res$status, "ok")
    expect_equal(res$threshold, 1e-4)
})

test_that("custom rejects non-numeric and non-positive values", {
    for (bad in list("abc", 0, -1)) {
        res <- suppressMessages(compute_pval_threshold(c(0.1, 0.2), "custom", bad))
        expect_identical(res$status, "too_few_tests")
        expect_true(is.na(res$threshold))
    }
})

# --- dispatch edges --------------------------------------------------------

test_that("an all-NA / empty input short-circuits to no_tests before dispatch", {
    res <- compute_pval_threshold(rep(NA_real_, 4), "bonf", 0.05)
    expect_identical(res$status, "no_tests")
    expect_identical(res$n_tested, 0L)
    expect_identical(res$n_na_dropped, 4L)
    expect_true(is.na(res$threshold))

    # Reached before the adjustment is even looked at, so an unknown method here
    # returns rather than erroring.
    expect_identical(compute_pval_threshold(numeric(0), "nonsense", 1)$status, "no_tests")
})

test_that("an unknown adjustment method errors", {
    expect_error(compute_pval_threshold(c(0.1, 0.2), "fdr_bh", 0.05),
                 "Unknown adjustment method")
})

# --- compute_qvalues_safe --------------------------------------------------

test_that("compute_qvalues_safe preserves length and stays in [0, 1]", {
    ps <- p_with_signal()
    qs <- compute_qvalues_safe(ps)
    expect_length(qs, length(ps))
    expect_true(all(qs >= 0 & qs <= 1))
})

test_that("compute_qvalues_safe DOES fall back to BH when qvalue() fails", {
    # Two p-values: far too few for pi0 estimation, so qvalue() throws and the
    # tryCatch takes over. This is the fallback the dossier's Tier 1 bullet meant
    # — compute_pval_threshold() deliberately has no equivalent.
    ps <- c(0.001, 0.5)
    qs <- expect_message(compute_qvalues_safe(ps), "using BH adjustment")
    expect_equal(qs, p.adjust(ps, method = "BH"))
})
