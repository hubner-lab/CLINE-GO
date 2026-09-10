# scripts/R/lib/rdadapt.R
#
# The robust-Mahalanobis RDA candidate test (Capblancq et al. 2018), shared by the
# GEA-mode RDA scan and preGEA's RDA setup so both call SNPs with identical math.

# --- fixtures --------------------------------------------------------------

# A small but real vegan RDA: 40 samples, 60 "SNPs", 3 predictors, with the first
# five SNPs genuinely loaded on e1 so the scan has something to find.
fit_rda <- function(seed = 7L, n = 40L, p = 60L) {
    set.seed(seed)
    X <- data.frame(e1 = rnorm(n), e2 = rnorm(n), e3 = rnorm(n))
    Y <- matrix(rnorm(n * p), n, p)
    Y[, 1:5] <- Y[, 1:5] + 3 * X$e1
    vegan::rda(Y ~ e1 + e2 + e3, data = X)
}

# --- qvalue_with_fallback --------------------------------------------------

test_that("qvalue_with_fallback uses Storey when pi0 estimation succeeds", {
    set.seed(42)
    out <- qvalue_with_fallback(c(runif(90), runif(10, 0, 1e-6)))
    expect_identical(out$method, "storey")
    expect_length(out$qvalues, 100L)
})

test_that("qvalue_with_fallback drops to lambda = 0 when the pi0 spline fails", {
    # An all-tiny p-vector leaves nothing in the null tail for pi0 to be estimated
    # from, so the default lambda sequence errors and the second branch takes over.
    out <- qvalue_with_fallback(rep(1e-10, 50))
    expect_identical(out$method, "storey_lambda0")

    # Too few p-values for the spline is the other route into the same branch.
    expect_identical(qvalue_with_fallback(c(0.01, 0.5, 0.9))$method, "storey_lambda0")
})

test_that("qvalue_with_fallback falls all the way back to BH", {
    # p-values outside [0, 1] are rejected by both qvalue() calls; p.adjust() is
    # the last resort and never errors.
    out <- qvalue_with_fallback(c(-1, 0.5, 2))
    expect_identical(out$method, "BH")
    expect_equal(out$qvalues, p.adjust(c(-1, 0.5, 2), "BH"))
})

test_that("qvalue_with_fallback always returns one q-value per p-value", {
    for (p in list(runif(50), rep(1e-10, 50), c(0.01, 0.5, 0.9))) {
        expect_length(qvalue_with_fallback(p)$qvalues, length(p))
    }
})

# --- rdadapt ---------------------------------------------------------------

test_that("rdadapt returns all five documented elements, one entry per SNP", {
    mod <- fit_rda()
    r   <- rdadapt(mod, 2)
    expect_named(r, c("p.values", "q.values", "gif_lambda",
                      "qvalue_method", "distance"))
    expect_length(r$p.values, 60L)
    expect_length(r$q.values, 60L)
    expect_length(r$distance, 60L)   # returned so callers never recompute covRob
    expect_length(r$gif_lambda, 1L)
})

test_that("gif_lambda is median(distance) / qchisq(0.5, df = K)", {
    mod <- fit_rda()
    for (K in 2:3) {
        r <- rdadapt(mod, K)
        expect_equal(r$gif_lambda, median(r$distance) / qchisq(0.5, df = K))
    }
})

test_that("p-values are genuine probabilities with no NA or Inf", {
    r <- rdadapt(fit_rda(), 2)
    expect_true(all(r$p.values >= 0 & r$p.values <= 1))
    expect_false(anyNA(r$p.values))
    expect_true(all(is.finite(r$p.values)))
})

test_that("the loaded SNPs come out more extreme than the null ones", {
    # Columns 1:5 were built with a real e1 effect; this is the sanity check that
    # the statistic points the right way, not just that it runs.
    r <- rdadapt(fit_rda(), 2)
    expect_lt(median(r$p.values[1:5]), median(r$p.values[6:60]))
})

test_that("K is honoured — a different K gives a different statistic", {
    mod <- fit_rda()
    expect_false(isTRUE(all.equal(rdadapt(mod, 2)$p.values,
                                  rdadapt(mod, 3)$p.values)))
})

test_that("K = 1 errors: covRob needs at least two loading columns", {
    # Documented at rdadapt.R:40-42 — the >= 2 floor is the CALLER's job, and this
    # pins the failure mode that makes it necessary.
    expect_error(rdadapt(fit_rda(), 1), "at least two columns")
})
