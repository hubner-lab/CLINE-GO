# scripts/R/lib/pregea_calibration.R
#
# lambda_GC itself, the band [deflation_floor, 1 + lambda_tol] preGEA picks rungs in, and the
# transfer guard's verdict. pregea_transfer_guard.R itself needs LEA and EMMAX, so the
# verdict is tested here and the script's use of it statically below.

FLOOR <- 0.90
TOL   <- 0.15

test_that("lambda_in_band is inclusive at both ends and treats NA as outside", {
    expect_identical(lambda_in_band(c(0.90, 1.00, 1.15, 0.8999, 1.1501, NA), FLOOR, TOL),
                     c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("transfer_verdict passes a shift only when the full-set lambda stays in the band", {
    # 961f35, SIMDATA EMMAX: 0.95 pruned -> 0.75 full, |delta| 0.2 < 0.3, was "transfers".
    expect_identical(transfer_verdict(0.95, 0.75, FLOOR, TOL), "diverges")
    expect_identical(transfer_verdict(0.95, 1.05, FLOOR, TOL), "transfers")
    expect_identical(transfer_verdict(1.00, 1.15, FLOOR, TOL), "transfers")
    # In the band but a shift of 0.3 or more still diverges.
    expect_identical(transfer_verdict(0.60, 0.95, FLOOR, TOL), "diverges")
    expect_identical(transfer_verdict(0.95, 1.30, FLOOR, TOL), "diverges")
    expect_identical(transfer_verdict(NA, 1.00, FLOOR, TOL), "unavailable")
    expect_identical(transfer_verdict(1.00, NA, FLOOR, TOL), "unavailable")
})

test_that("pregea_transfer_guard.R takes both verdicts from transfer_verdict()", {
    src <- readLines(file.path(getOption("clinego.repo_root", "/pipeline"), "scripts",
                               "pregea_transfer_guard.R"))
    expect_identical(sum(grepl("<- transfer_verdict\\(", src)), 2L)
    expect_false(any(grepl("abs\\(delta_(lfmm|emmax)\\) < 0.3", src)))
    expect_true(any(grepl("source\\(\"/pipeline/scripts/R/lib/pregea_calibration.R\"\\)", src)))
})

test_that("lambda_gc is 1 on uniform quantiles and recovers a known inflation", {
    p <- (1:999) / 1000                       # median p = 0.5 -> the null median exactly
    expect_equal(lambda_gc(p), 1)
    # Every 1-df statistic scaled by 1.5: the median scales with it.
    inflated <- pchisq(qchisq(p, df = 1, lower.tail = FALSE) * 1.5, df = 1, lower.tail = FALSE)
    expect_equal(lambda_gc(inflated), 1.5)
    expect_equal(lambda_gc(c(p, NA, NA)), 1)  # NA p-values are left out
    expect_identical(lambda_gc(c(0.2, NA, 0.4)), NA_real_)   # fewer than 3 tested
})

test_that("lambda_gc has one definition, shared by the transfer guard and the lambda tables", {
    root <- getOption("clinego.repo_root", "/pipeline")
    guard <- readLines(file.path(root, "scripts", "pregea_transfer_guard.R"))
    tables <- readLines(file.path(root, "scripts", "assoc_lambda.R"))
    expect_false(any(grepl("^\\s*lambda_gc\\s*<-\\s*function", guard)))
    expect_true(any(grepl("lambda_gc\\(", guard)))
    expect_true(any(grepl("source\\(\"/pipeline/scripts/R/lib/pregea_calibration.R\"\\)", tables)))
})
