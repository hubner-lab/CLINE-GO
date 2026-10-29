# scripts/R/lib/pregea_calibration.R
#
# The lambda_GC band [deflation_floor, 1 + lambda_tol] preGEA picks rungs in, and the
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
