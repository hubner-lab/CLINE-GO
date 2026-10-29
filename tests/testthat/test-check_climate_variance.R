# Invariant predictors and traits — CONTENT assertions, quick tier.
#
# check_climate_variance.R's wrapper rows assert exit 0 and a non-empty table, which
# an empty table passes. These read the rows:
#   * 932a8d: the climate call names its predictors. It used to grep "^bio_", which
#     matched nothing for custom climate, and a header-only table read as "every
#     predictor varies". A named predictor missing from the table is an error.
#   * 07991d: a column with ONE non-NA value is single_value. sd() of it is NA, and
#     `else if (NA)` aborted mode=climate and mode=traits.
#   * e83e54: the climate summary counts the invariant predictors among the selected
#     ones; the table itself keeps every predictor the app offers as a chip.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

fx_site_mixed <- function(d) {
    write_tsv(data.table::data.table(
        sample = sprintf("ID%03d", 1:4),
        bio_2  = c(4, 4, 4, 4),                  # zero variance
        bio_1  = c(1.5, 2.5, 3.5, 4.5),          # varies
        bio_3  = c(7, NA, NA, NA),               # one value
        bio_4  = c(NA, NA, NA, NA)),             # all NA
        file.path(d, "climate_present_site.tsv"))
}

test_that("named predictors are screened in table order with their reasons", {
    d <- withr::local_tempdir()
    out <- file.path(d, "invariant.tsv")
    res <- run_wrapper("check_climate_variance.R",
                       c(fx_site_mixed(d), out, "bio_1,bio_3,bio_4,bio_2"))
    expect_identical(res$status, 0L, info = res$output)
    got <- data.table::fread(out)
    expect_identical(got$predictor, c("bio_2", "bio_3", "bio_4"))
    expect_identical(got$reason, c("zero_variance", "single_value", "all_na"))
})

test_that("a named predictor the table lacks stops the check", {
    d <- withr::local_tempdir()
    res <- run_wrapper("check_climate_variance.R",
                       c(fx_site_mixed(d), file.path(d, "invariant.tsv"), "bio_1,tmin_jan"))
    expect_false(res$status == 0L)
    expect_match(res$output, "not in .*tmin_jan")
})

test_that("a trait with one non-NA value is single_value, not an abort", {
    d <- withr::local_tempdir()
    meta <- write_tsv(data.table::data.table(
        site = c("A", "A", "B"), sample = c("s1", "s2", "s3"),
        latitude = c(30, 30, 31), longitude = c(34, 34, 35),
        height = c(1.2, 1.4, 1.1), score = c(3, NA, NA)),
        file.path(d, "metadata.tsv"))
    out <- file.path(d, "invariant_traits.tsv")
    res <- run_wrapper("check_climate_variance.R", c(meta, out, "traits"))
    expect_identical(res$status, 0L, info = res$output)
    got <- data.table::fread(out)
    expect_identical(got$predictor, "score")
    expect_identical(got$reason, "single_value")
})

test_that("the climate summary counts invariant predictors among the selected ones", {
    d <- withr::local_tempdir()
    inv <- write_tsv(data.table::data.table(predictor = c("bio_12", "bio_14"),
                                            reason = "zero_variance"),
                     file.path(d, "climate_invariant_predictors.tsv"))
    count <- function(predictors) {
        out <- file.path(d, paste0("summary_", gsub(",", "_", predictors), ".tsv"))
        res <- run_wrapper("write_summary.R",
                           c("climate", out, inv, "NULL", "NULL", "NULL", "NULL", predictors))
        expect_identical(res$status, 0L, info = res$output)
        s <- data.table::fread(out, colClasses = "character")
        s[step == "climate" & metric == "n_invariant_predictors", value]
    }
    expect_identical(count("bio_1,bio_2,bio_3"), "0")   # SIMDATA: bio_14 is not selected
    expect_identical(count("bio_1,bio_12"), "1")
})
