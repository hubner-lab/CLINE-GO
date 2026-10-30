# scripts/R/lib/factor_correlation.R
#
# The matrix every correlogram draws and every correlation table lists (bffcfc): two blocks
# joined on `sample`, constant factors dropped, site-level factors correlated over one row
# per site. plot_correlation_heatmap.R and correlation_table.R both call it.

write_block <- function(dt) {
    path <- withr::local_tempfile(fileext = ".tsv", .local_envir = parent.frame())
    data.table::fwrite(dt, path, sep = "\t")
    path
}

# Climate-like table keyed by sample only: 4 sites, very unbalanced (5/1/1/1 samples), so a
# sample-weighted r differs from the site-level one.
site_rows <- data.table::data.table(bio_1 = c(10, 12, 15, 20), bio_2 = c(3, 1, 4, 9),
                                    bio_3 = c(7, 7, 7, 7))
n_per <- c(5L, 1L, 1L, 1L)
climate <- cbind(data.table::data.table(sample = sprintf("s%02d", 1:8)),
                 site_rows[rep(seq_len(4), n_per)])

test_that("site-level factors are correlated over one row per site", {
    fc <- quiet(factor_correlation(write_block(climate)))
    expect_equal(fc$r["bio_1", "bio_2"], cor(site_rows$bio_1, site_rows$bio_2))
    expect_false(isTRUE(all.equal(fc$r["bio_1", "bio_2"], cor(climate$bio_1, climate$bio_2))))
    expect_identical(fc$basis["bio_1", "bio_2"], "sites")
    expect_identical(fc$n_obs["bio_1", "bio_2"], 4)
})

test_that("constant factors are dropped and named", {
    fc <- quiet(factor_correlation(write_block(climate)))
    expect_identical(colnames(fc$r), c("bio_1", "bio_2"))
    expect_identical(fc$dropped, "bio_3")
})

test_that("a per-sample factor joined on sample is correlated over samples", {
    traits <- data.table::data.table(site = rep(c("A", "B", "C", "D"), n_per),
                                     sample = rev(climate$sample),      # another row order
                                     height = c(1, 4, 2, 8, 5, 7, 3, 6))
    fc <- quiet(factor_correlation(write_block(traits), write_block(climate)))
    joined <- merge(traits, climate, by = "sample")
    expect_equal(fc$r["height", "bio_1"], cor(joined$height, joined$bio_1))
    expect_identical(fc$basis["height", "bio_1"], "samples")
    expect_identical(fc$n_obs["height", "bio_1"], 8)
})

test_that("blocks without a shared key and different row counts are refused", {
    a <- data.table::data.table(x = 1:4, y = c(2, 1, 4, 3))
    b <- data.table::data.table(z = 1:3)
    expect_error(quiet(factor_correlation(write_block(a), write_block(b))),
                 "refusing to bind them positionally")
})

test_that("correlation_pairs lists each pair once, in the matrix's order", {
    dt <- data.table::data.table(sample = sprintf("s%d", 1:6), a = c(1, 2, 3, 4, 5, 7),
                                 b = c(2, 1, 4, 3, 6, 5), c = c(9, 7, 8, 5, 6, 4))
    fc <- quiet(factor_correlation(write_block(dt)))
    p <- correlation_pairs(fc)
    expect_identical(names(p), c("factor_1", "factor_2", "r", "basis", "n_obs"))
    expect_identical(paste(p$factor_1, p$factor_2), c("a b", "a c", "b c"))
    expect_equal(p$r, c(cor(dt$a, dt$b), cor(dt$a, dt$c), cor(dt$b, dt$c)))
    expect_identical(p$n_obs, rep(6L, 3))
})

test_that("the correlogram and the correlation tables take the matrix from one function", {
    root <- getOption("clinego.repo_root", "/pipeline")
    for (f in c("plot_correlation_heatmap.R", "correlation_table.R")) {
        src <- readLines(file.path(root, "scripts", f))
        expect_true(any(grepl("factor_correlation\\(TABLE_A, TABLE_B\\)", src)), info = f)
        expect_false(any(grepl("\\bcor\\(", src)), info = f)
    }
})
