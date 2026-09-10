test_that("chr_midpoints computes correctly", {
    coords <- list(
        chr_offsets = list("1" = 0,     "2" = 15000, "3" = 32000),
        chr_lengths = list("1" = 15000, "2" = 17000, "3" = 18000)
    )
    mids <- chr_midpoints(coords)
    expect_equal(mids, c(7500, 23500, 41000))
})

test_that("add_cum_pos adds cumulative positions", {
    coords <- list(
        chr_offsets = list("1" = 0, "2" = 100000),
        chr_lengths = list("1" = 100000, "2" = 80000),
        x_range = c(0, 180000),
        y_range = c(0, 15)
    )
    snps <- data.table::data.table(
        SNPID = c("s1", "s2"),
        chr   = c("1", "2"),
        pos   = c(50000L, 10000L),
        pvalue = c(1e-8, 1e-5)
    )
    result <- add_cum_pos(snps, coords)
    expect_equal(result$cum_pos, c(50000, 110000))
    expect_true("log10p" %in% names(result))
})

test_that("build_manhattan_plotly returns a plotly object", {
    coords <- list(
        chr_offsets = list("1" = 0),
        chr_lengths = list("1" = 100000),
        x_range = c(0, 101000),
        y_range = c(0, 15),
        bonferroni_y = 7.3
    )
    p <- build_manhattan_plotly(bg_uri = NULL, coords = coords)
    expect_s3_class(p, "plotly")
})
