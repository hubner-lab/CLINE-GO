# scripts/R/utils/manhattan_utils.R — the coordinate arithmetic every Manhattan,
# Miami and combined plot is drawn from, plus the palettes.
#
# prepare_manhattan_data() is the one that matters: its chr_info output is written
# verbatim into the *_coords.json the Shiny app reads (plot_manhattan.R:246-249),
# so a change here silently moves every interactive overlay marker. The app
# reimplements the consumer side of that arithmetic independently — asserted in
# tests/testthat/test-equivalence-app-pipeline.R, not here.
#
# parse_assoc_files_str() is NOT covered here: it already has four blocks at
# test-io.R:90-115.
#
# The zero-length palette defects (get_trait_colors(character(0)) and friends) are
# quarantined as correct-behaviour assertions in test-known-bugs.R, so this file
# deliberately never calls them with empty input.

mh_df <- function(chr = c("1", "1", "2", "2"),
                  pos = c(100L, 500L, 200L, 800L),
                  pvalue = c(0.1, 0.01, 0.001, 1e-4)) {
    data.frame(chr = chr, pos = pos, pvalue = pvalue, stringsAsFactors = FALSE)
}

test_that("prepare_manhattan_data orders chromosomes by their embedded digits", {
    # Input order is deliberately scrambled, and "10" must land after "2" — a
    # lexical sort would put it between "1" and "2".
    r <- prepare_manhattan_data(
        mh_df(chr = c("10", "2", "1"), pos = c(100L, 200L, 300L),
              pvalue = c(0.1, 0.2, 0.3)),
        pval_col = "pvalue")
    expect_identical(levels(r$data$chr_f), c("1", "2", "10"))
    expect_identical(as.character(r$chr_info$chr_f), c("1", "2", "10"))
})

test_that("prepare_manhattan_data puts a chromosome name with no digits LAST", {
    # str_extract(., "\\d+") is NA for "scaffoldA", and order() sends NA to the
    # end by default. Worth pinning: a scaffold-named assembly is ordinary, and
    # this is the behaviour the coords JSON inherits.
    r <- prepare_manhattan_data(
        mh_df(chr = c("scaffoldA", "2", "1"), pos = c(100L, 200L, 300L),
              pvalue = c(0.1, 0.2, 0.3)),
        pval_col = "pvalue")
    expect_identical(levels(r$data$chr_f), c("1", "2", "scaffoldA"))
    expect_identical(tail(as.character(r$chr_info$chr_f), 1L), "scaffoldA")
})

test_that("chr_len is the MAXIMUM OBSERVED position, not a true chromosome length", {
    # The pipeline never reads a .fai; the axis is built from the SNPs present.
    # Two datasets on the same genome therefore get different chr_len — which is
    # exactly why coords JSON is per-run and not comparable across runs.
    r <- prepare_manhattan_data(mh_df(), pval_col = "pvalue")
    expect_identical(r$chr_info$chr_len, c(500L, 800L))
})

test_that("the inter-chromosome gap is 2 percent of the MEAN chr_len", {
    r <- prepare_manhattan_data(mh_df(), pval_col = "pvalue")
    gap <- mean(c(500, 800)) * 0.02
    # chr 1 starts at 0; chr 2 starts one chr_len plus one gap later.
    expect_equal(r$chr_info$tot, c(0, 500 + gap))
})

test_that("tot, center and pos_cum are consistent with each other", {
    r <- prepare_manhattan_data(mh_df(), pval_col = "pvalue")
    ci <- r$chr_info
    expect_equal(ci$center, ci$tot + ci$chr_len / 2)

    d <- r$data
    tot_by_chr <- setNames(ci$tot, as.character(ci$chr_f))
    expect_equal(d$pos_cum, d$pos + unname(tot_by_chr[as.character(d$chr_f)]))
    # The first chromosome is never offset.
    expect_equal(ci$tot[1], 0)
})

test_that("log10p is -log10(pvalue), and a p-value of exactly 0 yields +Inf", {
    # This is the case add_scatter_layer:111-117 exists to filter: scattermore's
    # C routine aborts on a non-finite coordinate.
    r <- prepare_manhattan_data(
        mh_df(chr = c("1", "1"), pos = c(10L, 20L), pvalue = c(0.01, 0)),
        pval_col = "pvalue")
    expect_equal(r$data$log10p[1], 2)
    expect_true(is.infinite(r$data$log10p[2]))
    expect_gt(r$data$log10p[2], 0)
})

test_that("prepare_manhattan_data honours non-default column names", {
    df <- data.frame(CHROM = c("1", "2"), BP = c(10L, 20L), P = c(0.1, 0.01),
                     stringsAsFactors = FALSE)
    r <- prepare_manhattan_data(df, chr_col = "CHROM", pos_col = "BP", pval_col = "P")
    expect_true(all(c("pos_cum", "log10p", "chr_f") %in% names(r$data)))
    expect_equal(r$data$log10p, c(1, 2))
})

test_that("a single-chromosome dataset gets no offset and no gap", {
    r <- prepare_manhattan_data(
        mh_df(chr = c("1", "1"), pos = c(10L, 90L), pvalue = c(0.5, 0.05)),
        pval_col = "pvalue")
    expect_equal(nrow(r$chr_info), 1L)
    expect_equal(r$chr_info$tot, 0)
    expect_equal(r$data$pos_cum, c(10, 90))
})

test_that("get_chr_colors alternates two blues and respects the requested length", {
    expect_identical(get_chr_colors(1), "#2166AC")
    expect_identical(get_chr_colors(2), c("#2166AC", "#92C5DE"))
    expect_identical(get_chr_colors(5),
                     c("#2166AC", "#92C5DE", "#2166AC", "#92C5DE", "#2166AC"))
    # rep(length.out = 0) is correct here; the other three palette getters are
    # not (see test-known-bugs.R).
    expect_length(get_chr_colors(0), 0L)
})

test_that("get_trait_colors uses named Okabe-Ito up to 8 traits", {
    tc <- get_trait_colors(c("bio_1", "bio_2", "height"))
    expect_identical(names(tc), c("bio_1", "bio_2", "height"))
    expect_identical(unname(tc[1]), "#E69F00")
    expect_length(get_trait_colors(sprintf("t%d", 1:8)), 8L)
})

test_that("get_trait_colors switches to turbo above 8 traits", {
    tc <- get_trait_colors(sprintf("t%d", 1:9))
    expect_length(tc, 9L)
    expect_identical(names(tc), sprintf("t%d", 1:9))
    # Not the Okabe-Ito head any more.
    expect_false(identical(unname(tc[1]), "#E69F00"))
})

test_that("get_region_colors uses the fixed 10 below the cutoff and hue_pal above", {
    expect_identical(get_region_colors(1), "#E41A1C")
    expect_length(get_region_colors(10), 10L)
    expect_length(get_region_colors(11), 11L)
    expect_false(identical(get_region_colors(11)[1], "#E41A1C"))
})

test_that("get_method_shapes returns named ggplot2 shape codes", {
    ms <- get_method_shapes(c("EMMAX", "LFMM", "RDA"))
    expect_identical(names(ms), c("EMMAX", "LFMM", "RDA"))
    expect_identical(unname(ms), c(16, 17, 15))
})

test_that("theme_manhattan hides the legend and blanks the vertical gridlines", {
    th <- theme_manhattan()
    expect_s3_class(th, "theme")
    expect_identical(th$legend.position, "none")
    expect_s3_class(th$panel.grid.major.x, "element_blank")
})

test_that("add_scatter_layer drops non-finite points instead of aborting", {
    r <- prepare_manhattan_data(
        mh_df(chr = c("1", "1", "1"), pos = c(10L, 20L, 30L),
              pvalue = c(0.01, 0, NA)),
        pval_col = "pvalue")
    cols <- get_chr_colors(1)
    names(cols) <- "1"

    lyr <- add_scatter_layer(r$data, cols)
    expect_s3_class(lyr, "Layer")
    # Only the finite row survives: +Inf from pvalue 0 and NA from pvalue NA are
    # both removed before scattermore sees them.
    expect_identical(nrow(lyr$data), 1L)
    expect_true(all(is.finite(lyr$data$log10p)))
})

test_that("add_scatter_layer returns a geom_blank when every point is non-finite", {
    r <- prepare_manhattan_data(
        mh_df(chr = c("1", "1"), pos = c(10L, 20L), pvalue = c(0, 0)),
        pval_col = "pvalue")
    cols <- c("1" = "#2166AC")
    lyr <- add_scatter_layer(r$data, cols)
    expect_s3_class(lyr, "Layer")
    expect_true(inherits(lyr$geom, "GeomBlank"))
})
