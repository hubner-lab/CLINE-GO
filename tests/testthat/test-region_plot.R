# scripts/R/lib/region_plot.R — what plot_regionplot.R computes besides its topr call: each
# series' threshold (c34204), the window and the region's edges (ad9f15), readable gene labels
# (e9466a).

test_that("series_threshold is compute_pval_threshold's, NA rows not counted (c34204)", {
    p <- c(1e-9, 0.2, NA, 0.5, 0.01)
    expect_equal(series_threshold(p, "bonf_0.05"), 0.05 / 4)
    expect_identical(series_threshold(p, "bonf_0.05"),
                     compute_pval_threshold(p, "bonf", 0.05)$threshold)
    expect_identical(series_threshold(p, "top_2"), compute_pval_threshold(p, "top", 2)$threshold)
    expect_identical(series_threshold(p, "custom_1e-5"), 1e-5)
    expect_true(is.na(series_threshold(rep(NA_real_, 3), "bonf_0.05")))
})

test_that("each trait gets its own threshold, not the first trait column's", {
    expect_identical(series_threshold(c(1e-6, 0.5, 0.7), "top_1"), 1e-6)
    expect_identical(series_threshold(c(1e-3, 0.5, 0.7), "top_1"), 1e-3)
})

test_that("threshold_lines: every distinct threshold, in its only series' colour or shared", {
    out <- threshold_lines(c(1e-4, 2e-7, 2e-7, NA), c("red", "blue", "green", "grey"), "plum")
    expect_identical(out$thresh, c(2e-7, 1e-4))
    expect_identical(out$color, c("plum", "red"))
    one <- threshold_lines(c(2e-7, 2e-7), c("blue", "blue"), "plum")
    expect_identical(one$color, "blue")
    none <- threshold_lines(c(NA_real_, NA_real_), c("a", "b"), "plum")
    expect_length(none$thresh, 0)
})

test_that("region_window adds the flank and region_edge_lines marks the region inside it (ad9f15)", {
    expect_identical(region_window(1000, 2000, 500), c(500, 2500))
    expect_identical(region_window(300, 2000, 500), c(1, 2500))
    expect_identical(region_window(1000, 2000), c(1000, 2000))
    expect_null(region_edge_lines(1000, 2000, 1000, 2000))
    expect_identical(region_edge_lines(1000, 2000, 500, 2500), c(1000, 2000))
    expect_identical(region_edge_lines(1, 2000, 1, 2500), 2000)
})

test_that("repel_gene_labels puts each name right of its box and keeps every one (e9466a)", {
    genes <- data.frame(gene_start = c(100, 150, 300), gene_end = c(140, 190, 400), y = 1,
                        gene_symbol = c("RCC1", "RCC1", ""))
    gp <- ggplot2::ggplot(genes) +
        ggplot2::geom_segment(ggplot2::aes(x = gene_start, xend = gene_end, y = y, yend = y)) +
        ggplot2::geom_text(data = genes, ggplot2::aes(x = gene_start, y = y, label = gene_symbol),
                           size = 2.5)
    out <- repel_gene_labels(gp, "black")
    l <- out$layers[[2]]
    expect_s3_class(l$geom, "GeomTextRepel")
    expect_identical(l$data$gene_symbol, c("RCC1", "RCC1"))       # the unnamed gene pushes nothing
    expect_equal(l$data$gene_start, c(140, 190) + 0.004 * 300)   # just right of each box
    expect_identical(l$aes_params$size, 2.5)
    expect_s3_class(out$layers[[1]]$geom, "GeomSegment")
})

test_that("repel_gene_labels refuses a text layer without topr's gene columns", {
    gp <- ggplot2::ggplot() +
        ggplot2::geom_text(data = data.frame(x = 1, y = 1, lab = "a"),
                           ggplot2::aes(x = x, y = y, label = lab))
    expect_error(repel_gene_labels(gp, "black"), "gene_symbol")
})
