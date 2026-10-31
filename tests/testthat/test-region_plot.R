# scripts/R/lib/region_plot.R — what plot_regionplot.R computes besides its topr call: each
# series' threshold (c34204) and its line, one red with a dash pattern and a legend entry each
# (F7), the window and the region's edges (ad9f15), readable gene labels (e9466a).

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

test_that("threshold_lines: one line per distinct threshold, named after its series in series order (F7)", {
    out <- threshold_lines(c(1e-4, 2e-7, 2e-7, NA),
                           c("bio_2 (RDA)", "bio_10 (EMMAX)", "bio_1 (EMMAX)", "bio_3 (LFMM)"))
    expect_identical(out$thresh, c(2e-7, 1e-4))                 # increasing: the highest line first
    # The caller's series order, not text order: bio_10 stays before bio_1 here (9e5dd0 orders
    # the series; a sort would undo it).
    expect_identical(out$label, c("bio_10 (EMMAX), bio_1 (EMMAX)", "bio_2 (RDA)"))
    expect_identical(out$linetype, c("dashed", "dotdash"))
    # with the methods: a line that carries all of a method's series is named by the method
    by_m <- threshold_lines(c(1e-4, 2e-7, 2e-7, NA),
                            c("bio_2 (RDA)", "bio_10 (EMMAX)", "bio_1 (EMMAX)", "bio_3 (LFMM)"),
                            c("RDA", "EMMAX", "EMMAX", "LFMM"))
    expect_identical(by_m$label, c("EMMAX", "RDA"))
    # d150 region 12: nine EMMAX traits at one Bonferroni line, RDA at its q-value line
    nine <- sprintf("bio_%d (EMMAX)", c(1, 2, 4, 6, 7, 8, 9, 10, 11))
    d150 <- threshold_lines(c(rep(2.821e-7, 9), 1.029e-5), c(nine, "climate_multivariate (RDA)"),
                            c(rep("EMMAX", 9), "RDA"))
    expect_identical(d150$label, c("EMMAX", "RDA"))
    # a method split over two lines names its series on each
    split <- threshold_lines(c(1e-5, 2e-5, 2e-5), c("bio_1 (EMMAX)", "bio_2 (EMMAX)", "bio_2 (LFMM)"),
                             c("EMMAX", "EMMAX", "LFMM"))
    expect_identical(split$label, c("bio_1 (EMMAX)", "bio_2 (EMMAX), LFMM"))
    none <- threshold_lines(c(NA_real_, NA_real_), c("a", "b"))
    expect_identical(nrow(none), 0L)
})

test_that("wrap_series_names breaks between series names, never inside one", {
    nine <- sprintf("bio_%d (EMMAX)", c(1, 2, 4, 6, 7, 8, 9, 10, 11))   # d150 region 12's EMMAX line
    expect_identical(wrap_series_names(nine), paste(
        "bio_1 (EMMAX), bio_2 (EMMAX),", "bio_4 (EMMAX), bio_6 (EMMAX),",
        "bio_7 (EMMAX), bio_8 (EMMAX),", "bio_9 (EMMAX), bio_10 (EMMAX),", "bio_11 (EMMAX)",
        sep = "\n"))
    expect_identical(wrap_series_names("climate_multivariate (RDA)", width = 10),
                     "climate_multivariate (RDA)")
})

test_that("threshold_linetypes: a dash pattern of its own per line, never solid or the edges' dotted", {
    lt <- threshold_linetypes(40)
    expect_identical(lt[1:4], c("dashed", "dotdash", "longdash", "twodash"))
    expect_identical(anyDuplicated(lt), 0L)
    # "13", "44", "73": dotted, dashed and longdash spelled in hex
    expect_false(any(lt %in% c("solid", "dotted", "13", "44", "73")))
    expect_length(threshold_linetypes(0), 0L)
})

test_that("add_threshold_lines: one red layer, each line's pattern keyed to its series (F7)", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = c(1, 5, 8), s = c("a", "b", "c")),
                         ggplot2::aes(x, y, colour = s)) + ggplot2::geom_point()
    lines <- threshold_lines(c(1e-3, 1e-5, 1e-5), c("bio_1 (EMMAX)", "bio_2 (EMMAX)", "clim (RDA)"))
    out <- add_threshold_lines(p, lines, linewidth = 0.5)
    is_h <- vapply(out$layers, function(l) inherits(l$geom, "GeomHline"), logical(1))
    expect_identical(sum(is_h), 1L)
    h <- out$layers[[which(is_h)]]
    expect_identical(rlang::as_label(h$mapping$linetype), "label")
    expect_identical(h$aes_params$colour, CLINEGO_SIGNIFICANCE)
    expect_equal(h$data$yintercept, c(5, 3))
    expect_identical(h$data$label, c("bio_2 (EMMAX), clim (RDA)", "bio_1 (EMMAX)"))
    # NA, not TRUE: TRUE would draw the line across every key of the series legend too
    expect_true(is.na(h$show.legend))
    built <- ggplot2::ggplot_build(out)
    drawn <- built$data[[which(is_h)]]
    expect_identical(drawn$linetype, c("dashed", "dotdash"))
    expect_identical(unique(drawn$colour), CLINEGO_SIGNIFICANCE)
    sc <- built$plot$scales$get_scales("linetype")
    expect_identical(sc$name, "Significance threshold")
    expect_identical(sc$get_labels(), lines$label)                  # the legend's entries, in line order
    expect_identical(unname(sc$map(lines$label)), c("dashed", "dotdash"))

    # The overview: the same lines without a legend, mirrored below zero past topr's ntop.
    ov <- add_threshold_lines(p, lines, linewidth = 0.2, legend = FALSE, mirror = TRUE)
    hv <- ov$layers[[which(vapply(ov$layers, function(l) inherits(l$geom, "GeomHline"),
                                  logical(1)))]]
    expect_equal(hv$data$yintercept, c(5, 3, -5, -3))
    expect_false(hv$show.legend)
    expect_identical(ov$scales$get_scales("linetype")$guide, "none")

    expect_identical(add_threshold_lines(p, threshold_lines(NA_real_, "a"), linewidth = 0.5), p)
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
