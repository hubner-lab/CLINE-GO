# scripts/R/utils/theme_clinego.R — the shared plot look and the semantic palette.
#
# Worth testing despite being "just styling": the semantic constants are a
# CONTRACT (CLAUDE.md Code Structure Rule 9 — "always use these constants, never
# raw hex, so meaning stays consistent everywhere"), and the file's own comments
# claim two convergences that nothing checked: CLINEGO_CATEGORICAL is asserted to
# be identical to manhattan_utils.R's Okabe-Ito vector, and theme_clinego_grid()
# is asserted to blank the axis line that theme_clinego() draws.

test_that("the semantic color constants hold their documented meanings", {
    expect_identical(CLINEGO_RETAINED,  "#0B775E")   # green  — passed a filter
    expect_identical(CLINEGO_REMOVED,   "#F2300F")   # red    — dropped
    expect_identical(CLINEGO_THRESHOLD, "#35274A")   # plum   — cutoff lines
    # Aliased on purpose: "no more grey dots" (user preference, documented :34-35).
    expect_identical(CLINEGO_NEUTRAL, CLINEGO_RETAINED)
})

test_that("CLINEGO_CATEGORICAL is the same Okabe-Ito palette manhattan_utils.R uses", {
    # The file header claims these two palettes were converged. get_trait_colors()
    # holds the other copy (manhattan_utils.R:82-83); if either is edited alone,
    # a trait changes color between a Manhattan and a processing plot.
    from_manhattan <- unname(get_trait_colors(sprintf("t%d", 1:8)))
    expect_identical(CLINEGO_CATEGORICAL, from_manhattan)
    expect_length(CLINEGO_CATEGORICAL, 8L)
})

test_that("clinego_cluster_palette returns the first n colors within the palette", {
    expect_identical(clinego_cluster_palette(1), "#0072B2")
    expect_identical(clinego_cluster_palette(3), CLINEGO_CLUSTERS[1:3])
    expect_identical(clinego_cluster_palette(length(CLINEGO_CLUSTERS)),
                     CLINEGO_CLUSTERS)
})

test_that("clinego_cluster_palette interpolates beyond the palette length", {
    n <- length(CLINEGO_CLUSTERS) + 5L
    got <- clinego_cluster_palette(n)
    expect_length(got, n)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}", got)))
})

test_that("clinego_cluster_palette contains no grey", {
    # CLINEGO_NEUTRAL owns "no flag"; a grey cluster would read as one.
    expect_false("#999999" %in% CLINEGO_CLUSTERS)
})

test_that("theme_clinego draws axis lines and puts the legend at the bottom", {
    th <- theme_clinego()
    expect_s3_class(th, "theme")
    expect_s3_class(th$axis.line, "element_line")
    expect_identical(th$legend.position, "bottom")
})

test_that("theme_clinego scales its type sizes off base_size", {
    th <- theme_clinego(base_size = 20)
    expect_identical(th$plot.title$size, 22)      # base + 2
    expect_identical(th$plot.subtitle$size, 18)   # base - 2
    expect_identical(th$axis.title$size, 20)
})

test_that("theme_clinego_grid swaps the axis line for a panel border", {
    # This is the whole difference between the two themes, and it is easy to
    # break by reordering the + composition.
    base <- theme_clinego()
    grid <- theme_clinego_grid()
    expect_s3_class(base$axis.line, "element_line")
    expect_s3_class(grid$axis.line, "element_blank")
    expect_s3_class(grid$panel.border, "element_rect")
    expect_s3_class(grid$strip.background, "element_rect")
})

test_that("scale_color_clinego and scale_fill_clinego carry the shared palette", {
    sc <- scale_color_clinego()
    sf <- scale_fill_clinego()
    expect_s3_class(sc, "ggproto")
    expect_s3_class(sf, "ggproto")
    expect_identical(sc$palette(8), CLINEGO_CATEGORICAL)
    expect_identical(sf$palette(8), CLINEGO_CATEGORICAL)
})

test_that("scale_*_clinego forwards its arguments to the underlying scale", {
    sc <- scale_color_clinego(name = "Trait")
    expect_identical(sc$name, "Trait")
})

test_that("clinego_empty_plot builds a real plot carrying the message", {
    p <- clinego_empty_plot("Not enough samples")
    expect_s3_class(p, "ggplot")
    # The message must survive into the layer data, not just the call.
    labels <- unlist(lapply(p$layers, function(l) l$aes_params$label))
    expect_true("Not enough samples" %in% labels)
})

test_that("clinego_save_both writes a non-empty PNG and SVG pair", {
    # Pins the documented defense at :120-125: the old file.create() stub pattern
    # produced 0-byte PNGs that render as broken images.
    d <- withr::local_tempdir()
    stem <- file.path(d, "figure")
    p <- clinego_empty_plot("x")
    quiet(clinego_save_both(stem, p, w = 4, h = 3))

    for (ext in c(".png", ".svg")) {
        f <- paste0(stem, ext)
        expect_true(file.exists(f))
        expect_gt(file.size(f), 0)
    }
})
