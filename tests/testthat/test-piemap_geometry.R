# scripts/R/lib/piemap_geometry.R — where and how plot_piemap.R draws its pies: one pie per site
# at its centroid (cdb1a1), and round on a longitude/latitude map (f3fcce).

test_that("site_centroids gives one row per site at the mean of its samples", {
    s <- data.frame(site      = c("A", "A", "A", "B", "B", "C"),
                    latitude  = c(31.0, 31.2, NA, 32, 32, NA),
                    longitude = c(35.0, 35.4, NA, 34, 34, NA))
    out <- site_centroids(s)
    expect_identical(out$site, c("A", "B"))       # C has no coordinate: not drawn
    expect_equal(out$latitude, c(31.1, 32))
    expect_equal(out$longitude, c(35.2, 34))
})

test_that("pie_x_stretch is 1/cos of the middle latitude on a lon/lat map, 1 otherwise", {
    expect_equal(pie_x_stretch(c(30, 36)), 1 / cos(33 * pi / 180))
    expect_equal(pie_x_stretch(c(36, 30, 33)), 1 / cos(33 * pi / 180))
    expect_identical(pie_x_stretch(c(30, 36), lonlat = FALSE), 1)
    expect_identical(pie_x_stretch(c(30, 36), lonlat = NA), 1)
})

test_that("pie_wedges: slices in proportion, clockwise from 12 o'clock, longitude radius stretched", {
    pies <- data.frame(x = 35, y = 32, r = 0.1, group = "A", C1 = 0.25, C2 = 0.75)
    w <- pie_wedges(pies, c("C1", "C2"), x_stretch = 1.5)
    expect_identical(levels(w$slice), c("C1", "C2"))
    c1 <- w[w$slice == "C1", ]
    # The first quarter runs from north to east: right of and above the centre.
    expect_true(all(c1$x >= 35 - 1e-12) && all(c1$y >= 32 - 1e-12))
    expect_equal(max(c1$x) - 35, 0.1 * 1.5)
    expect_equal(max(c1$y) - 32, 0.1)
    expect_equal(diff(range(w$x)), 2 * 0.1 * 1.5)
    expect_equal(diff(range(w$y)), 2 * 0.1)
})

test_that("pie_wedges leaves out empty and NA slices, and pies with nothing to draw", {
    pies <- data.frame(x = c(35, 36), y = c(32, 33), r = 0.1, group = c("A", "B"),
                       C1 = c(1, NA), C2 = c(0, 0))
    w <- pie_wedges(pies, c("C1", "C2"))
    expect_identical(unique(w$pie), "A")
    expect_identical(unique(as.character(w$slice)), "C1")
    expect_false(any(w$x == 35 & w$y == 32))   # a whole pie is the circle, no centre vertex
    expect_identical(nrow(pie_wedges(pies[2, ], c("C1", "C2"))), 0L)
})

test_that("a pie drawn on a lon/lat raster comes out round (f3fcce)", {
    r <- terra::rast(nrows = 20, ncols = 10, xmin = 34, xmax = 36, ymin = 30, ymax = 34,
                     crs = "EPSG:4326")
    terra::values(r) <- seq_len(terra::ncell(r))
    pies <- data.frame(x = 35, y = 33.5, r = 0.1, group = "A", C1 = 0.5, C2 = 0.5)
    y_range <- range(c(30, 34, pies$y - pies$r, pies$y + pies$r))
    stretch <- pie_x_stretch(y_range, terra::is.lonlat(r))
    expect_gt(stretch, 1.15)                    # 1/cos(32 N): a plain circle would be ovalled
    w <- pie_wedges(pies, c("C1", "C2"), stretch)
    p <- ggplot2::ggplot() +
        ggspatial::layer_spatial(r, ggplot2::aes(fill = ggplot2::after_stat(band1))) +
        ggplot2::geom_polygon(data = w, ggplot2::aes(x = x, y = y, group = wedge))
    b  <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]
    aspect <- b$layout$coord$aspect(pp)         # panel height / width, as coord_sf draws it
    width  <- diff(range(w$x)) / diff(pp$x_range)
    height <- diff(range(w$y)) / diff(pp$y_range) * aspect
    expect_equal(width / height, 1, tolerance = 1e-6)
})
