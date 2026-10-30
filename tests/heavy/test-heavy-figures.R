# Release 6 figures that need a raster or LEA, so the heavy tier: what the saved plot objects
# hold, beyond the exit-0 rows in heavy_wrappers.R.

test_that("plot_piemap.R draws one round pie per site, at the site's centroid (f3fcce, cdb1a1)", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- heavy_piemap_inputs(d)
    res <- run_wrapper("plot_piemap.R",
                       c(f[1], "bio_1", "bio_1", f[2], f[3], "NULL", "NULL", "0.6", "FALSE", "4",
                         f[4], f[5], "piemap_bio_1", "NULL", "1", "Clusters", "FALSE"), wd = d)
    expect_identical(res$status, 0L, info = res$output)

    g <- qs::qread(file.path(f[5], "piemap_bio_1.qs"))
    is_poly <- vapply(g$layers, function(l) inherits(l$geom, "GeomPolygon"), logical(1))
    expect_identical(sum(is_poly), 1L)
    w <- g$layers[[which(is_poly)]]$data
    expect_setequal(unique(w$pie), c("A", "B", "C"))     # three sites, three pies
    a <- w[w$pie == "A", ]
    # Every slice's polygon starts at its pie's centre: site A's centroid.
    expect_equal(unique(round(a$x[!duplicated(a$wedge)], 9)), 34.7)
    expect_equal(unique(round(a$y[!duplicated(a$wedge)], 9)), 31.1)
    # Round on the map: coord_sf scales latitude by 1/cos(panel mid-latitude); the pie's
    # longitude radius carries the same factor. The arcs are sampled (120 vertices a circle),
    # so the extents are measured to 0.1 %; the old pies were 18 % taller than wide here.
    mid <- mean(range(c(30, 34, w$y)))
    expect_equal(diff(range(a$x)) / diff(range(a$y)), 1 / cos(mid * pi / 180), tolerance = 1e-3)
    expect_match(res$output, "Pie longitude radius x1.1", fixed = TRUE)
})

test_that("plot_pca.R labels each axis with its variance and marks the Tracy-Widom axes (576d53, fff847)", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    # Two populations of 12 with diverged allele frequencies: PC1 carries real structure.
    set.seed(11)
    p1 <- runif(80, 0.1, 0.9); p2 <- pmin(pmax(p1 + rnorm(80, 0, 0.35), 0.02), 0.98)
    g <- rbind(t(replicate(12, rbinom(80, 2, p1))), t(replicate(12, rbinom(80, 2, p2))))
    g <- g[, apply(g, 2, stats::var) > 0]   # LEA's scaled PCA stops (and segfaults) on a constant SNP
    lfmm <- file.path(d, "geno.lfmm")
    write.table(g, lfmm, sep = " ", row.names = FALSE, col.names = FALSE)
    meta <- write_tsv(data.table::data.table(site = rep(c("P1", "P2"), each = 12),
                                             sample = sprintf("I%02d", 1:24),
                                             latitude = 31, longitude = 35),
                      file.path(d, "metadata.tsv"))
    inter <- fx_inter_dir(d, "inter")
    # LEA::pca() writes its project into the working directory, so the run's cwd is d.
    res <- run_wrapper("plot_pca.R", c(lfmm, meta, file.path(d, "pca.png"),
                                       file.path(d, "tracy_widom.png"), inter), wd = d)
    expect_identical(res$status, 0L, info = res$output)
    expect_gt(file.size(file.path(d, "tracy_widom.png")), 0)

    pca <- qs::qread(file.path(inter, "pca.qs"))
    expect_match(pca$labels$x, "^PC1 \\([0-9.]+%\\)$")
    expect_match(pca$labels$y, "^PC2 \\([0-9.]+%\\)$")
    tw <- qs::qread(file.path(inter, "tracy_widom.qs"))
    expect_identical(levels(tw$data$significant), c("p < 0.05", "p ≥ 0.05"))
    expect_identical(as.character(tw$data$significant[1]), "p < 0.05")
    pc1 <- as.numeric(sub("^PC1 \\(([0-9.]+)%\\)$", "\\1", pca$labels$x))
    expect_equal(tw$data$variance[1], pc1, tolerance = 0.05)   # one number, shown twice
    expect_identical(tw$labels$y, "Variance explained (%)")
})
