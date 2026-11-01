# geometric_offset.R — present and future raster cells paired by ID (2c9c67). Heavy tier:
# the script attaches LEA. Same contract as the quick tier's test-id_binding.R: the future
# table permuted must give the unpermuted output byte for byte.

test_that("geometric_offset.R pairs present and future cells by ID", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    w <- fx_offset_world(d)
    run <- function(tag, all_fut) {
        o <- file.path(d, paste0(tag, c("_site.tsv", "_map.tsv", "_offset.tif", "_imp.png",
                                        "_diag.tsv")))
        res <- run_wrapper("geometric_offset.R",
                           c(w$lfmm, w$vcfsnp, w$removed, w$snps, w$site_pres, w$site_fut,
                             w$all_pres, all_fut, w$raster, w$meta, "bio_1,bio_12", "2", "TRUE",
                             o[1], o[2], o[3], o[4], "s1", o[5]))
        list(res = res, site = o[1], map = o[2])
    }
    base <- run("base", w$all_fut)
    expect_identical(base$res$status, 0L, info = base$res$output)
    shuffled <- fx_reorder_rows(w$all_fut, c(51:100, 1:50), file.path(d, "fut_shuffled.tsv"))
    perm <- run("perm", shuffled)
    expect_identical(perm$res$status, 0L, info = perm$res$output)
    expect_identical(readLines(perm$map), readLines(base$map))
    expect_identical(readLines(perm$site), readLines(base$site))
})

# With more than one scenario the LFMM2 fit runs once and scenarios 2.. are projected by the
# script's own §7a, not by LEA::genetic.gap(). Scenario 1 still calls genetic.gap() and the
# job stops if the two disagree. Two identical scenarios: the job passes that guard, and
# scenario 2 (the script's projection) carries scenario 1's (genetic.gap()'s) offsets.
test_that("geometric_offset.R projects later scenarios as genetic.gap() does", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    w <- fx_offset_world(d)
    for (scale in c("TRUE", "FALSE")) {
        o <- function(s, x) file.path(d, paste0(scale, "_", s, x))
        res <- run_wrapper("geometric_offset.R",
                           c(w$lfmm, w$vcfsnp, w$removed, w$snps, w$site_pres,
                             paste(w$site_fut, w$site_fut, sep = ","),
                             w$all_pres, paste(w$all_fut, w$all_fut, sep = ","),
                             w$raster, w$meta, "bio_1,bio_12", "2", scale,
                             paste(o("s1", "_site.tsv"), o("s2", "_site.tsv"), sep = ","),
                             paste(o("s1", "_map.tsv"),  o("s2", "_map.tsv"),  sep = ","),
                             paste(o("s1", "_offset.tif"), o("s2", "_offset.tif"), sep = ","),
                             o("", "imp.png"), "s1,s2", o("", "diag.tsv")))
        expect_identical(res$status, 0L, info = paste(scale, res$output))
        expect_true(any(grepl("hoisted LFMM2 fit vs LEA::genetic.gap", res$output)),
                    info = "the two-scenario call must take the hoisted path")
        s1 <- data.table::fread(o("s1", "_site.tsv")); s2 <- data.table::fread(o("s2", "_site.tsv"))
        expect_identical(s2[, !"genetic_offset"], s1[, !"genetic_offset"], info = scale)
        expect_gt(sum(is.finite(s1$genetic_offset)), 0)
        expect_equal(s2$genetic_offset, s1$genetic_offset, tolerance = 1e-10, info = scale)
        # The map is the raster as a headerless grid (NA cells empty).
        m1 <- as.matrix(data.table::fread(o("s1", "_map.tsv"), header = FALSE))
        m2 <- as.matrix(data.table::fread(o("s2", "_map.tsv"), header = FALSE))
        expect_identical(is.na(m2), is.na(m1), info = scale)
        expect_gt(sum(!is.na(m1)), 0)
        expect_equal(m2, m1, tolerance = 1e-10, info = scale)
    }
})
