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
