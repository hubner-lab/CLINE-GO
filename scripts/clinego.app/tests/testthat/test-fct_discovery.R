test_that("find_projects finds _results directories", {
    tmp <- withr::local_tempdir()
    dir.create(file.path(tmp, "SIMDATA_results"))
    dir.create(file.path(tmp, "TEST_results"))
    dir.create(file.path(tmp, "other_dir"))

    projects <- find_projects(tmp)
    expect_setequal(projects, c("SIMDATA", "TEST"))
})

test_that("find_k_values returns sorted integers", {
    tmp <- withr::local_tempdir()
    # Build the tree through MOD_PRESTRUCT rather than a literal. This test used
    # to hardcode "structure/plots/K*", which stopped being where find_k_values()
    # looks when the module directories were reorganised (PreStructure/plots/K*),
    # so it silently asserted against an empty directory. The find_k_range tests
    # below already use the constant, which is why they were unaffected.
    for (k in c("K2", "K3", "K5")) {
        dir.create(file.path(tmp, "SIMDATA_results", MOD_PRESTRUCT, "plots", k),
                   recursive = TRUE)
    }
    withr::local_options(clinego.pipeline_path = tmp)
    ks <- find_k_values("SIMDATA")
    expect_equal(ks, c(2L, 3L, 5L))
})

test_that("find_k_range resolves by config, then mtime — never alphabetically", {
    tmp <- withr::local_tempdir()
    d <- file.path(tmp, "SIMDATA_results", MOD_PRESTRUCT, "plots")
    dir.create(d, recursive = TRUE)
    # K2-6 sorts first but is the stale one; K2-7.0 is the float spelling the
    # Shiny config writer used to produce.
    file.create(file.path(d, "cross_entropy_K2-6.png"))
    Sys.setFileTime(file.path(d, "cross_entropy_K2-6.png"), Sys.time() - 3600)
    file.create(file.path(d, "cross_entropy_K2-7.0.png"))
    withr::local_options(clinego.pipeline_path = tmp)

    cur <- find_k_range("SIMDATA", list(sNMF = list(k_start = 2, k_end = 7)))
    expect_equal(basename(cur$path), "cross_entropy_K2-7.0.png")
    expect_equal(c(cur$k_start, cur$k_end), c(2L, 7L))  # greedy regex used to give NA
    expect_equal(cur$stale, "cross_entropy_K2-6.png")

    back <- find_k_range("SIMDATA", list(sNMF = list(k_start = 2, k_end = 6)))
    expect_equal(basename(back$path), "cross_entropy_K2-6.png")
    expect_equal(back$stale, "cross_entropy_K2-7.0.png")

    none <- find_k_range("SIMDATA", NULL)  # falls back to newest mtime
    expect_equal(basename(none$path), "cross_entropy_K2-7.0.png")
})

test_that("input_to_config_value writes whole numbers as integers", {
    expect_identical(input_to_config_value(7, "numeric"), 7L)
    expect_identical(input_to_config_value("100", "numeric"), 100L)
    expect_identical(input_to_config_value(0.05, "numeric"), 0.05)
    expect_identical(input_to_config_value(2.5, "numeric"), 2.5)
    expect_null(input_to_config_value("", "numeric"))
})

test_that("find_k_range breaks a config tie by mtime, not alphabetically", {
    # Both spellings of one k_end match the same config numerically: a pre-fix run
    # left cross_entropy_K2-7.0.png, the next run writes cross_entropy_K2-7.png
    # beside it. Alphabetically "K2-7.0.png" < "K2-7.png", so taking the first hit
    # would show the stale plot and call the fresh one superseded.
    tmp <- withr::local_tempdir()
    d <- file.path(tmp, "SIMDATA_results", MOD_PRESTRUCT, "plots")
    dir.create(d, recursive = TRUE)
    for (f in c("cross_entropy_K2-6.png", "cross_entropy_K2-7.0.png")) {
        file.create(file.path(d, f))
        Sys.setFileTime(file.path(d, f), Sys.time() - 3600)
    }
    file.create(file.path(d, "cross_entropy_K2-7.png"))
    withr::local_options(clinego.pipeline_path = tmp)

    r <- find_k_range("SIMDATA", list(sNMF = list(k_start = 2, k_end = 7.0)))
    expect_equal(basename(r$path), "cross_entropy_K2-7.png")
    expect_setequal(r$stale, c("cross_entropy_K2-6.png", "cross_entropy_K2-7.0.png"))
    expect_equal(c(r$k_start, r$k_end), c(2L, 7L))
})

# --- maladaptation runs (0dc91e, 0d9efe, a0a580) ------------------------------------

mala_tree <- function(tmp) {
    base <- file.path(tmp, "P_results", "Maladaptation")
    # a gradient_forest run with plots, a zoom folder beside its offset maps
    gf <- file.path(base, "plots", "gradient_forest", "setB_spatial", "year2061-2080_ssp585")
    dir.create(file.path(gf, "zoom", "34p5_35p5_31_32p5"), recursive = TRUE)
    file.create(file.path(gf, "genetic_offset_piemap_notrait.png"))
    file.create(file.path(gf, "zoom", "34p5_35p5_31_32p5", "genetic_offset_piemap_notrait.png"))
    dir.create(file.path(base, "tables", "gradient_forest", "setB_spatial", "year2061-2080_ssp585"),
               recursive = TRUE)
    # emit_plots: false -- tables only, and under the registry KEY, not the engine
    for (m in c("rda_offset", "rda_offset_corrected")) {
        t <- file.path(base, "tables", m, "setA_nospatial", "year2061-2080_ssp585")
        dir.create(t, recursive = TRUE)
        file.create(file.path(t, "genetic_offset_site.tsv"))
    }
    file.create(file.path(base, "tables", "snp_sets_summary.tsv"))
    base
}

test_that("find_mala_models lists table-only runs and every method folder by its key", {
    tmp <- withr::local_tempdir()
    mala_tree(tmp)
    withr::local_options(clinego.pipeline_path = tmp)
    m <- find_mala_models("P")
    expect_identical(unname(m), c("gradient_forest:::setB_spatial",
                                  "rda_offset:::setA_nospatial",
                                  "rda_offset_corrected:::setA_nospatial"))
    expect_identical(names(m)[3], "RDAoff corrected / setA_nospatial")
    # a key parses to itself: rda_offset is not a prefix match of rda_offset_corrected
    expect_identical(parse_model_key(m[[3]])$method, "rda_offset_corrected")
    expect_identical(parse_model_key(m[[2]])$method, "rda_offset")
    expect_identical(mala_engine("rda_offset_corrected"), "rda_offset")
    expect_identical(mala_engine("gradient_forest"), "gradient_forest")
})

test_that("find_gf_zooms finds the zoom folders beside a run's offset maps", {
    tmp <- withr::local_tempdir()
    mala_tree(tmp)
    withr::local_options(clinego.pipeline_path = tmp)
    expect_identical(find_gf_zooms("P", "setB_spatial"), "34p5_35p5_31_32p5")
    expect_identical(find_gf_zooms("P", "setA_nospatial", method = "rda_offset"), character(0))
    p <- gf_offset_zoom_piemap_path("P", "setB_spatial", "34p5_35p5_31_32p5")
    expect_true(file.exists(p))
    expect_identical(basename(p), "genetic_offset_piemap_notrait.png")
})
