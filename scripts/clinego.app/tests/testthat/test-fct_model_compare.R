# R/fct_model_compare.R — the app side of scripts/compare_offsets.R.

test_that("compare_scenario_label builds the same label the launcher and the novelty tab share", {
    # Until 2026-09-13 the launcher wrote the ExDet cache under "ssp585_2061-2080" while
    # the novelty tab read "future" (audit B7); one helper now serves both.
    cfg <- list(Future = list(ssp = list("585"), year = list("2061-2080")))
    s <- compare_scenario_label(cfg)
    expect_equal(s$ssp, "585")
    expect_equal(s$year, "2061-2080")
    expect_equal(s$label, "ssp585_2061-2080")
})

test_that("compare_scenario_label falls back to ssp585 / 2080 when Future is absent", {
    s <- compare_scenario_label(list())
    expect_equal(s$label, "ssp585_2080")
})

test_that("compare_scenario_label takes the FIRST scenario when several are configured", {
    cfg <- list(Future = list(ssp = list("126", "585"), year = list("2041-2060", "2061-2080")))
    expect_equal(compare_scenario_label(cfg)$label, "ssp126_2041-2060")
})

test_that("load_gf_diagnostics returns the long table as a named list, and list() when absent", {
    d <- withr::local_tempdir()
    withr::local_options(clinego.pipeline_path = d)
    expect_identical(load_gf_diagnostics("P", "s"), list())
    p <- gf_diagnostics_path("P", "s")
    dir.create(dirname(p), recursive = TRUE)
    writeLines(c("quantity\tvalue", "extrap\tTRUE", "pct_future_cells_outside_range__year2061-2080_ssp585\t65.7"), p)
    got <- load_gf_diagnostics("P", "s")
    expect_equal(got[["extrap"]], "TRUE")
    expect_equal(got[["pct_future_cells_outside_range__year2061-2080_ssp585"]], "65.7")
})
