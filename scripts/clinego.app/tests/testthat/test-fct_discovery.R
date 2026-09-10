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
    dir.create(file.path(tmp, "SIMDATA_results", "structure", "plots", "K2"),
               recursive = TRUE)
    dir.create(file.path(tmp, "SIMDATA_results", "structure", "plots", "K3"),
               recursive = TRUE)
    dir.create(file.path(tmp, "SIMDATA_results", "structure", "plots", "K5"),
               recursive = TRUE)
    # Override pipeline path
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
