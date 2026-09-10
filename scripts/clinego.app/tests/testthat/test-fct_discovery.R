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
