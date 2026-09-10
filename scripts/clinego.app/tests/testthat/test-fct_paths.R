test_that("project_base constructs path correctly", {
    p <- project_base("SIMDATA", "/pipeline")
    expect_equal(p, "/pipeline/SIMDATA_results")
})

test_that("mod_path constructs nested paths", {
    p <- mod_path("SIMDATA", MOD_ASSOC, "plots", "manhattan", "EMMAX")
    expect_equal(p, "/pipeline/SIMDATA_results/association/plots/manhattan/EMMAX")
})

test_that("manhattan_bg_path uses correct filename format", {
    p <- manhattan_bg_path("SIMDATA", MOD_ASSOC, "EMMAX", "bio_1", 3, "bonf_0.05")
    expect_true(grepl("manhattan_bio_1_K3_bonf_0.05_background.png$", p))
})

test_that("combined_manhattan_coords_path is correct", {
    p <- combined_manhattan_coords_path("SIMDATA", MOD_ASSOC, 3)
    expect_true(grepl("manhattan_combined_K3_coords.json$", p))
})

test_that("miami_bg_path is correct", {
    p <- miami_bg_path("SIMDATA", 3)
    expect_true(grepl("overlapping/plots/miami_combined_K3_background.png$", p))
})

test_that("pipeline_summary_path is correct", {
    p <- pipeline_summary_path("SIMDATA")
    expect_equal(p, "/pipeline/SIMDATA_results/pipeline_summary.tsv")
})
