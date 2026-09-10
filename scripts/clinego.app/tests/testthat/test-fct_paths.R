# Module constants were renamed: the old single MOD_ASSOC ("association") split
# into MOD_GEA ("GEA") and MOD_GWAS ("GWAS"), and the overlap module became
# MOD_GEAXGWAS ("GEAxGWAS"). See MOD_* in R/fct_paths.R.

test_that("project_base constructs path correctly", {
    p <- project_base("SIMDATA", "/pipeline")
    expect_equal(p, "/pipeline/SIMDATA_results")
})

test_that("mod_path constructs nested paths", {
    p <- mod_path("SIMDATA", MOD_GEA, "plots", "manhattan", "EMMAX")
    expect_equal(p, "/pipeline/SIMDATA_results/GEA/plots/manhattan/EMMAX")
})

test_that("manhattan_bg_path uses correct filename format", {
    p <- manhattan_bg_path("SIMDATA", MOD_GEA, "EMMAX", "bio_1", 3, "bonf_0.05")
    expect_true(grepl("manhattan_bio_1_K3_bonf_0.05_background.png$", p))
})

test_that("combined_manhattan_coords_path is correct", {
    p <- combined_manhattan_coords_path("SIMDATA", MOD_GEA, 3)
    expect_true(grepl("manhattan_combined_K3_coords.json$", p))
})

test_that("miami_bg_path is correct", {
    p <- miami_bg_path("SIMDATA", 3)
    expect_true(grepl("GEAxGWAS/plots/miami_combined_K3_background.png$", p))
})

test_that("pipeline_summary_path is correct", {
    p <- pipeline_summary_path("SIMDATA")
    expect_equal(p, "/pipeline/SIMDATA_results/pipeline_summary.tsv")
})
