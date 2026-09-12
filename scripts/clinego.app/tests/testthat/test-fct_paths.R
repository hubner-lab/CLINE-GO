# Module constants were renamed: the old single MOD_ASSOC ("association") split
# into MOD_GEA ("GEA") and MOD_GWAS ("GWAS"), and the overlap module became
# MOD_GEAXGWAS ("GEAxGWAS"). See MOD_* in R/fct_paths.R.

test_that("project_base constructs path correctly", {
    # Explicit root argument, so this one is independent of the option.
    p <- project_base("SIMDATA", "/pipeline")
    expect_equal(p, "/pipeline/SIMDATA_results")
})

test_that("mod_path constructs nested paths", {
    # Pin the root explicitly: helper-pipeline-path.R redirects
    # getOption("clinego.pipeline_path") to a tempdir for the whole suite, and this
    # block asserts the {project}_results/<module>/... convention as a WHOLE literal
    # string, which is worth keeping. Path builders are pure, so nothing is created.
    withr::local_options(clinego.pipeline_path = "/pipeline")
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
    withr::local_options(clinego.pipeline_path = "/pipeline")
    p <- pipeline_summary_path("SIMDATA")
    expect_equal(p, "/pipeline/SIMDATA_results/pipeline_summary.tsv")
})

test_that("the path builders follow getOption(\"clinego.pipeline_path\")", {
    # The mechanism the app suite's hermeticity now rests on: every
    # {PROJECT}_results/ path resolves through get_pipeline_path(), which honours
    # this option ahead of the golem config, the env var and the /pipeline
    # fallback (R/app_config.R:33-42). If that priority is ever reordered, the
    # suite silently starts writing into the mounted repo again.
    withr::local_options(clinego.pipeline_path = "/tmp/somewhere-else")
    expect_equal(project_base("P"), "/tmp/somewhere-else/P_results")
    expect_equal(pipeline_summary_path("P"),
                 "/tmp/somewhere-else/P_results/pipeline_summary.tsv")
})
