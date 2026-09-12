# R/fct_config_writer.R — prepare_config_for_yaml(), which every config the app
# writes passes through on its way to disk.
#
# Its job is to make an R list safe for yaml::write_yaml(): drop NULLs (a NULL key
# must be absent, not `~`), turn NA into NULL, and — the part that is easy to break
# — treat NAMED and UNNAMED lists differently. A named list has NULL entries
# removed (:96); an unnamed one keeps every entry (:101), because association.configs
# is a YAML SEQUENCE whose positions carry meaning. Collapsing one entry out of that
# sequence silently drops a method from the run.

test_that("prepare_config_for_yaml drops NULL entries from a named list", {
    got <- prepare_config_for_yaml(list(a = 1, b = NULL, c = "x"))
    expect_named(got, c("a", "c"))
    expect_identical(got$a, 1)
    expect_identical(got$c, "x")
})

test_that("prepare_config_for_yaml converts a scalar NA to NULL", {
    got <- prepare_config_for_yaml(list(a = NA, b = 2))
    # An NA would round-trip as YAML `.na` and then read back as a string in some
    # parsers; absence is the intended spelling for "not set".
    expect_named(got, "b")
})

test_that("prepare_config_for_yaml recurses into nested named lists", {
    got <- prepare_config_for_yaml(
        list(sNMF = list(k_start = 2L, k_best = NULL),
             Filter = list(maf = 0.05, snp_miss = NA)))
    expect_named(got$sNMF, "k_start")
    expect_named(got$Filter, "maf")
})

test_that("prepare_config_for_yaml removes a named list that empties out", {
    # Every child dropped means the parent is an empty named list, which write_yaml
    # renders as `{}`. It survives as an empty list rather than vanishing — pinned
    # so the behaviour is known rather than assumed.
    got <- prepare_config_for_yaml(list(Climate = list(a = NULL, b = NA), keep = 1))
    expect_true("keep" %in% names(got))
    expect_true(is.list(got$Climate))
    expect_length(got$Climate, 0L)
})

test_that("prepare_config_for_yaml keeps every entry of an UNNAMED list", {
    # association.configs is the case this branch exists for: a YAML sequence whose
    # entries are positional. Dropping one would remove a method from the run.
    configs <- list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"),
        list(method = "LFMM",  adjust = "qval", threshold = "0.05"))
    got <- prepare_config_for_yaml(list(GEA = list(configs = configs)))
    expect_length(got$GEA$configs, 2L)
    expect_identical(got$GEA$configs[[1]]$method, "EMMAX")
    expect_identical(got$GEA$configs[[2]]$method, "LFMM")
})

test_that("an unnamed list keeps its length even when an entry becomes empty", {
    # The asymmetry between :96 and :101 in one assertion: the inner named list is
    # pruned to nothing, but the SEQUENCE still has two positions.
    configs <- list(list(method = "EMMAX"), list(method = NULL, adjust = NA))
    got <- prepare_config_for_yaml(configs)
    expect_length(got, 2L)
    expect_identical(got[[1]]$method, "EMMAX")
    expect_length(got[[2]], 0L)
})

test_that("prepare_config_for_yaml leaves scalars and vectors alone", {
    expect_identical(prepare_config_for_yaml(5L), 5L)
    expect_identical(prepare_config_for_yaml("SIMDATA"), "SIMDATA")
    expect_identical(prepare_config_for_yaml(TRUE), TRUE)
    # A multi-element vector must NOT hit the is.na() branch, which is guarded on
    # length == 1 — c(1, NA) is a real value, not "unset".
    expect_identical(prepare_config_for_yaml(c(1, NA)), c(1, NA))
})

test_that("prepare_config_for_yaml preserves a character vector of predictors", {
    got <- prepare_config_for_yaml(list(Climate = list(predictors = c("bio_1", "bio_12"))))
    expect_identical(got$Climate$predictors, c("bio_1", "bio_12"))
})

test_that("config_file_path names config_{project}.yaml under the pipeline root", {
    # Resolves through get_pipeline_path(), which helper-pipeline-path.R redirects,
    # so pin the root here to assert the naming convention itself.
    withr::local_options(clinego.pipeline_path = "/pipeline")
    expect_identical(config_file_path("SIMDATA"), "/pipeline/config_SIMDATA.yaml")
})

test_that("a prepared config round-trips through yaml", {
    # End-to-end: the point of the function is that write_yaml accepts its output
    # and read_yaml gives the same structure back.
    d <- withr::local_tempdir()
    f <- file.path(d, "config_TEST.yaml")
    cfg <- list(
        project_name = "TEST",
        cpu = 4L,
        Input = list(dir = "data/", vcf = "x.vcf", gff = NULL),
        Filter = list(maf = 0.05, snp_miss = NA),
        GEA = list(configs = list(
            list(method = "EMMAX", adjust = "bonf", threshold = "0.05"),
            list(method = "LFMM",  adjust = "qval", threshold = "0.05"))))

    yaml::write_yaml(prepare_config_for_yaml(cfg), f)
    back <- yaml::read_yaml(f)

    expect_identical(back$project_name, "TEST")
    expect_identical(back$cpu, 4L)
    expect_false("gff" %in% names(back$Input))
    expect_false("snp_miss" %in% names(back$Filter))
    expect_length(back$GEA$configs, 2L)
    expect_identical(back$GEA$configs[[2]]$method, "LFMM")
})
