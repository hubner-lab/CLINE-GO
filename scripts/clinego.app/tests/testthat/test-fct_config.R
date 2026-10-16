test_that("config_get returns default when key missing", {
    cfg <- list(snmf = list(k_best = 3))
    expect_equal(config_get(cfg, "snmf", "k_best"), 3)
    expect_equal(config_get(cfg, "snmf", "missing_key", default = 10L), 10L)
    expect_null(config_get(cfg, "nonexistent"))
})

test_that("config_k_best extracts k_best", {
    # Group key is "sNMF" (config_k_best reads config_get(config, "sNMF", ...)),
    # matching config_SIMDATA.yaml and inst/config_default.yaml. Lowercase "snmf"
    # is not a config shape the app ever sees.
    cfg <- list(sNMF = list(k_best = 5L))
    expect_equal(config_k_best(cfg), 5L)
})

test_that("config_k_best returns NA on missing", {
    expect_true(is.na(config_k_best(list())))
})

test_that("null coalescing operator works", {
    expect_equal(NULL %||% "default", "default")
    expect_equal("value" %||% "default", "value")
    expect_equal(character(0) %||% "default", "default")
})

test_that("config_mala_predictors falls back to Climate.predictors", {
    clim <- list(Climate = list(predictors = "bio_1,bio_3,bio_5,bio_12"))
    expect_equal(config_mala_predictors(clim), c("bio_1", "bio_3", "bio_5", "bio_12"))
    for (unset in list(NULL, "", "auto", " auto ")) {
        cfg <- clim
        cfg$Maladaptation <- list(predictors = unset)
        expect_equal(config_mala_predictors(cfg), c("bio_1", "bio_3", "bio_5", "bio_12"))
    }
    expect_equal(config_mala_predictors(list()), character(0))
})

test_that("config_mala_predictors reads Maladaptation.predictors as string or list", {
    cfg <- list(Climate = list(predictors = "bio_1,bio_3,bio_5,bio_12"),
                Maladaptation = list(predictors = "bio_1, bio_5"))
    expect_equal(config_mala_predictors(cfg), c("bio_1", "bio_5"))
    cfg$Maladaptation$predictors <- list("bio_3", "bio_12")
    expect_equal(config_mala_predictors(cfg), c("bio_3", "bio_12"))
    cfg$Maladaptation$predictors <- "bio_12"    # R's yaml writer collapses a one-item list
    expect_equal(config_mala_predictors(cfg), "bio_12")
})
