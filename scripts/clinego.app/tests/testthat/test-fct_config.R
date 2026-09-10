test_that("config_get returns default when key missing", {
    cfg <- list(snmf = list(k_best = 3))
    expect_equal(config_get(cfg, "snmf", "k_best"), 3)
    expect_equal(config_get(cfg, "snmf", "missing_key", default = 10L), 10L)
    expect_null(config_get(cfg, "nonexistent"))
})

test_that("config_k_best extracts k_best", {
    cfg <- list(snmf = list(k_best = 5L))
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
