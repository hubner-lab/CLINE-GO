# Per-method significance-rule helpers (R/fct_threshold_rules.R) + the
# additive `overrides` param on compute_method_thresholds() /
# compute_method_sigsnps_cached() (R/fct_data_loading.R).

test_that("threshold_overrides_key is order-independent (sorted)", {
    a <- list(RDA = list(type = "bonf", value = 0.01), LFMM = list(type = "qval", value = 0.1))
    b <- list(LFMM = list(type = "qval", value = 0.1), RDA = list(type = "bonf", value = 0.01))
    expect_equal(threshold_overrides_key(a), threshold_overrides_key(b))
})

test_that("threshold_overrides_key is empty string for no overrides", {
    expect_equal(threshold_overrides_key(list()), "")
})

test_that("threshold_overrides_key differs when a value differs", {
    a <- list(RDA = list(type = "bonf", value = 0.01))
    b <- list(RDA = list(type = "bonf", value = 0.05))
    expect_false(identical(threshold_overrides_key(a), threshold_overrides_key(b)))
})

test_that("effective_rule_for follows master when method is absent from overrides", {
    r <- effective_rule_for("EMMAX", overrides = list(RDA = list(type = "bonf", value = 0.01)),
                            master_type = "qval", master_value = 0.1)
    expect_equal(r$type, "qval")
    expect_equal(r$value, 0.1)
    expect_equal(r$source, "master")
})

test_that("effective_rule_for uses the override when present", {
    r <- effective_rule_for("RDA", overrides = list(RDA = list(type = "bonf", value = 0.01)),
                            master_type = "qval", master_value = 0.1)
    expect_equal(r$type, "bonf")
    expect_equal(r$value, 0.01)
    expect_equal(r$source, "override")
})

test_that("normalize_threshold_overrides handles NULL/empty/legacy shapes", {
    expect_equal(normalize_threshold_overrides(NULL), list())
    expect_equal(normalize_threshold_overrides(list()), list())
    raw <- list(RDA = list(type = "bonf", value = 0.01))
    out <- normalize_threshold_overrides(raw)
    expect_equal(out$RDA$type, "bonf")
    expect_equal(out$RDA$value, 0.01)
})

test_that("normalize_threshold_overrides drops malformed entries without error", {
    raw <- list(RDA = list(type = "nope", value = 0.01))   # invalid type/value pair for "nope"
    expect_equal(normalize_threshold_overrides(raw), list())

    raw2 <- list(RDA = list(type = "bonf"))   # missing value
    expect_equal(normalize_threshold_overrides(raw2), list())
})

# ── compute_method_thresholds() — additive overrides param ─────────────────

test_that("compute_method_thresholds with no overrides arg is identical to pre-rework output", {
    pv <- list(EMMAX = data.table::data.table(
        SNPID = paste0("s", 1:20), chr = "1", pos = 1:20,
        bio_1 = c(rep(0.001, 5), runif(15, 0.1, 0.9))
    ))
    out_default  <- compute_method_thresholds(pv, "bonf", 0.05)
    out_explicit <- compute_method_thresholds(pv, "bonf", 0.05, overrides = list())
    expect_equal(out_default, out_explicit)
    expect_equal(out_default[["bio_1::EMMAX"]], 0.05 / 20)
})

test_that("compute_method_thresholds override changes only the targeted method", {
    pv <- list(
        EMMAX = data.table::data.table(SNPID = paste0("s", 1:20), chr = "1", pos = 1:20,
                                       bio_1 = runif(20, 0, 1)),
        RDA   = data.table::data.table(SNPID = paste0("s", 1:20), chr = "1", pos = 1:20,
                                       climate_multivariate = runif(20, 0, 1))
    )
    out <- compute_method_thresholds(pv, "bonf", 0.05,
                                     overrides = list(RDA = list(type = "bonf", value = 0.01)))
    expect_equal(out[["bio_1::EMMAX"]], 0.05 / 20)
    expect_equal(out[["climate_multivariate::RDA"]], 0.01 / 20)
})

# ── compute_method_sigsnps_cached() — R4 regression: NA-cutoff fallback must
#    use the METHOD'S OWN rule, not the master, once overrides exist.
#    (Unreachable before per-method rules existed — this is the bug that
#    becomes live the moment overrides are non-empty.)

test_that("compute_method_sigsnps_cached NA-cutoff fallback uses the per-method rule", {
    skip_if_not(exists("compute_pval_threshold", mode = "function"),
               "compute_pval_threshold not sourced (pval_threshold.R not on path)")

    set.seed(1)
    pv <- list(
        # <10 tests -> qval is "too_few_tests" (NA) at the MASTER rule; but the
        # override for this method requests bonf, which is always computable.
        RDA = data.table::data.table(SNPID = paste0("s", 1:5), chr = "1", pos = 1:5,
                                     climate_multivariate = c(0.001, 0.2, 0.4, 0.6, 0.8))
    )
    project <- paste0("test_threshold_rules_", as.integer(stats::runif(1, 1, 1e6)))

    # cutoffs = NULL forces the cold path (no precomputed combo_thresholds),
    # which is where the R4 bug lived.
    out <- compute_method_sigsnps_cached(
        pvalues_list = pv, type = "qval", value = 0.1,
        k = 3, regime = "snp", project = project, module = "GEA",
        cutoffs = NULL,
        overrides = list(RDA = list(type = "bonf", value = 0.05))
    )
    # Master (qval on 5 tests) would yield NA -> empty. The override (bonf
    # 0.05/5 = 0.01) should instead select the one SNP with p=0.001.
    expect_true("RDA" %in% names(out))
    expect_equal(nrow(out$RDA), 1L)
    expect_equal(out$RDA$SNPID, "s1")

    # Clean up the disk cache this test created.
    cache_dir <- interactive_sigsnps_dir(project, "GEA")
    if (dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE)
})

test_that("compute_method_sigsnps_cached with no overrides arg matches pre-rework single-global behaviour", {
    skip_if_not(exists("compute_pval_threshold", mode = "function"),
               "compute_pval_threshold not sourced (pval_threshold.R not on path)")

    pv <- list(EMMAX = data.table::data.table(
        SNPID = paste0("s", 1:20), chr = "1", pos = 1:20,
        bio_1 = c(0.0001, runif(19, 0.1, 0.9))
    ))
    project <- paste0("test_threshold_rules_", as.integer(stats::runif(1, 1, 1e6)))

    out_default  <- compute_method_sigsnps_cached(pv, "bonf", 0.05, k = 3, regime = "snp",
                                                  project = project, module = "GEA", cutoffs = NULL)
    cache_dir <- interactive_sigsnps_dir(project, "GEA")
    if (dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE)

    out_explicit <- compute_method_sigsnps_cached(pv, "bonf", 0.05, k = 3, regime = "snp",
                                                  project = project, module = "GEA", cutoffs = NULL,
                                                  overrides = list())
    if (dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE)

    expect_equal(out_default, out_explicit)
})
