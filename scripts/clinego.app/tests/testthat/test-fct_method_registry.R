# R/fct_method_registry.R — the pure half: sentinel resolution and override merge.
#
# The "@k_best" sentinel is the mechanism that keeps existing configs behaving
# identically across an upgrade (workflow/methods/gea.py's header says so): LFMM K,
# EMMAX #PCs and GAPIT PCA.total all equal sNMF.k_best today, and the sentinel
# preserves that without hardcoding a literal that would diverge from a project's
# own k_best. If it stops resolving, the GEA tab silently offers 3 to a project
# whose k_best is 5 — a different model, no error.
#
# These defaults are ALSO injected into the row-editor JS as PARAM_DEFAULTS
# (:117-122), so an R-created row and a JS-created row must agree. That makes this
# a contract between two languages, tested on the R side here.
#
# gea_method_registry()/gwas_method_registry() are NOT tested: they shell out to
# python3 and cache the result — including the degraded fallback — under a
# pipeline_path + mtime key. Deferred to the testing track.

fake_registry <- function() list(
    LFMM = list(
        engine = "lfmm",
        adjust_default = "qval", threshold_default = 0.05,
        significance_family = "univariate_pvalue",
        params = list(
            K       = list(type = "int", default = "@k_best", min = 1),
            n_boot  = list(type = "int", default = 100))),
    EMMAX = list(
        engine = "emmax",
        adjust_default = "bonf", threshold_default = "0.05",
        significance_family = "univariate_pvalue",
        params = list(n_pcs = list(type = "int", default = "@k_best"))),
    RDA = list(
        engine = "rda",
        significance_family = "multivariate_pvalue",
        params = list())
)

test_that("gea_param_default resolves the @k_best sentinel", {
    spec <- list(type = "int", default = "@k_best")
    expect_identical(gea_param_default(spec, k_best = 5), 5)
    expect_identical(gea_param_default(spec, k_best = 2), 2)
})

test_that("gea_param_default falls back to 3 when k_best is unknown", {
    # 3 is the documented fallback. A NULL leaking through instead would put an
    # empty numericInput in front of the user.
    spec <- list(type = "int", default = "@k_best")
    expect_identical(gea_param_default(spec, k_best = NULL), 3)
    expect_identical(gea_param_default(spec), 3)
})

test_that("gea_param_default passes a literal default through untouched", {
    expect_identical(gea_param_default(list(type = "int", default = 100), k_best = 5), 100)
    expect_identical(gea_param_default(list(type = "bool", default = FALSE), k_best = 5),
                     FALSE)
    # Only the exact sentinel string is special.
    expect_identical(gea_param_default(list(default = "k_best"), k_best = 5), "k_best")
    expect_identical(gea_param_default(list(default = "@K_BEST"), k_best = 5), "@K_BEST")
})

test_that("gea_method_param_defaults resolves every method's params", {
    got <- gea_method_param_defaults(fake_registry(), k_best = 4)
    expect_named(got, c("LFMM", "EMMAX", "RDA"))
    expect_identical(got$LFMM$K, 4)          # sentinel resolved
    expect_identical(got$LFMM$n_boot, 100)   # literal preserved
    expect_identical(got$EMMAX$n_pcs, 4)
})

test_that("gea_method_param_defaults gives a paramless method an empty list", {
    got <- gea_method_param_defaults(fake_registry(), k_best = 4)
    expect_type(got$RDA, "list")
    expect_length(got$RDA, 0L)
})

test_that("resolve_row_params starts from the registry defaults", {
    got <- resolve_row_params("LFMM", NULL, fake_registry(), k_best = 6)
    expect_identical(got$K, 6)
    expect_identical(got$n_boot, 100)
})

test_that("resolve_row_params lets a config row override a default", {
    got <- resolve_row_params("LFMM", list(K = 9), fake_registry(), k_best = 6)
    expect_identical(got$K, 9)
    expect_identical(got$n_boot, 100)   # untouched
})

test_that("resolve_row_params drops unknown param names silently", {
    # common.smk's resolve_method_params() is the authoritative validator; this side
    # is for display, and must not invent a row the editor cannot render.
    got <- resolve_row_params("LFMM", list(K = 9, not_a_param = 1),
                              fake_registry(), k_best = 6)
    expect_named(got, c("K", "n_boot"))
    expect_false("not_a_param" %in% names(got))
})

test_that("resolve_row_params returns an empty list for an unknown method", {
    expect_length(resolve_row_params("NOPE", list(K = 1), fake_registry(), k_best = 6), 0L)
})

test_that("resolve_row_params and gea_method_param_defaults agree with no overrides", {
    # The R row renderer uses one, the injected JS uses the other. They must not
    # drift, which is exactly the failure that would produce a row whose displayed
    # default differs from the value a JS-created row gets.
    reg <- fake_registry()
    per_method <- gea_method_param_defaults(reg, k_best = 7)
    for (m in names(reg)) {
        expect_identical(resolve_row_params(m, NULL, reg, k_best = 7), per_method[[m]],
                         info = paste0("defaults disagree for ", m))
    }
})

test_that("gea_method_significance_defaults reads the registry's own fields", {
    got <- gea_method_significance_defaults(fake_registry())
    expect_identical(got$LFMM$adjust, "qval")
    expect_identical(got$LFMM$family, "univariate_pvalue")
    expect_identical(got$RDA$family, "multivariate_pvalue")
})

test_that("gea_method_significance_defaults always returns threshold as CHARACTER", {
    # The config YAML stores thresholds as strings and the value becomes a filename
    # component (the `adjust` wildcard in _assoc_downstream.smk), so a numeric here
    # would render as "0.05" in one place and "0.050000" in another.
    got <- gea_method_significance_defaults(fake_registry())
    expect_type(got$LFMM$threshold, "character")   # was numeric 0.05 in the registry
    expect_identical(got$LFMM$threshold, "0.05")
    expect_type(got$EMMAX$threshold, "character")
})

test_that("gea_method_significance_defaults falls back for a degraded registry", {
    # The python shell-out can fail, leaving entries without the significance keys.
    degraded <- list(EMMAX = list(engine = "emmax"))
    got <- gea_method_significance_defaults(degraded)
    expect_identical(got$EMMAX, list(adjust = "bonf", threshold = "0.05",
                                     family = "univariate_pvalue"))
})
