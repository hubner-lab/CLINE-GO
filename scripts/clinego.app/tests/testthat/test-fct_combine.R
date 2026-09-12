# R/fct_combine.R — the significance-rule helpers.
#
# These four decide what the GEA/GWAS tabs offer, pre-fill and display, and the
# threshold VALUE means something different in each mode: an alpha, an FDR target,
# a count of SNPs per trait, a raw p cutoff. A validator that accepts 100 as a
# Bonferroni alpha, or a badge that prints "Top 0.05 SNPs", is wrong in a way no
# error surfaces — the user just gets a selection they did not ask for.
#
# combine_sigsnps() itself is NOT tested here: it is compared against the
# pipeline's implementation in tests/testthat/test-equivalence-app-pipeline.R,
# which is the stronger assertion.

test_that("default_threshold reads the FIRST configs entry", {
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "qval", threshold = "0.1"),
        list(method = "LFMM",  adjust = "bonf", threshold = "0.05"))))
    got <- default_threshold(cfg, MOD_GEA)
    expect_identical(got$type, "qval")
    expect_identical(got$value, 0.1)
})

test_that("default_threshold reads GWAS.configs for the GWAS module", {
    cfg <- list(
        GEA  = list(configs = list(list(method = "EMMAX", adjust = "qval", threshold = "0.1"))),
        GWAS = list(configs = list(list(method = "EMMAX", adjust = "top",  threshold = "500"))))
    expect_identical(default_threshold(cfg, MOD_GWAS)$type, "top")
    expect_identical(default_threshold(cfg, MOD_GWAS)$value, 500)
})

test_that("default_threshold falls back to bonf 0.05", {
    expect_identical(default_threshold(list(), MOD_GEA), list(type = "bonf", value = 0.05))
    expect_identical(default_threshold(list(GEA = list(configs = list())), MOD_GEA),
                     list(type = "bonf", value = 0.05))
})

test_that("default_threshold rejects a non-positive or unparseable threshold", {
    mk <- function(v) list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = v))))
    expect_identical(default_threshold(mk("0"), MOD_GEA)$value, 0.05)
    expect_identical(default_threshold(mk("-1"), MOD_GEA)$value, 0.05)
    expect_identical(default_threshold(mk("abc"), MOD_GEA)$value, 0.05)
})

test_that("threshold_value_default_for_type gives a sane default per mode", {
    expect_identical(threshold_value_default_for_type("bonf"), 0.05)
    expect_identical(threshold_value_default_for_type("qval"), 0.05)
    expect_identical(threshold_value_default_for_type("top"), 100)
    expect_identical(threshold_value_default_for_type("custom"), 1e-5)
    # An unknown mode must not return NULL: it feeds a numericInput's value.
    expect_identical(threshold_value_default_for_type("something_new"), 0.05)
})

test_that("threshold_value_valid_for_type bounds the probability modes to (0, 1]", {
    for (ty in c("bonf", "qval")) {
        expect_true(threshold_value_valid_for_type(ty, 0.05))
        expect_true(threshold_value_valid_for_type(ty, 1))
        expect_false(threshold_value_valid_for_type(ty, 0))
        expect_false(threshold_value_valid_for_type(ty, -0.01))
        # 100 as an alpha is the mistake this guard exists for: switching mode from
        # "top" to "bonf" leaves the old value in the input.
        expect_false(threshold_value_valid_for_type(ty, 100))
    }
})

test_that("threshold_value_valid_for_type requires top >= 1 and custom > 0", {
    expect_true(threshold_value_valid_for_type("top", 1))
    expect_true(threshold_value_valid_for_type("top", 5000))
    expect_false(threshold_value_valid_for_type("top", 0.5))
    expect_false(threshold_value_valid_for_type("top", 0))

    expect_true(threshold_value_valid_for_type("custom", 1e-8))
    # custom is a raw p cutoff and is deliberately NOT capped at 1.
    expect_true(threshold_value_valid_for_type("custom", 2))
    expect_false(threshold_value_valid_for_type("custom", 0))
})

test_that("threshold_value_valid_for_type rejects NULL, NA and non-numeric", {
    # The NULL check must come first: is.na(NULL) is logical(0), which would make
    # `if` error rather than return FALSE. Upstream, fct_threshold_rules.R:104-105
    # also collapses a length > 1 value to NA before calling in.
    expect_false(threshold_value_valid_for_type("bonf", NULL))
    expect_false(threshold_value_valid_for_type("bonf", NA))
    expect_false(threshold_value_valid_for_type("bonf", NA_real_))
    expect_false(threshold_value_valid_for_type("bonf", "abc"))
    # A numeric string IS accepted — the sidebar can hand over either.
    expect_true(threshold_value_valid_for_type("bonf", "0.05"))
})

test_that("threshold_value_valid_for_type rejects an unknown mode outright", {
    expect_false(threshold_value_valid_for_type("something_new", 0.05))
})

test_that("format_threshold_rule spells out the UNIT, not just the number", {
    expect_identical(format_threshold_rule("bonf", 0.05), "Bonferroni α 0.05")
    expect_identical(format_threshold_rule("qval", 0.05), "FDR q ≤ 0.05")
    expect_identical(format_threshold_rule("top", 100), "Top 100 SNPs/trait")
    expect_identical(format_threshold_rule("custom", 1e-5), "raw p < 1e-05")
})

test_that("format_threshold_rule rounds a fractional top-N count", {
    # "Top 99.6 SNPs/trait" is not a thing the pipeline can do.
    expect_identical(format_threshold_rule("top", 99.6), "Top 100 SNPs/trait")
})

test_that("format_threshold_rule degrades to the mode NAME with no value", {
    expect_identical(format_threshold_rule("bonf", NULL), "Bonferroni")
    expect_identical(format_threshold_rule("qval", NA), "FDR (qval)")
    # The unknown-mode-with-no-value case renders "NA" instead of the mode name:
    # the `%||% type` fallback at :290 is dead code (names(which(...))[1] is
    # NA_character_, not NULL). Filed; correct behaviour asserted in
    # test-known-bugs.R.
})

test_that("format_threshold_rule defaults a NULL type to bonf", {
    expect_identical(format_threshold_rule(NULL, 0.05), "Bonferroni α 0.05")
})

test_that("every THRESHOLD_TYPE_CHOICES mode is renderable and has a default", {
    # The choices vector drives the selectInput; a mode offered in the UI with no
    # validator branch or no default would be selectable and then broken.
    for (ty in unname(THRESHOLD_TYPE_CHOICES)) {
        expect_true(nzchar(format_threshold_rule(ty, 0.05)))
        dflt <- threshold_value_default_for_type(ty)
        expect_true(threshold_value_valid_for_type(ty, dflt),
                    info = paste0("mode ", ty, " rejects its own default ", dflt))
    }
})
