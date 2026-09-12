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

# ═══════════════════════════════════════════════════════════════════════════════
# compute_interactive_sigsnps() — fct_combine.R:690
#
# The trait x method matrix's selection arrives as a JSON string from JS, and the
# contract at :706-708 is that an EMPTY array means "the matrix rendered with no
# active cells", i.e. un-initialised, NOT a deliberate deselect. So NULL, "",
# "[]" and unparseable JSON all mean SHOW EVERYTHING — there is no way to express
# "select nothing" through this argument at all. That is the kind of rule that
# inverts silently under a refactor, so every one of those four inputs is pinned.
#
# The other sharp edge is the return type, which is three different things:
#   NULL                      when the method list is empty (:694)
#   7-column empty table      when nothing survives filtering (:731, :736)
#   8 columns on success      combine_sigsnps stamps min_pvalue at :135
# save_snp_set() consumes this and needs min_pvalue, so the difference is load
# bearing rather than cosmetic.
# ═══════════════════════════════════════════════════════════════════════════════

cis_methods <- function() {
    list(
        EMMAX = data.table::data.table(
            SNPID = c("1:100", "2:300"), chr = c("1", "2"), pos = c(100L, 300L),
            pvalue = c(1e-8, 5e-7), method = "EMMAX", trait = c("bio_1", "bio_2")),
        LFMM = data.table::data.table(
            SNPID = "1:100", chr = "1", pos = 100L,
            pvalue = 3e-9, method = "LFMM", trait = "bio_1")
    )
}

cis_counts <- function() list(`bio_1::EMMAX` = 1L, `bio_2::EMMAX` = 1L, `bio_1::LFMM` = 1L)

call_cis <- function(tm_json, strategy = "Union", counts = cis_counts(),
                     known = c("bio_1", "bio_2")) {
    compute_interactive_sigsnps(
        all_method_sigsnps = cis_methods(), tm_selection_json = tm_json,
        combo_counts = counts, known_traits = known, strategy = strategy,
        clumping_distance = 1000L, project_name = "P", module = MOD_GEA)
}

test_that("compute_interactive_sigsnps returns NULL only for an empty method list", {
    expect_null(compute_interactive_sigsnps(
        all_method_sigsnps = list(), tm_selection_json = NULL,
        combo_counts = list(), known_traits = character(0), strategy = "Union",
        clumping_distance = 1000L, project_name = "P", module = MOD_GEA))
})

test_that("NULL, \"\", \"[]\" and unparseable JSON all mean SHOW EVERYTHING", {
    # The line that makes "[]" mean ALL rather than none is :716, and the
    # tryCatch at :712-715 sends malformed JSON down the same path.
    baseline <- call_cis(NULL)
    expect_identical(nrow(baseline), 3L)
    for (json in list("", "[]", "not json at all", "{oops")) {
        got <- call_cis(json)
        expect_identical(nrow(got), nrow(baseline), info = paste("input:", json))
        expect_setequal(got$SNPID, baseline$SNPID)
    }
})

test_that("there is no way to express \"select nothing\" through tm_selection_json", {
    # A direct consequence of the rule above, stated as its own assertion because
    # it is the surprising half: a caller cannot deselect everything.
    expect_identical(nrow(call_cis("[]")), 3L)
})

test_that("an explicit selection narrows to the named trait::method pairs", {
    got <- call_cis(jsonlite::toJSON(c("bio_1::EMMAX")))
    expect_identical(unique(got$method), "EMMAX")
    expect_identical(unique(got$trait), "bio_1")
})

test_that("a selection naming one method drops the other method entirely", {
    got <- call_cis(jsonlite::toJSON(c("bio_1::LFMM")))
    expect_identical(unique(got$method), "LFMM")
})

test_that("combo_counts of zero and unknown traits are both excluded from the default", {
    # Two separate filters at :698-703: a non-positive count, and a trait absent
    # from known_traits (which prevents cross-tab bleed).
    got <- call_cis(NULL, counts = list(`bio_1::EMMAX` = 1L, `bio_2::EMMAX` = 0L,
                                        `bio_1::LFMM` = 1L))
    expect_false("bio_2" %in% got$trait)

    got2 <- call_cis(NULL, known = "bio_1")
    expect_false("bio_2" %in% got2$trait)
})

test_that("compute_interactive_sigsnps returns the 7-column empty table when nothing survives", {
    got <- call_cis(jsonlite::toJSON(c("bio_9::EMMAX")))
    expect_identical(nrow(got), 0L)
    expect_identical(names(got), c("SNPID", "chr", "pos", "pvalue",
                                   "method", "trait", "region_id"))
    expect_false("min_pvalue" %in% names(got))
})

test_that("the success shape carries min_pvalue, which the empty shape does not", {
    # save_snp_set() aggregates min(min_pvalue), so a caller handed the empty
    # shape errors rather than writing an empty set.
    got <- call_cis(NULL)
    expect_true("min_pvalue" %in% names(got))
    expect_identical(got[SNPID == "1:100"][1L]$min_pvalue, 3e-9)
})

test_that("region_id is returned as NA, to be stamped from live regions downstream", {
    got <- call_cis(NULL)
    expect_true(all(is.na(got$region_id)))
})

test_that("legacy strategy aliases are normalised before combining", {
    for (alias in c("All", "Sum")) {
        expect_identical(nrow(call_cis(NULL, strategy = alias)), 3L)
    }
})

test_that("tm_selection_json = character(0) errors rather than falling back", {
    # nzchar() at :709 is not length-guarded, so a zero-length argument is an
    # error where every other degenerate input silently means ALL.
    expect_error(call_cis(character(0)))
})

# ═══════════════════════════════════════════════════════════════════════════════
# default_threshold() vs resolve_adjust() — fct_combine.R:202 and
# utils_helpers.R:61.
#
# Both read the SAME config block, {GEA|GWAS}.configs, by two different rules:
# default_threshold takes configs[[1]] and is method-agnostic; resolve_adjust
# matches on cfg$method. The pipeline has no first-entry notion at all —
# common.smk:584-598 builds a strictly per-method dict and raises on duplicate
# methods — so resolve_adjust agrees with it and default_threshold does not.
#
# Consequence on one render: the app opens a method's p-value file at the
# threshold resolve_adjust returns, while seeding the interactive threshold
# control from the FIRST config entry.
# ═══════════════════════════════════════════════════════════════════════════════

two_method_config <- function() {
    list(GEA = list(configs = list(
        list(method = "LFMM",  adjust = "qval", threshold = 0.1),
        list(method = "EMMAX", adjust = "bonf", threshold = 0.05))))
}

test_that("resolve_adjust matches on method, as the pipeline's per-method dict does", {
    cfg <- two_method_config()
    expect_identical(resolve_adjust(cfg, "LFMM",  "GEA"), "qval_0.1")
    expect_identical(resolve_adjust(cfg, "EMMAX", "GEA"), "bonf_0.05")
})

test_that("resolve_adjust returns NULL for a method absent from the configs", {
    expect_null(resolve_adjust(two_method_config(), "RDA", "GEA"))
})

test_that("default_threshold ignores the method and always takes the FIRST entry", {
    # The divergence, pinned. For this config the app opens EMMAX's file at
    # bonf_0.05 while seeding the threshold control with qval/0.1.
    got <- default_threshold(two_method_config(), MOD_GEA)
    expect_identical(got$type, "qval")
    expect_identical(got$value, 0.1)
    expect_false(identical(got$type, "bonf"))
})

test_that("default_threshold falls back to bonf/0.05 when the module has no configs", {
    expect_identical(default_threshold(list(), MOD_GEA), list(type = "bonf", value = 0.05))
})

test_that("the two helpers disagree on an empty GWAS block: GEA fallback vs hardcoded default", {
    # resolve_adjust falls back to GEA.configs (utils_helpers.R:65);
    # default_threshold does not (fct_combine.R:205). With GWAS.configs empty and
    # GEA populated, the path and the threshold come from different sources.
    cfg <- two_method_config()
    expect_identical(resolve_adjust(cfg, "LFMM", "GWAS"), "qval_0.1")
    expect_identical(default_threshold(cfg, MOD_GWAS), list(type = "bonf", value = 0.05))
})

test_that("default_threshold repairs a non-numeric or non-positive threshold", {
    bad <- list(GEA = list(configs = list(list(method = "LFMM", adjust = "qval",
                                               threshold = "not a number"))))
    expect_identical(default_threshold(bad, MOD_GEA)$value, 0.05)
    neg <- list(GEA = list(configs = list(list(method = "LFMM", adjust = "qval",
                                               threshold = -1))))
    expect_identical(default_threshold(neg, MOD_GEA)$value, 0.05)
})

test_that("MOD_GEAXGWAS resolves to the GEA block on both paths", {
    cfg <- two_method_config()
    expect_identical(default_threshold(cfg, MOD_GEAXGWAS)$type, "qval")
    expect_identical(resolve_adjust(cfg, "EMMAX", MOD_GEAXGWAS), "bonf_0.05")
})
