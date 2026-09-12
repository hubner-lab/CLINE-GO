# Correct-behaviour assertions for defects in the Shiny app that are KNOWN and
# DELIBERATELY NOT FIXED yet.
#
# Same convention as the pipeline root's tests/testthat/test-known-bugs.R: every
# block starts with skip(), so the suite stays green and usable as a merge gate,
# and every block is written as the behaviour the code SHOULD have — never as an
# assertion of the buggy output. Fixing a bug means deleting one skip() line and
# getting a real regression test for free. Do not "fix" a test here by weakening it
# to match current output.
#
# Each skip() names where the defect is filed: docs/pipeline_improvement_requests.md
# and/or the ADAPTOGENE testing track's ## Findings.
#
# NOTE ON SCOPE. Most app-side quarantines live in the ROOT suite's
# tests/testthat/test-equivalence-app-pipeline.R, because they are divergences
# between the app and the pipeline (load_gff_genes' mixed-attribute GFF, the
# Cross-method min_pvalue shift, fct_regions.R:51's dropped setorder, the two
# meanings of promoter_length, qvalue-not-in-Imports). This file is only for
# defects internal to the app, where there is no pipeline counterpart to compare
# against.
#
# DELIBERATELY NOT HERE: .apply_bounds()'s silent `union` fallback for an
# unrecognised strategy (R/fct_overlap.R:119-122). test-fct_overlap.R:10-14 already
# pins it as a characterization contract, and
# docs/pipeline_improvement_requests.md:218 records that pinning as the AGREED
# mitigation. A skip()ped counter-assertion here would contradict an active test in
# a sibling file and re-open a closed decision.

# ---------------------------------------------------------------------------
# 1. fct_combine.R:290 — format_threshold_rule()'s no-value branch ends in
#    `names(which(THRESHOLD_TYPE_CHOICES == type))[1] %||% type`. For a mode that
#    is not in THRESHOLD_TYPE_CHOICES, which() is empty, names() is character(0),
#    and [1] is NA_character_ — NOT NULL — so %||% never fires and the `type`
#    fallback is unreachable dead code. The summary badge then reads "NA".
#    Reachable from config: a method's `adjust` string comes from the project YAML
#    (or a registry adjust_default), so a typo like `fdr` produces it.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-12).
#    Fix: use a length check, e.g.
#      nm <- names(which(THRESHOLD_TYPE_CHOICES == type)); if (length(nm)) nm[1] else type
# ---------------------------------------------------------------------------

test_that("format_threshold_rule names an unrecognised mode instead of NA", {
    skip("known bug: fct_combine.R:290 %||% cannot catch NA_character_ — filed 2026-09-12")

    expect_identical(format_threshold_rule("something_new", NULL), "something_new")
    expect_identical(format_threshold_rule("fdr", NA), "fdr")
})

# ---------------------------------------------------------------------------
# 2. DESCRIPTION — `sass` is used ::-qualified at R/app_theme.R:32,35 but appears
#    in neither Imports nor Suggests; it arrives transitively via bslib.
#
#    This is NOT the same mechanism as the filed qvalue defect, and the difference
#    matters: a ::-qualified call resolves whenever the package is INSTALLED,
#    regardless of DESCRIPTION, so there is no wrong behaviour today (sass 0.4.10
#    is present in the image). It becomes real only on a pruned install — which is
#    the mode the Dockerfile uses (dependencies = FALSE). So the defect is the
#    undeclared dependency itself, pinned structurally, with no behavioural
#    counterpart to assert.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-12). Severity: low.
# ---------------------------------------------------------------------------

test_that("sass is declared in Imports, since app_theme.R calls it", {
    skip("known bug: sass used at app_theme.R:32,35 but undeclared — filed 2026-09-12")

    imports <- utils::packageDescription("clinego.app", fields = "Imports")
    expect_true(grepl("\\bsass\\b", imports))
})
