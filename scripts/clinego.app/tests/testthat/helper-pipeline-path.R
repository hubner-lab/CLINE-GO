# Redirect every {PROJECT}_results/ write this suite makes into a session tempdir.
#
# THE DEFECT THIS FIXES. Before this helper, the app suite left
# test_threshold_rules_<random>_results/ directories in the REPO ROOT, dated
# 2026-07-28 and accumulating invisibly (.gitignore:23 `*_results/` hides them, so
# `git status` stayed clean). Mechanism: compute_method_sigsnps_cached() reaches
# fct_data_loading.R:909 `dir.create(cache_dir, recursive = TRUE)`, where
# cache_dir resolves through interactive_sigsnps_dir() -> mod_path() ->
# project_base() -> get_pipeline_path(). recursive = TRUE creates the WHOLE chain
# including <project>_results/, while the tests' own cleanup
# (test-fct_threshold_rules.R:197-199,213-216) unlinks only the leaf .../GEA.
#
# WHY ONE OPTION IS ENOUGH. get_pipeline_path() (R/app_config.R:33-42) checks
# getOption("clinego.pipeline_path") FIRST, ahead of the golem config
# (inst/golem-config.yml:5 = /pipeline = the bind-mounted repo), the PIPELINE_PATH
# env var, and the "/pipeline" fallback. Every dir-creating site in the package
# goes through project_base(), so redirecting that one option covers
# fct_data_loading.R:746,806,909, fct_snp_sets.R:59-60,137,
# fct_region_params.R:39 and fct_pipeline.R:50 at once — and unblocks testing the
# file-IO functions at all, which is the same fix.
#
# plain options(), NOT withr::local_options(). A withr::local_* call at the top
# level of a helper unwinds the moment the helper finishes sourcing, so it would
# look correct and do nothing.
#
# NOT touched: CLINEGO_SHARED_LIBS (R/zzz.R:39-42) and app_sys()
# (R/app_config.R:11) hardcode /pipeline and never consult get_pipeline_path().
# .onLoad() has already run at library(clinego.app) in testthat.R before any
# helper is sourced, so shared-lib resolution is unaffected.

CLINEGO_TEST_ROOT <- file.path(tempdir(), "clinego-test-root")
dir.create(CLINEGO_TEST_ROOT, recursive = TRUE, showWarnings = FALSE)
options(clinego.pipeline_path = CLINEGO_TEST_ROOT)

# Reset the package's in-memory cache.
#
# .cache_env (fct_data_loading.R:1068) is a 200MB cachem::cache_mem reached via
# get_app_cache() and NEVER cleared, so every load_cached()-wrapped function
# carries state across tests within a session. Two keys are especially sharp:
#   * load_gff_genes()      — key "gff_genes_<project>" with NO fingerprint, so it
#                             is sticky for the whole session;
#   * assign_region_ids()   — fingerprinted on the regions file's MTIME, so two
#                             writes inside the same second collide.
# The standing discipline is therefore a UNIQUE `project` string per test block;
# use this only where a genuinely cold cache is required, via
# withr::defer(reset_app_cache()).
reset_app_cache <- function() {
    cache <- get_app_cache()
    if (!is.null(cache)) cache$reset()
    invisible(NULL)
}
