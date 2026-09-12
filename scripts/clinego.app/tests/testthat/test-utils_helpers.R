# R/utils_helpers.R — three small functions whose failure mode is silence.
#
# keep_merged_method_files() is the sharp one. The discovery layer globs
# methods/<M>/*_pvalues_K*.tsv, which matches BOTH the merged per-method table and
# the per-trait tables written beside it. It separates them by a path-token rule:
# the basename must start with the name of the `methods/` child directory. Get that
# wrong and the app reads a single trait's p-values as if they were the method's
# whole table — no error, a smaller SNP set, and nothing says which file was used.

mk_path <- function(...) paste(c(...), collapse = .Platform$file.sep)

test_that("keep_merged_method_files keeps the merged per-method table", {
    files <- c(
        mk_path("SIMDATA_results", "GEA", "tables", "methods", "EMMAX",
                "EMMAX_pvalues_K3.tsv"),
        mk_path("SIMDATA_results", "GEA", "tables", "methods", "LFMM",
                "LFMM_pvalues_K3.tsv"))
    expect_identical(keep_merged_method_files(files), files)
})

test_that("keep_merged_method_files drops the per-trait siblings", {
    dir <- c("SIMDATA_results", "GEA", "tables", "methods", "EMMAX")
    merged  <- mk_path(dir, "EMMAX_pvalues_K3.tsv")
    pertrait <- c(mk_path(dir, "bio_1_pvalues_K3.tsv"),
                  mk_path(dir, "height_pvalues_K3.tsv"))
    got <- keep_merged_method_files(c(pertrait[1], merged, pertrait[2]))
    expect_identical(got, merged)
})

test_that("keep_merged_method_files rejects a path with no methods/ segment", {
    # Without the segment there is no method name to compare the basename against,
    # so the only safe answer is "not a merged table".
    files <- mk_path("SIMDATA_results", "GEA", "tables", "EMMAX_pvalues_K3.tsv")
    expect_length(keep_merged_method_files(files), 0L)
})

test_that("keep_merged_method_files uses the WZA token in wza regime", {
    dir <- c("SIMDATA_results", "GEA", "tables", "methods", "EMMAX")
    snp <- mk_path(dir, "EMMAX_pvalues_K3.tsv")
    wza <- mk_path(dir, "EMMAX_wza_K3.tsv")

    expect_identical(keep_merged_method_files(c(snp, wza), regime = "snp"), snp)
    expect_identical(keep_merged_method_files(c(snp, wza), regime = "wza"), wza)
})

test_that("keep_merged_method_files handles an empty input", {
    expect_length(keep_merged_method_files(character(0)), 0L)
})

test_that("keep_merged_method_files is not fooled by a method-name PREFIX", {
    # A directory named EMMAX must not claim EMMAX2_pvalues_K3.tsv. startsWith on
    # "EMMAX" plus the token is what prevents it: the token has to follow the name
    # immediately.
    dir <- c("SIMDATA_results", "GEA", "tables", "methods", "EMMAX")
    expect_length(
        keep_merged_method_files(mk_path(dir, "EMMAX2_pvalues_K3.tsv")), 0L)
})

test_that("resolve_adjust finds the method's own rule", {
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"),
        list(method = "LFMM",  adjust = "qval", threshold = "0.1"))))
    expect_identical(resolve_adjust(cfg, "EMMAX"), "bonf_0.05")
    expect_identical(resolve_adjust(cfg, "LFMM"), "qval_0.1")
})

test_that("resolve_adjust returns NULL for a method that is not configured", {
    # NULL, not a default: the caller uses it to decide whether the method ran at
    # all, so inventing "bonf_0.05" here would fabricate an output path.
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))))
    expect_null(resolve_adjust(cfg, "RDA"))
    expect_null(resolve_adjust(list(), "EMMAX"))
})

test_that("resolve_adjust prefers the module's own configs block", {
    cfg <- list(
        GEA  = list(configs = list(list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))),
        GWAS = list(configs = list(list(method = "EMMAX", adjust = "top",  threshold = "500"))))
    expect_identical(resolve_adjust(cfg, "EMMAX", module = "GEA"),  "bonf_0.05")
    expect_identical(resolve_adjust(cfg, "EMMAX", module = "GWAS"), "top_500")
})

test_that("resolve_adjust falls back to GEA.configs when the module has none", {
    # The documented inheritance: GWAS inherits GEA's rules unless overridden.
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))))
    expect_identical(resolve_adjust(cfg, "EMMAX", module = "GWAS"), "bonf_0.05")
})

test_that("tab_to_mode maps each tab to its Snakemake mode string", {
    expect_identical(tab_to_mode("processing"), "processing")
    expect_identical(tab_to_mode("prestructure"), "prestructure")
    expect_identical(tab_to_mode("climate"), "climate")
    expect_identical(tab_to_mode("traits"), "traits")
    expect_identical(tab_to_mode("pregea"), "pregea")
    expect_identical(tab_to_mode("structure"), "structure")
    expect_identical(tab_to_mode("gea"), "gea")
    expect_identical(tab_to_mode("gwas"), "gwas")
})

test_that("tab_to_mode returns NULL for a tab with no single mode", {
    # Haplotype is driven by two separate runner instances, and an unknown tab must
    # not silently become a mode string that Snakemake would reject at parse time.
    expect_null(tab_to_mode("haplotype"))
    expect_null(tab_to_mode("not_a_tab"))
})

test_that("file_ok requires a single existing non-empty path", {
    d <- withr::local_tempdir()
    full  <- file.path(d, "full.tsv");  writeLines("x", full)
    empty <- file.path(d, "empty.tsv"); file.create(empty)

    expect_true(file_ok(full))
    expect_false(file_ok(empty))
    expect_false(file_ok(file.path(d, "missing.tsv")))
    expect_false(file_ok(NULL))
    expect_false(file_ok(NA_character_))
    expect_false(file_ok(c(full, full)))
})
