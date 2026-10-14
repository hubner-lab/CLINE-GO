# scripts/build_snp_set.R — the wrapper that builds a config-declared SNP set
# (GEA.snp_sets / GWAS.snp_sets) into the store mode=maladaptation reads.
#
# The exit-0 smoke row lives in test-cli-wrappers.R. This file runs the script
# as a subprocess and asserts CONTENT: that the file it writes is byte-for-byte
# what select_snp_set() returns in-process (so the wrapper adds no selection
# logic of its own), the inclusive cutoff, the zero-SNP contract (header-only TSV,
# n_snps: 0, exit 0), the loud failures, and that a record it wrote is accepted
# back as input — the mechanism the Shiny round trip will use.
#
# Fixture: fx_snp_set_inputs() in tests/lib/wrapper_harness.R.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

build_set <- function(f, decl = f$decl, tsv = f$tsv, record = f$record,
                      pvals = f$pvals, d = dirname(f$decl)) {
    sandbox <- file.path(d, "cwd")
    dir.create(sandbox, showWarnings = FALSE)
    res <- run_wrapper("build_snp_set.R", c(decl, tsv, record, "3", pvals), wd = sandbox)
    expect_identical(list.files(sandbox, all.files = TRUE, no.. = TRUE), character(0))
    res
}

read_set <- function(path) {
    data.table::fread(path, sep = "\t", colClasses = c(SNPID = "character", chr = "character"))
}

test_that("the TSV is exactly what select_snp_set() writes in-process", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, strategy = "Cross-method per-trait")
    res <- build_set(f)
    expect_identical(res$status, 0L, info = res$output)

    decl  <- read_snp_set_decl(f$decl)
    pvals <- lapply(f$pvals, data.table::fread, sep = "\t",
                    colClasses = c(chr = "character", SNPID = "character"))
    names(pvals) <- basename(dirname(f$pvals))
    want <- quiet(select_snp_set(pvals, decl$cells, decl$pairs,
                                 decl$strategy, decl$clumping_distance))$snps
    ref <- file.path(d, "ref.tsv")
    data.table::fwrite(want, ref, sep = "\t", quote = FALSE)
    expect_identical(readLines(f$tsv), readLines(ref))
    expect_identical(read_set(f$tsv)$SNPID, c("1:100", "1:200"))
})

test_that("the cutoff is inclusive", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, pairs = "bio_1::EMMAX")
    expect_identical(build_set(f)$status, 0L)
    expect_identical(read_set(f$tsv)$SNPID, c("1:100", "2:300"))   # 2:300 p = 1e-4
})

test_that("the record carries the declaration plus k_best, n_snps and created", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d)
    expect_identical(build_set(f)$status, 0L)
    rec <- read_snp_set_decl(f$record)
    decl <- read_snp_set_decl(f$decl)
    expect_identical(rec[SNP_SET_DECL_FIELDS], decl[SNP_SET_DECL_FIELDS])
    expect_identical(rec$k_best, 3L)
    expect_identical(rec$n_snps, 3L)
    expect_match(rec$created, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}$")
})

test_that("zero selected SNPs: header-only TSV, n_snps 0, a warning, exit 0", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, value = 1e-30)
    res <- build_set(f)
    expect_identical(res$status, 0L, info = res$output)
    expect_identical(readLines(f$tsv), "SNPID\tchr\tpos\tmin_pvalue")
    expect_identical(read_snp_set_decl(f$record)$n_snps, 0L)
    expect_match(res$output, "WARNING: SNP set 'gea.fx' selected 0 SNPs")
})

test_that("a record it wrote is accepted back as input and rebuilds the same set", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, strategy = "Cross-method")
    expect_identical(build_set(f)$status, 0L)
    tsv2 <- file.path(d, "again", "selected_snps.tsv")
    rec2 <- file.path(d, "again", "snp_set.yaml")
    res <- build_set(f, decl = f$record, tsv = tsv2, record = rec2)
    expect_identical(res$status, 0L, info = res$output)
    expect_identical(readLines(tsv2), readLines(f$tsv))
})

test_that("a pair naming a missing trait column stops, writing nothing", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, pairs = c("bio_1::EMMAX", "bio_9::LFMM"))
    res <- build_set(f)
    expect_gt(res$status, 0L)
    expect_match(res$output, "trait column missing for pair\\(s\\): bio_9::LFMM")
    expect_false(file.exists(f$tsv))
    expect_false(file.exists(f$record))
})

test_that("a method with no p-value file stops", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, methods = "EMMAX")
    res <- build_set(f)
    expect_gt(res$status, 0L)
    expect_match(res$output, "names method\\(s\\) with no p-value file: LFMM")
})

test_that("a non-snp regime is refused", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d, regime = "wza")
    res <- build_set(f)
    expect_gt(res$status, 0L)
    expect_match(res$output, "regime 'wza' cannot be built from config")
})

test_that("too few arguments print the usage and exit non-zero", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_snp_set_inputs(d)
    res <- run_wrapper("build_snp_set.R", c(f$decl, f$tsv, f$record, "3"))
    expect_gt(res$status, 0L)
    expect_match(res$output, "Usage: build_snp_set.R")
})
