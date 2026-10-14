# filter_imputed_maf.R — CONTENT assertions, quick tier.
#
# The script decides which SNPs every GEA/GWAS method tests (sNMF.imputed_maf on the
# imputed matrix), and it writes three files that downstream rules pair BY POSITION:
# column j of the .lfmm is line j of the .vcfsnp is body row j of the .vcf. A smoke row
# (exit 0, non-empty) cannot see a wrong keep set, a misaligned triple or a boundary SNP
# lost to floating point, so these read the outputs back.
#
# Every failure case also asserts NO output is left behind: Snakemake would delete a
# failed rule's outputs, but a run by hand (or a later caller) must not find a partial
# triple that looks finished.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

out_paths <- function(d) {
    c(lfmm   = file.path(d, "out.lfmm"),
      vcfsnp = file.path(d, "out.vcfsnp"),
      vcf    = file.path(d, "out.vcf"),
      report = file.path(d, "out_report.tsv"))
}

run_filter <- function(d, f, threshold) {
    o <- out_paths(d)
    res <- run_wrapper("filter_imputed_maf.R",
                       c(f$lfmm, f$vcfsnp, f$vcf, as.character(threshold),
                         o[["lfmm"]], o[["vcfsnp"]], o[["vcf"]], o[["report"]]))
    c(res, list(out = o))
}

read_lfmm <- function(p) {
    as.matrix(utils::read.table(p, header = FALSE, colClasses = "integer"))
}

expect_no_outputs <- function(r) {
    for (p in r$out) {
        expect_false(file.exists(p), info = paste0("failed run left ", p, "\n", r$output))
    }
}

test_that("keeps exactly the SNPs at or above the threshold, as one aligned triple", {
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d)
    r <- run_filter(d, f, 0.14)
    expect_identical(r$status, 0L, info = r$output)

    keep <- c(1L, 2L, 4L, 6L)
    m <- read_lfmm(r$out[["lfmm"]])
    dimnames(m) <- NULL
    expect_identical(m, f$g[, keep])

    expect_identical(readLines(r$out[["vcfsnp"]]), readLines(f$vcfsnp)[keep])

    vin  <- readLines(f$vcf)
    vout <- readLines(r$out[["vcf"]])
    head_in <- vin[startsWith(vin, "#")]
    expect_identical(vout[startsWith(vout, "#")], head_in)
    expect_identical(vout[!startsWith(vout, "#")], vin[!startsWith(vin, "#")][keep])

    rep <- data.table::fread(r$out[["report"]])
    expect_identical(names(rep), c("threshold", "n_samples", "n_snps_before",
                                   "n_snps_after", "n_snps_removed"))
    expect_equal(rep$threshold, 0.14)
    expect_equal(c(rep$n_samples, rep$n_snps_before, rep$n_snps_after, rep$n_snps_removed),
                 c(25, 6, 4, 2))
})

test_that("a SNP exactly on the threshold is kept despite floating point", {
    # 0.14 * 50 == 7.000000000000001 in double precision: mac 7 is MAF 0.14 exactly.
    expect_gt(0.14 * 50, 7)
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d, ac = c(7L, 43L, 6L))
    r <- run_filter(d, f, 0.14)
    expect_identical(r$status, 0L, info = r$output)
    expect_identical(ncol(read_lfmm(r$out[["lfmm"]])), 2L)
    expect_identical(readLines(r$out[["vcfsnp"]]), readLines(f$vcfsnp)[1:2])
})

test_that("minor allele is counted from whichever side is rarer", {
    # ac 44 of 50 = the COUNTED allele is the major one; mac is 6, MAF 0.12.
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d, ac = c(44L, 25L))
    r <- run_filter(d, f, 0.13)
    expect_identical(r$status, 0L, info = r$output)
    expect_identical(readLines(r$out[["vcfsnp"]]), readLines(f$vcfsnp)[2])
})

test_that("threshold 0 copies all three inputs byte for byte", {
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d, ac = c(0L, 1L, 25L, 50L, 7L))   # incl. both monomorphic shapes
    r <- run_filter(d, f, 0)
    expect_identical(r$status, 0L, info = r$output)
    md5 <- function(p) unname(tools::md5sum(p))
    expect_identical(md5(r$out[["lfmm"]]),   md5(f$lfmm))
    expect_identical(md5(r$out[["vcfsnp"]]), md5(f$vcfsnp))
    expect_identical(md5(r$out[["vcf"]]),    md5(f$vcf))
    rep <- data.table::fread(r$out[["report"]])
    expect_equal(c(rep$n_snps_before, rep$n_snps_after, rep$n_snps_removed), c(5, 5, 0))
})

test_that("removing every SNP fails with no output", {
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d, ac = c(1L, 2L, 49L))
    r <- run_filter(d, f, 0.3)
    expect_gt(r$status, 0L)
    expect_match(r$output, "no SNP has imputed MAF")
    expect_no_outputs(r)
})

test_that("malformed or misaligned inputs fail with no output", {
    cases <- list(
        "ragged lfmm row" = function(f) {
            x <- readLines(f$lfmm); x[3] <- sub(" [0-9]$", "", x[3]); writeLines(x, f$lfmm)
        },
        "missing call (9) left in the lfmm" = function(f) {
            x <- readLines(f$lfmm); x[2] <- sub("^[0-9]", "9", x[2]); writeLines(x, f$lfmm)
        },
        "dosage above 2" = function(f) {
            x <- readLines(f$lfmm); x[2] <- sub("^[0-9]", "3", x[2]); writeLines(x, f$lfmm)
        },
        "vcfsnp one line short" = function(f) {
            x <- readLines(f$vcfsnp); writeLines(x[-length(x)], f$vcfsnp)
        },
        "VCF POS disagrees with vcfsnp" = function(f) {
            x <- readLines(f$vcf); i <- length(x) - 1L
            x[i] <- sub("^([^\t]+)\t[0-9]+", "\\1\t999999", x[i]); writeLines(x, f$vcf)
        },
        "VCF one body row short" = function(f) {
            x <- readLines(f$vcf); writeLines(x[-length(x)], f$vcf)
        },
        "VCF one body row extra" = function(f) {
            x <- readLines(f$vcf); writeLines(c(x, x[length(x)]), f$vcf)
        }
    )
    for (label in names(cases)) {
        d <- withr::local_tempdir()
        f <- fx_imputed_triple(d)
        cases[[label]](f)
        r <- run_filter(d, f, 0.14)
        expect_gt(r$status, 0L, label = label)
        for (p in r$out) {
            expect_false(file.exists(p), info = paste0(label, ": left ", p, "\n", r$output))
        }
    }
})

test_that("bad argv fails before reading anything", {
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d)
    o <- out_paths(d)
    short <- run_wrapper("filter_imputed_maf.R", c(f$lfmm, f$vcfsnp, f$vcf, "0.1"))
    expect_gt(short$status, 0L)
    expect_match(short$output, "expected 8 arguments")
    for (bad in c("0.5", "-0.1", "abc")) {
        r <- run_wrapper("filter_imputed_maf.R",
                         c(f$lfmm, f$vcfsnp, f$vcf, bad, o[["lfmm"]], o[["vcfsnp"]],
                           o[["vcf"]], o[["report"]]))
        expect_gt(r$status, 0L, label = bad)
        expect_false(any(file.exists(o)), label = bad)
    }
})

test_that("write_summary.R turns the report into the three summary rows", {
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d)
    r <- run_filter(d, f, 0.14)
    expect_identical(r$status, 0L, info = r$output)
    summ <- file.path(d, "pipeline_summary.tsv")
    ld <- fx_ld_decay_table(d)
    res <- run_wrapper("write_summary.R",
                       c("structure", summ, "3", "NULL", "NULL", ld, "site", "both", "NULL",
                         r$out[["report"]]))
    expect_identical(res$status, 0L, info = res$output)
    s <- data.table::fread(summ, colClasses = "character")
    got <- s[step == "structure" & grepl("imputed_maf", metric)]
    expect_identical(setNames(got$value, got$metric),
                     c(imputed_maf_threshold = "0.14", snps_before_imputed_maf = "6",
                       snps_after_imputed_maf = "4"))
})
