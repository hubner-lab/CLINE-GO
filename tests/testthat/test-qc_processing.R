# scripts/R/lib/qc_processing.R
#
# The numbers behind Processing's QC tables, taken from the step that decided them (release
# Phase 4c): plink's own report of filter_vcf's --geno and --maf steps (7e11b6), the lists the
# sample filters wrote and the het filter's own bounds (7ccf66), and missing calls and depth
# per sampling site (787ffe). plot_qc_processing.R calls all of it.

# plink 1.90b7.2's log of filter_vcf's flags on SIMDATA's raw VCF and 47 final samples, from
# the 4c probe: --maf 0.3 --geno 0.02, and --maf 0.05 --geno 1 (plink then skips --geno and
# prints no line for it). Progress and header lines trimmed; the parser must ignore the rest.
plink_log_strict <- c(
    "PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/",
    "354 variants loaded from .bim file.",
    "51 samples (0 males, 0 females, 51 ambiguous) loaded from .fam.",
    "--keep: 47 samples remaining.",
    "Total genotyping rate in remaining samples is 0.960332.",
    "300 variants removed due to missing genotype data (--geno).",
    "1 variant removed due to minor allele threshold(s)",
    "(--maf/--max-maf/--mac/--max-mac).",
    "53 variants and 47 samples pass filters and QC.",
    "Note: No phenotypes present.")
plink_log_geno1 <- c(
    "354 variants loaded from .bim file.",
    "--keep: 47 samples remaining.",
    "0 variants removed due to minor allele threshold(s)",
    "(--maf/--max-maf/--mac/--max-mac).",
    "354 variants and 47 samples pass filters and QC.")

test_that("plink's filter report is read step by step, --geno before --maf", {
    p <- parse_plink_filter_log(plink_log_strict, geno = 0.02)
    expect_identical(p, list(n_loaded = 354, n_removed_geno = 300, n_removed_maf = 1,
                             n_pass = 53, n_samples_pass = 47))
})

test_that("a missing --geno line reads as 0 only at --geno 1, where plink skips the filter", {
    expect_identical(parse_plink_filter_log(plink_log_geno1, geno = 1)$n_removed_geno, 0)
    expect_error(parse_plink_filter_log(plink_log_geno1, geno = 0.2), "no --geno line")
})

test_that("a log that is not a filtering run, or does not add up, stops", {
    no_maf <- plink_log_strict[!grepl("minor allele", plink_log_strict)]
    expect_error(parse_plink_filter_log(no_maf, geno = 0.02), "no --maf line")
    expect_error(parse_plink_filter_log(c("PLINK v1.90b7.2", "Error: something"), geno = 0.2),
                 "pass filters and QC")
    off <- sub("^53 variants", "54 variants", plink_log_strict)
    expect_error(parse_plink_filter_log(off, geno = 0.02), "does not add up")
})

test_that("older plink wording ('people') is read too", {
    old <- sub("47 samples pass", "47 people pass", plink_log_strict)
    expect_identical(parse_plink_filter_log(old, geno = 0.02)$n_samples_pass, 47)
})

test_that("filtering_attrition rows are the steps that ran, and every count telescopes", {
    p <- parse_plink_filter_log(plink_log_strict, geno = 0.02)
    a <- filtering_attrition(51, 354,
                             list("After sample missingness filter" = 50L,
                                  "After relatedness filter" = 47L),
                             p, 30)
    expect_identical(as.character(a$stage),
                     c("Raw VCF", "After sample missingness filter", "After relatedness filter",
                       "After SNP missingness filter", "After MAF + SNP missingness filter",
                       "After LD pruning"))
    expect_identical(a$n_samples, c(51, 50, 47, 47, 47, 47))
    # 354 - 300 (--geno) = 54, then - 1 (--maf) = 53: what plink reported, in its order
    expect_identical(a$n_snps, c(354, 354, 354, 54, 53, 30))
    expect_identical(a$pct_snps, round(100 * c(354, 354, 354, 54, 53, 30) / 354, 1))
    expect_identical(nrow(check_filtering_monotone(a)), 0L)
})

test_that("het_outlier_bounds uses the het filter's population SD", {
    x <- c(0.30, 0.31, 0.32, 0.33, 0.50)
    b <- het_outlier_bounds(x, 2)
    pop_sd <- sqrt(mean((x - mean(x))^2))
    expect_equal(unname(b[c("mean", "sd")]), c(mean(x), pop_sd))
    expect_equal(unname(b[c("lo", "hi")]), mean(x) + c(-2, 2) * pop_sd)
    expect_false(isTRUE(all.equal(pop_sd, sd(x))))   # sd() divides by n - 1
    # 0.50 lies beyond 1.9 population SDs but inside 1.9 sample SDs: the formula decides
    expect_true(0.50 > het_outlier_bounds(x, 1.9)[["hi"]])
    expect_true(0.50 < mean(x) + 1.9 * sd(x))
})

test_that("read_sample_list takes the IID of space- and tab-separated FID IID lines", {
    path <- withr::local_tempfile()
    writeLines(c("0 NEG01", "0\tTAV02", "", "0 88"), path)
    expect_identical(read_sample_list(path), c("NEG01", "TAV02", "88"))
    writeLines(c("0 NEG01", "NEG02"), path)
    expect_error(read_sample_list(path), "without an IID")
})

# 3 sites: A loses one sample to the missingness filter, B loses one later (relatedness), C a
# whole site of weak samples, every one dropped. Numeric-looking IDs stay character.
site_meta  <- data.table::data.table(site = c("A", "A", "A", "B", "B", "C", "C"),
                                     sample = c("a1", "a2", "a3", "88", "b2", "c1", "c2"))
site_smiss <- data.table::data.table(IID = c("c2", "c1", "b2", "88", "a3", "a2", "a1"),
                                     F_MISS = c(0.9, 0.8, 0.1, 0.2, 0.7, 0.1, 0.1))
passed <- c("a1", "a2", "88", "b2")
final  <- c("a1", "a2", "88")

test_that("site_qc_table lists every site, dropped samples included, in the metadata's order", {
    s <- site_qc_table(site_meta, site_smiss, passed, final)
    expect_identical(s$site, c("A", "B", "C"))
    expect_identical(s$n_samples, c(3L, 2L, 2L))
    expect_identical(s$n_failed_missingness, c(1L, 0L, 2L))
    expect_identical(s$n_retained, c(2L, 1L, 0L))
    expect_identical(s$f_miss_mean, round(c(mean(c(0.1, 0.1, 0.7)), 0.15, 0.85), 4))
    expect_identical(s$f_miss_median, c(0.1, 0.15, 0.85))
    expect_identical(s$f_miss_max, c(0.7, 0.2, 0.9))
    expect_true(all(is.na(s$depth_mean)) && all(is.na(s$depth_median)))
})

test_that("site_qc_table summarises per-sample depth, leaving out samples without one", {
    depth <- data.table::data.table(INDV = c("a1", "a2", "a3", "88", "b2", "c1"),
                                    MEAN_DEPTH = c(10, 20, NaN, 6, 8, 1.5))
    s <- site_qc_table(site_meta, site_smiss, passed, final, depth = depth)
    expect_identical(s$depth_mean, c(15, 7, 1.5))
    expect_identical(s$depth_median, c(15, 7, 1.5))
})

test_that("a sample of the missingness table without a site stops the table", {
    stray <- rbind(site_smiss, data.table::data.table(IID = "zz9", F_MISS = 0.1))
    expect_error(site_qc_table(site_meta, stray, passed, final), "no site in the metadata: zz9")
})
