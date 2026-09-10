# scripts/R/lib/sig_snps.R
#
# Per-trait significance filtering plus cross-trait overlap annotation. The
# selection rule is INCLUSIVE (`p <= threshold`, line 67) because
# compute_pval_threshold() returns a cutoff that is itself a member of the
# intended call set for 'top' and 'qval'.

# --- fixtures --------------------------------------------------------------

pvals <- function(chr = "1", pos = c(100, 200, 300),
                  SNPID = paste0(chr, ":", pos), ...) {
    dt <- data.table::data.table(SNPID = SNPID, chr = as.character(chr),
                                 pos = as.integer(pos))
    traits <- list(...)
    for (nm in names(traits)) dt[[nm]] <- traits[[nm]]
    dt[]
}

sig_table <- function(chr, pos, trait, SNPID = paste0(chr, ":", pos),
                      pvalue = 1e-8) {
    data.table::data.table(
        SNPID = SNPID, chr = as.character(chr), pos = as.integer(pos),
        pvalue = pvalue, pval_threshold = 1e-6, trait = trait
    )
}

# --- trait discovery -------------------------------------------------------

test_that("traits are every column except SNPID/chr/pos", {
    dt <- pvals(bio_1 = c(1e-9, 0.5, 0.5), bio_2 = c(0.5, 1e-9, 0.5))
    res <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_setequal(unique(res$sig_snps$trait), c("bio_1", "bio_2"))
})

test_that("exclude_cols keeps WZA metadata columns out of the trait list", {
    dt <- pvals(bio_1 = c(1e-9, 0.5, 0.5),
                n_snps = c(10, 20, 30), mean_maf = c(0.1, 0.2, 0.3))
    res <- quiet(find_significant_snps_per_trait(
        dt, "custom", 1e-6, 1, exclude_cols = c("n_snps", "mean_maf")))
    expect_identical(unique(res$sig_snps$trait), "bio_1")
    expect_false(any(c("n_snps", "mean_maf") %in% res$diagnostics$trait))
})

# --- the inclusive boundary ------------------------------------------------

test_that("a SNP sitting exactly ON the threshold is selected, not dropped", {
    # This is the whole point of `<=` at sig_snps.R:67. A strict `<` here would
    # silently drop the boundary SNP and every SNP tied with it.
    dt  <- pvals(bio_1 = c(1e-6, 1e-9, 0.5))
    res <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_identical(nrow(res$sig_snps), 2L)
    expect_true("1:100" %in% res$sig_snps$SNPID)
})

test_that("top N with `p <= threshold` returns all N, including the Nth", {
    dt  <- pvals(pos = 1:5, bio_1 = c(0.5, 0.01, 0.3, 0.02, 0.4))
    res <- quiet(find_significant_snps_per_trait(dt, "top", 3, 1))
    expect_identical(nrow(res$sig_snps), 3L)
    expect_setequal(res$sig_snps$pvalue, c(0.01, 0.02, 0.3))
})

test_that("NA p-values never select a SNP", {
    dt  <- pvals(bio_1 = c(1e-9, NA, 0.5))
    res <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_identical(nrow(res$sig_snps), 1L)
    expect_identical(res$sig_snps$SNPID, "1:100")
})

# --- non-ok statuses -------------------------------------------------------

test_that("a trait whose threshold cannot be computed contributes zero rows", {
    # qval needs >=10 tests; with 3 it returns too_few_tests and the trait is skipped.
    dt  <- pvals(bio_1 = c(1e-9, 1e-9, 1e-9))
    res <- quiet(find_significant_snps_per_trait(dt, "qval", 0.1, 1))
    expect_identical(nrow(res$sig_snps), 0L)
})

test_that("when nothing is significant the empty schema is returned, not an error", {
    dt  <- pvals(bio_1 = c(0.5, 0.6, 0.7))
    res <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_identical(nrow(res$sig_snps), 0L)
    expect_identical(colnames(res$sig_snps),
                     c("SNPID", "chr", "pos", "pvalue", "pval_threshold", "trait"))
    # Documented schema quirk: the empty path declares pos as numeric, while the
    # populated path inherits integer from the input.
    expect_true(is.numeric(res$sig_snps$pos))
})

test_that("chr is coerced to character so numeric-looking chromosomes survive", {
    dt <- data.table::data.table(SNPID = "1:100", chr = 1L, pos = 100L,
                                 bio_1 = 1e-9)
    res <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_true(is.character(res$sig_snps$chr))
    expect_identical(res$sig_snps$chr, "1")
})

test_that("the caller's p-value table is not mutated", {
    dt     <- pvals(bio_1 = c(1e-9, 0.5, 0.5))
    before <- copy(dt)
    quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_identical(colnames(dt), colnames(before))
    expect_identical(dt$chr, before$chr)
})

# --- diagnostics -----------------------------------------------------------

test_that("diagnostics carries one row per trait at cpu = 1", {
    dt  <- pvals(bio_1 = c(1e-9, 0.5, 0.5), bio_2 = c(0.5, 1e-9, 0.5))
    res <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    expect_identical(nrow(res$diagnostics), 2L)
    expect_setequal(res$diagnostics$trait, c("bio_1", "bio_2"))
    expect_identical(colnames(res$diagnostics),
                     c("trait", "status", "threshold", "n_tested", "n_na_dropped"))
    expect_true(all(res$diagnostics$status == "ok"))
})

test_that("diagnostics records n_na_dropped and a non-ok status", {
    dt  <- pvals(bio_1 = c(1e-9, NA, 0.5))
    res <- quiet(find_significant_snps_per_trait(dt, "bonf", 0.05, 1))
    expect_identical(res$diagnostics$n_tested, 2L)
    expect_identical(res$diagnostics$n_na_dropped, 1L)
})

# --- annotate_cross_trait_overlaps -----------------------------------------

test_that("annotate_cross_trait_overlaps adds its three columns and echoes the distance", {
    dt <- rbind(sig_table("1", 100, "bio_1"), sig_table("1", 110, "bio_2"))
    r  <- annotate_cross_trait_overlaps(dt, 1000L)
    expect_true(all(c("overlap_traits", "overlap_snps", "overlap_distance")
                    %in% colnames(r)))
    expect_true(all(r$overlap_distance == 1000L))
})

test_that("SNPs of different traits further apart than the distance do not overlap", {
    dt <- rbind(sig_table("1", 100, "bio_1"), sig_table("1", 500000, "bio_2"))
    r  <- annotate_cross_trait_overlaps(dt, 1000L)
    expect_true(all(is.na(r$overlap_traits)))
    expect_true(all(is.na(r$overlap_snps)))
})

test_that("a SNP of the same trait nearby is not counted as a cross-trait overlap", {
    dt <- rbind(sig_table("1", 100, "bio_1"), sig_table("1", 110, "bio_1"))
    r  <- annotate_cross_trait_overlaps(dt, 1000L)
    expect_true(all(is.na(r$overlap_traits)))
})

test_that("cross-trait overlap does not reach across chromosomes", {
    dt <- rbind(sig_table("1", 100, "bio_1"), sig_table("2", 100, "bio_2"))
    r  <- annotate_cross_trait_overlaps(dt, 1000000L)
    expect_true(all(is.na(r$overlap_traits)))
})

test_that("annotate_cross_trait_overlaps does not mutate the caller's table", {
    dt <- rbind(sig_table("1", 100, "bio_1"), sig_table("1", 110, "bio_2"))
    annotate_cross_trait_overlaps(dt, 1000L)
    expect_false("overlap_traits" %in% colnames(dt))
})

test_that("a zero-row table gets the three columns added rather than erroring", {
    # Note: this branch mutates the caller's table by reference — there is no
    # copy() before it, unlike the populated path.
    dt <- sig_table("1", 100, "bio_1")[0]
    r  <- annotate_cross_trait_overlaps(dt, 1000L)
    expect_identical(nrow(r), 0L)
    expect_true(all(c("overlap_traits", "overlap_snps", "overlap_distance")
                    %in% colnames(r)))
})
