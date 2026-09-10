# scripts/R/lib/regions.R
#
# Single-linkage region clustering. The criterion is stated in the file header:
# two consecutive SNPs merge when their gap <= 2*distance (equivalent to the
# GRanges extend-by-dist + reduce this replaced), and boundaries are
# [max(1, min(pos) - dist), max(pos) + dist].

# --- fixtures --------------------------------------------------------------

# BUILDER, not a shared object. cluster_snps_to_regions() copies its input, but
# regions.R:47 does `sig_snps[, min_pvalue := NA_real_]` and data.table `:=`
# mutates by reference — a shared fixture would leak a min_pvalue column into
# later tests.
# min_pvalue defaults to a real value rather than being absent: when the column is
# missing, regions.R:47 fills it with NA and regions.R:94's min(na.rm = TRUE) then
# warns "no non-missing arguments to min" on every single cluster. That is known
# bug #5 (Inf instead of NA), pinned once in test-known-bugs.R — it should not
# spray 35 warnings across the rest of the suite.
snps <- function(chr, pos, SNPID = paste0(chr, ":", pos),
                 min_pvalue = 1e-6, ...) {
    dt <- data.table::data.table(SNPID = SNPID, chr = as.character(chr),
                                 pos = as.integer(pos))
    if (!is.null(min_pvalue)) dt[, min_pvalue := min_pvalue]
    extra <- list(...)
    for (nm in names(extra)) dt[[nm]] <- extra[[nm]]
    dt[]
}

# --- .get_dist -------------------------------------------------------------

test_that(".get_dist passes an unnamed scalar through", {
    expect_identical(.get_dist(50000L, "1"), 50000L)
    expect_identical(.get_dist(50000, "anything"), 50000L)
})

test_that(".get_dist looks up a named per-chromosome vector", {
    spec <- c("1" = 1000L, "2" = 9000L)
    expect_identical(.get_dist(spec, "1"), 1000L)
    expect_identical(.get_dist(spec, 2), 9000L)   # numeric chr coerced to character
})

test_that(".get_dist warns and uses the max when the chromosome is absent", {
    spec <- c("1" = 1000L, "2" = 9000L)
    d <- expect_message(.get_dist(spec, "7"), "not in per-chr distance map")
    expect_identical(d, 9000L)
})

# --- cluster_snps_to_regions: boundaries -----------------------------------

test_that("a single SNP yields one region spanning [pos - dist, pos + dist]", {
    r <- cluster_snps_to_regions(snps("1", 5000), 1000L)
    expect_identical(nrow(r), 1L)
    expect_identical(r$start, 4000L)
    expect_identical(r$end, 6000L)
    expect_identical(r$length, 2000L)
    expect_identical(r$snp_count, 1L)
})

test_that("the start is floored at 1 when pos < dist", {
    r <- cluster_snps_to_regions(snps("1", 500), 1000L)
    expect_identical(r$start, 1L)
    expect_identical(r$end, 1500L)
})

test_that("distance 0 keeps every SNP in its own zero-length region", {
    r <- cluster_snps_to_regions(snps("1", c(100, 200, 300)), 0L)
    expect_identical(nrow(r), 3L)
    expect_identical(r$start, c(100L, 200L, 300L))
    expect_identical(r$end, c(100L, 200L, 300L))
    expect_true(all(r$snp_count == 1L))
})

# --- cluster_snps_to_regions: the merge criterion --------------------------

test_that("two SNPs merge when the gap is exactly 2*dist, and split at 2*dist + 1", {
    d <- 1000L
    merged <- cluster_snps_to_regions(snps("1", c(10000, 10000 + 2 * d)), d)
    expect_identical(nrow(merged), 1L)
    expect_identical(merged$snp_count, 2L)

    split <- cluster_snps_to_regions(snps("1", c(10000, 10000 + 2 * d + 1)), d)
    expect_identical(nrow(split), 2L)
    expect_true(all(split$snp_count == 1L))
})

test_that("clustering is single-linkage: a chain of small gaps merges end to end", {
    d <- 1000L
    # Consecutive gaps are all 2*d, so the whole chain is one region even though
    # the outer two SNPs are 6000 bp apart.
    r <- cluster_snps_to_regions(snps("1", c(10000, 12000, 14000, 16000)), d)
    expect_identical(nrow(r), 1L)
    expect_identical(r$snp_count, 4L)
    expect_identical(r$start, 9000L)
    expect_identical(r$end, 17000L)
})

test_that("input order does not matter — SNPs are sorted before clustering", {
    d <- 1000L
    forward <- cluster_snps_to_regions(snps("1", c(10000, 12000, 20000)), d)
    shuffled <- cluster_snps_to_regions(snps("1", c(20000, 10000, 12000)), d)
    expect_identical(forward$region_id, shuffled$region_id)
    expect_identical(forward$snp_count, shuffled$snp_count)
})

# --- cluster_snps_to_regions: structural invariants ------------------------

test_that("regions never span chromosomes even at identical positions", {
    r <- cluster_snps_to_regions(snps(c("1", "2"), c(10000, 10000)), 1000000L)
    expect_identical(nrow(r), 2L)
    expect_identical(sort(r$chr), c("1", "2"))
})

test_that("output is sorted by chr then start, and regions do not overlap within a chr", {
    r <- cluster_snps_to_regions(
        snps(c("2", "1", "1", "2"), c(50000, 90000, 10000, 10000)), 1000L
    )
    expect_identical(r$chr, sort(r$chr))
    for (this_chr in unique(r$chr)) {
        sub <- r[chr == this_chr]
        expect_identical(sub$start, sort(sub$start))
        if (nrow(sub) > 1) expect_true(all(sub$start[-1] > head(sub$end, -1)))
    }
})

test_that("region count never exceeds SNP count and start <= end always", {
    r <- cluster_snps_to_regions(snps("1", c(100, 5000, 900000)), 1000L)
    expect_lte(nrow(r), 3L)
    expect_true(all(r$start <= r$end))
    expect_identical(sum(r$snp_count), 3L)
})

# --- cluster_snps_to_regions: identity and content -------------------------

test_that("region_id is {chr}_{start}-{end}, with the trait appended when labelled", {
    plain <- cluster_snps_to_regions(snps("3", 5000), 1000L)
    expect_identical(plain$region_id, "3_4000-6000")
    expect_identical(plain$trait, "")

    labelled <- cluster_snps_to_regions(snps("3", 5000), 1000L, trait_label = "bio_1")
    expect_identical(labelled$region_id, "3_4000-6000_bio_1")
    expect_identical(labelled$trait, "bio_1")
})

test_that("snp_ids lists the member SNPs and methods lists only non-empty method columns", {
    dt <- snps("1", c(10000, 10500), EMMAX = c("bio_1", ""), LFMM = c("", ""))
    r <- cluster_snps_to_regions(dt, 1000L)
    expect_identical(r$snp_ids, "1:10000,1:10500")
    expect_identical(r$methods, "EMMAX")     # LFMM is empty throughout, so inactive
})

test_that("min_pvalue is the smallest p in the cluster", {
    dt <- snps("1", c(10000, 10500), min_pvalue = c(1e-3, 1e-8))
    expect_equal(cluster_snps_to_regions(dt, 1000L)$min_pvalue, 1e-8)
})

# --- schemas ---------------------------------------------------------------

test_that("the empty and populated schemas differ — 12 columns vs 10", {
    # Not a defect to fix here, but a contract downstream code depends on: the
    # empty constructor carries other_traits/other_snp_count that the populated
    # path does not, so rbind()ing the two produces NA-filled columns.
    empty <- cluster_snps_to_regions(NULL, 1000L)
    expect_identical(nrow(empty), 0L)
    expect_identical(colnames(empty),
                     c("region_id", "trait", "chr", "start", "end", "length",
                       "snp_count", "snp_ids", "methods", "min_pvalue",
                       "other_traits", "other_snp_count"))

    populated <- cluster_snps_to_regions(snps("1", 5000), 1000L)
    expect_identical(colnames(populated),
                     c("region_id", "trait", "chr", "start", "end", "length",
                       "snp_count", "snp_ids", "methods", "min_pvalue"))
})

test_that("a zero-row input returns the empty table, not an error", {
    expect_identical(nrow(cluster_snps_to_regions(snps("1", integer(0)), 1000L)), 0L)
})

# --- build_per_trait_regions -----------------------------------------------

test_that("build_per_trait_regions explodes comma-separated traits into one set per trait", {
    dt <- snps("1", c(10000, 10500, 900000),
               EMMAX = c("bio_1", "bio_1,bio_2", "bio_2"))
    r <- suppressMessages(build_per_trait_regions(dt, 1000L))
    expect_setequal(unique(r$trait), c("bio_1", "bio_2"))
    # bio_1 hits the two adjacent SNPs -> one region; bio_2 hits 10500 and 900000
    # -> two regions.
    expect_identical(nrow(r[trait == "bio_1"]), 1L)
    expect_identical(nrow(r[trait == "bio_2"]), 2L)
})

test_that("build_per_trait_regions reports cross-trait evidence for the same interval", {
    dt <- snps("1", c(10000, 10500), EMMAX = c("bio_1", "bio_2"))
    r <- suppressMessages(build_per_trait_regions(dt, 1000L))
    expect_identical(nrow(r), 2L)
    expect_identical(r[trait == "bio_1"]$other_traits, "bio_2")
    expect_identical(r[trait == "bio_1"]$other_snp_count, 1L)
    expect_identical(r[trait == "bio_2"]$other_traits, "bio_1")
})

test_that("build_per_trait_regions reports no cross-trait evidence when traits are far apart", {
    dt <- snps("1", c(10000, 900000), EMMAX = c("bio_1", "bio_2"))
    r <- suppressMessages(build_per_trait_regions(dt, 1000L))
    expect_true(all(r$other_traits == ""))
    expect_true(all(r$other_snp_count == 0L))
})

test_that("build_per_trait_regions returns the 12-column empty table on empty input", {
    r <- build_per_trait_regions(NULL, 1000L)
    expect_identical(nrow(r), 0L)
    expect_true(all(c("other_traits", "other_snp_count") %in% colnames(r)))
})

# --- build_combined_regions ------------------------------------------------

test_that("build_combined_regions clusters all traits together and drops the trait column", {
    dt <- snps("1", c(10000, 10500), EMMAX = c("bio_1", "bio_2"))
    r <- suppressMessages(build_combined_regions(dt, 1000L))
    expect_identical(nrow(r), 1L)                       # both traits, one interval
    expect_false("trait" %in% colnames(r))
    expect_true("traits" %in% colnames(r))
    expect_setequal(strsplit(r$traits, ",")[[1]], c("bio_1", "bio_2"))
})

test_that("build_combined_regions emits data-dependent <trait>_snps and <method>_snps columns", {
    dt <- snps("1", c(10000, 10500),
               EMMAX = c("bio_1", "bio_2"), LFMM = c("bio_1", ""))
    r <- suppressMessages(build_combined_regions(dt, 1000L))
    # Names are derived from the data, never hardcoded by the function.
    expect_true(all(c("bio_1_snps", "bio_2_snps", "EMMAX_snps", "LFMM_snps")
                    %in% colnames(r)))
    expect_identical(r$bio_1_snps, 1L)
    expect_identical(r$bio_2_snps, 1L)
    expect_identical(r$EMMAX_snps, 2L)                  # both SNPs carry an EMMAX value
    expect_identical(r$LFMM_snps, 1L)
})

test_that("build_combined_regions returns its own 10-column empty table", {
    r <- build_combined_regions(NULL, 1000L)
    expect_identical(nrow(r), 0L)
    expect_identical(colnames(r),
                     c("region_id", "chr", "start", "end", "length", "snp_count",
                       "snp_ids", "traits", "methods", "min_pvalue"))
})
