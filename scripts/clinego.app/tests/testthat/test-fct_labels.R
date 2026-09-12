# R/fct_labels.R — the display-string builders. Every failure here is cosmetic
# until it is not: build_region_labels() produces the NAMES of a selectInput's
# choices, so a collision silently removes a region from the picker, and
# format_si() is read as a count in a value box.
#
# The sharp one is build_region_labels()'s return TYPE, which differs by input:
# character(0) for no regions (:10) but a named LIST otherwise (:42). Anything
# treating the empty case as "a list with no elements" is fine; anything doing
# names() arithmetic on it is not.

mk_regions <- function(ids, snp_count = 1L) {
    data.table::data.table(region_id = ids,
                           snp_count = rep_len(snp_count, length(ids)))
}

# ---------------------------------------------------------------- build_region_labels

test_that("build_region_labels returns character(0) for no regions, not an empty list", {
    expect_identical(build_region_labels(NULL), character(0))
    expect_identical(build_region_labels(mk_regions(character(0))), character(0))
})

test_that("build_region_labels maps display label -> region_id, not the reverse", {
    got <- build_region_labels(mk_regions("1_100-200", snp_count = 3L))
    # selectInput wants choices = list(<label> = <value>), so the region id must
    # be the VALUE and the human-readable string the NAME.
    expect_type(got, "list")
    expect_identical(unname(unlist(got)), "1_100-200")
    expect_identical(names(got), "1_100-200 (3 SNPs)")
})

test_that("build_region_labels omits the E/P segment only when both counts are zero", {
    r <- mk_regions("1_100-200", snp_count = 5L)
    genes_zero <- data.table::data.table(region_id = "1_100-200",
                                         exon_snp_count = 0L, promoter_snp_count = 0L)
    expect_identical(names(build_region_labels(r, genes = genes_zero)),
                     "1_100-200 (5 SNPs)")

    # 0E is reachable on its own — the guard is (xe > 0 || yp > 0), not (xe > 0).
    genes_p <- data.table::data.table(region_id = "1_100-200",
                                      exon_snp_count = 0L, promoter_snp_count = 4L)
    expect_identical(names(build_region_labels(r, genes = genes_p)),
                     "1_100-200 (5 SNPs, 0E 4P)")
})

test_that("build_region_labels sums exon/promoter counts across a region's genes", {
    r <- mk_regions("1_100-200", snp_count = 9L)
    genes <- data.table::data.table(
        region_id = c("1_100-200", "1_100-200", "2_1-2"),
        exon_snp_count = c(1L, 2L, 99L),
        promoter_snp_count = c(3L, 4L, 99L))
    # The other region's row must not leak in.
    expect_identical(names(build_region_labels(r, genes = genes)),
                     "1_100-200 (9 SNPs, 3E 7P)")
})

test_that("build_region_labels ignores a genes table with no region_id column", {
    r <- mk_regions("1_100-200", snp_count = 2L)
    expect_identical(names(build_region_labels(r, genes = data.table::data.table(x = 1L))),
                     "1_100-200 (2 SNPs)")
})

test_that("build_region_labels counts GO rows per region and omits a zero", {
    r <- mk_regions(c("1_100-200", "2_300-400"), snp_count = 1L)
    enrich <- data.table::data.table(region_id = c("1_100-200", "1_100-200"))
    expect_identical(names(build_region_labels(r, enrich = enrich)),
                     c("1_100-200 (1 SNPs, 2 GO)", "2_300-400 (1 SNPs)"))
})

test_that("build_region_labels appends [HAP] only for listed regions", {
    r <- mk_regions(c("1_100-200", "2_300-400"), snp_count = 1L)
    got <- build_region_labels(r, hap_regions = "2_300-400")
    expect_identical(names(got), c("1_100-200 (1 SNPs)", "2_300-400 (1 SNPs [HAP])"))
})

test_that("build_region_labels assembles every segment in a fixed order", {
    r <- mk_regions("1_100-200", snp_count = 7L)
    genes  <- data.table::data.table(region_id = "1_100-200",
                                     exon_snp_count = 2L, promoter_snp_count = 1L)
    enrich <- data.table::data.table(region_id = rep("1_100-200", 3L))
    expect_identical(
        names(build_region_labels(r, genes = genes, enrich = enrich,
                                  hap_regions = "1_100-200")),
        "1_100-200 (7 SNPs, 2E 1P, 3 GO [HAP])")
})

test_that("build_region_labels can emit duplicate labels for distinct regions", {
    # Two different region_ids with identical summaries collide. vapply does not
    # uniquify, and a named list with repeated names silently loses entries in a
    # picker. Characterisation, not endorsement — it is the reason a caller must
    # not rely on names() being a key.
    r <- mk_regions(c("1_100-200", "1_100-200"), snp_count = 1L)
    got <- build_region_labels(r)
    expect_length(got, 2L)
    expect_identical(length(unique(names(got))), 1L)
})

# ---------------------------------------------------------------- format_hap_tag

test_that("format_hap_tag title-cases and slash-joins the underscore parts", {
    expect_identical(format_hap_tag("site_association"), "Site / Association")
    expect_identical(format_hap_tag("cluster_phenotype"), "Cluster / Phenotype")
})

test_that("format_hap_tag passes a single token through with no separator", {
    expect_identical(format_hap_tag("site"), "Site")
})

test_that("format_hap_tag is NOT vectorized", {
    # strsplit(tag, "_")[[1]] takes only the FIRST element, so a length-2 input
    # silently formats just the first tag rather than erroring.
    expect_identical(format_hap_tag(c("site_association", "cluster_x")),
                     "Site / Association")
})

test_that("format_hap_tag returns an empty string for an empty tag", {
    expect_identical(format_hap_tag(""), "")
})

# ---------------------------------------------------------------- format_si

test_that("format_si switches units exactly at the thresholds", {
    expect_identical(format_si(999),     "999")
    expect_identical(format_si(1000),    "1k")
    expect_identical(format_si(999999),  "1000k")   # not "1.0M": the M branch is >= 1e6
    expect_identical(format_si(1e6),     "1.0M")
    expect_identical(format_si(2.5e6),   "2.5M")
})

test_that("format_si renders NA as an em dash", {
    expect_identical(format_si(NA), "—")
    expect_identical(format_si(NA_integer_), "—")
})

test_that("format_si sends negatives down the as.character branch", {
    # Only >= comparisons are made, so a negative never reaches the k/M branches
    # and inherits R's default formatting, scientific notation included.
    expect_identical(format_si(-5), "-5")
    expect_identical(format_si(-5e6), "-5e+06")
})

test_that("format_si is NOT vectorized", {
    # The if() chain needs a length-1 condition.
    expect_error(format_si(c(1, 2)))
})
