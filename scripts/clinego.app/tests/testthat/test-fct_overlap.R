# Unit tests for R/fct_overlap.R — the GEA x GWAS region overlap computation.
#
# This file is data.table-only: no shiny, no file I/O, no subprocess, and no
# dependency on the shared pipeline libraries that zzz.R sources from /pipeline.
# It therefore runs anywhere, with no container and no mount.
#
# Two contracts below are pinned deliberately because they are surprising, not
# because they are obviously right:
#   * .apply_bounds() silently falls back to "union" for an unrecognised
#     strategy (fct_overlap.R:119-122) — a typo is not an error today.
#   * assign_region_ids_from_regions() mutates its argument BY REFERENCE, on
#     every path including the early returns. Callers defend with
#     data.table::copy() (mod_gea.R:547-552, mod_gwas.R:502-507).
# If either is ever fixed, this file is where the change must be made explicit.

regs <- function(region_id, chr, start, end) {
    data.table::data.table(
        region_id = as.character(region_id),
        chr       = as.character(chr),
        start     = as.integer(start),
        end       = as.integer(end)
    )
}

sig <- function(SNPID, chr, pos, pvalue = 1e-8, method = "EMMAX", trait = "bio_1") {
    data.table::data.table(
        SNPID  = as.character(SNPID),
        chr    = as.character(chr),
        pos    = as.integer(pos),
        pvalue = as.numeric(pvalue),
        method = as.character(method),
        trait  = as.character(trait)
    )
}

# ── compute_region_overlaps() ─────────────────────────────────────────────────

test_that("compute_region_overlaps pairs only genuinely overlapping regions", {
    gea  <- regs(c("g1", "g2"), c("1", "2"), c(100L, 100L), c(200L, 200L))
    gwas <- regs(c("w1", "w2"), c("1", "2"), c(150L, 500L), c(250L, 600L))

    out <- compute_region_overlaps(gea, gwas)

    expect_equal(nrow(out), 1L)
    expect_equal(out$gea_region_id,  "g1")
    expect_equal(out$gwas_region_id, "w1")
    expect_equal(out$overlap_start, 150L)
    expect_equal(out$overlap_end,   200L)
    expect_equal(out$union_start,   100L)
    expect_equal(out$union_end,     250L)
})

test_that("compute_region_overlaps keeps the GEA and GWAS sides straight", {
    # foverlaps(x = gwas, y = gea) returns the y side UNPREFIXED and the x side
    # as i.*  (fct_overlap.R:52-61). This file reads them the right way round;
    # sig_snps.R:152 and genes_in_regions.R:174 read the wrong one. Pin it here
    # so the correct call site stays correct.
    gea  <- regs("g1", "1", 100L, 200L)
    gwas <- regs("w1", "1", 120L, 300L)

    out <- compute_region_overlaps(gea, gwas)

    expect_equal(out$gea_region_id,  "g1")
    expect_equal(out$gea_start,      100L)
    expect_equal(out$gea_end,        200L)
    expect_equal(out$gwas_region_id, "w1")
    expect_equal(out$gwas_start,     120L)
    expect_equal(out$gwas_end,       300L)
})

test_that("compute_region_overlaps finds the same pair set when the arguments are swapped", {
    a <- regs(c("a1", "a2"), c("1", "1"), c(100L, 900L), c(200L, 950L))
    b <- regs(c("b1", "b2"), c("1", "2"), c(150L, 100L), c(250L, 200L))

    fwd <- compute_region_overlaps(a, b)
    rev <- compute_region_overlaps(b, a)

    expect_equal(nrow(fwd), nrow(rev))
    # Same unordered pairs, with the roles exchanged.
    expect_setequal(paste(fwd$gea_region_id, fwd$gwas_region_id),
                    paste(rev$gwas_region_id, rev$gea_region_id))
    # The symmetric quantities must agree exactly.
    expect_equal(sort(fwd$overlap_start), sort(rev$overlap_start))
    expect_equal(sort(fwd$union_end),     sort(rev$union_end))
})

test_that("compute_region_overlaps counts a shared boundary base exactly once", {
    gea  <- regs("g1", "1", 100L, 200L)
    gwas <- regs("w1", "1", 200L, 300L)

    out <- compute_region_overlaps(gea, gwas)

    expect_equal(nrow(out), 1L)
    expect_equal(out$overlap_start, 200L)
    expect_equal(out$overlap_end,   200L)
})

test_that("compute_region_overlaps never pairs across chromosomes", {
    gea  <- regs("g1", "1", 100L, 200L)
    gwas <- regs("w1", "2", 100L, 200L)

    expect_equal(nrow(compute_region_overlaps(gea, gwas)), 0L)
})

test_that("compute_region_overlaps returns the typed empty skeleton, not NULL or NA", {
    gea  <- regs("g1", "1", 100L, 200L)
    gwas <- regs("w1", "1", 900L, 950L)

    disjoint <- compute_region_overlaps(gea, gwas)
    empty_in <- compute_region_overlaps(gea, regs(character(), character(), integer(), integer()))
    null_in  <- compute_region_overlaps(NULL, gwas)

    for (out in list(disjoint, empty_in, null_in)) {
        expect_s3_class(out, "data.table")
        expect_equal(nrow(out), 0L)
        expect_equal(names(out), names(.empty_overlap_pairs()))
    }
})

test_that("compute_region_overlaps output is sorted by chr then gea_start", {
    gea  <- regs(c("g1", "g2", "g3"), c("2", "1", "1"), c(100L, 900L, 100L), c(200L, 950L, 200L))
    gwas <- regs(c("w1", "w2", "w3"), c("2", "1", "1"), c(100L, 900L, 100L), c(200L, 950L, 200L))

    out <- compute_region_overlaps(gea, gwas)

    expect_equal(out$chr,       c("1", "1", "2"))
    expect_equal(out$gea_start, c(100L, 900L, 100L))
})

# ── overlap_pct arithmetic ────────────────────────────────────────────────────

test_that("overlap_pct is 100 for identical regions", {
    out <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                   regs("w1", "1", 100L, 200L))
    expect_equal(out$overlap_pct, 100)
})

test_that("overlap_pct stays finite for a zero-length region", {
    # min_len is 0 here; the pmax(min_len, 1L) guard at fct_overlap.R:71 is what
    # keeps this from being 0/0 = NaN.
    out <- compute_region_overlaps(regs("g1", "1", 100L, 100L),
                                   regs("w1", "1", 100L, 100L))
    expect_equal(nrow(out), 1L)
    expect_true(is.finite(out$overlap_pct))
})

test_that("overlap_pct is relative to the SMALLER region", {
    # GEA 100-200 (len 100), GWAS 150-1150 (len 1000), overlap 150-200 = 50.
    # Relative to the smaller region that is 50%, not 5%.
    out <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                   regs("w1", "1", 150L, 1150L))
    expect_equal(out$overlap_pct, 50)
})

# ── .apply_bounds() ───────────────────────────────────────────────────────────

test_that(".apply_bounds honours each named strategy", {
    pair <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                    regs("w1", "1", 150L, 300L))

    expect_equal(.apply_bounds(pair, "union"),        list(chr = "1", start = 100L, end = 300L))
    expect_equal(.apply_bounds(pair, "intersection"), list(chr = "1", start = 150L, end = 200L))
    expect_equal(.apply_bounds(pair, "gea_only"),     list(chr = "1", start = 100L, end = 200L))
    expect_equal(.apply_bounds(pair, "gwas_only"),    list(chr = "1", start = 150L, end = 300L))
})

test_that(".apply_bounds silently falls back to union for an unknown strategy", {
    # CURRENT BEHAVIOUR, not an endorsement: the unnamed final switch() branch
    # (fct_overlap.R:119-122) makes a typo'd strategy return plausible wrong
    # bounds with no error. Filed in docs/pipeline_improvement_requests.md.
    pair <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                    regs("w1", "1", 150L, 300L))

    expect_equal(.apply_bounds(pair, "intersecton"), .apply_bounds(pair, "union"))
    expect_equal(.apply_bounds(pair, ""),            .apply_bounds(pair, "union"))
})

# ── compute_all_overlap_regions() ─────────────────────────────────────────────

test_that("compute_all_overlap_regions mints a region_id that round-trips", {
    pairs <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                     regs("w1", "1", 150L, 300L))
    out <- compute_all_overlap_regions(pairs, NULL, NULL, strategy = "union")

    expect_equal(out$region_id, "1_100-300")
    parsed <- .parse_region_id(out$region_id)
    expect_equal(parsed$chr,   "1")
    expect_equal(parsed$start, 100)
    expect_equal(parsed$end,   300)
})

test_that("compute_all_overlap_regions length is end - start, matching regions.R:90", {
    pairs <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                     regs("w1", "1", 150L, 300L))
    out <- compute_all_overlap_regions(pairs, NULL, NULL, strategy = "union")
    expect_equal(out$length, out$end - out$start)
    expect_equal(out$length, 200L)
})

test_that("compute_all_overlap_regions counts unique SNPIDs across both sources", {
    pairs <- compute_region_overlaps(regs("g1", "1", 100L, 300L),
                                     regs("w1", "1", 100L, 300L))
    gea_snps  <- sig(c("1:150", "1:160"), "1", c(150L, 160L), trait = "bio_1", method = "LFMM")
    gwas_snps <- sig(c("1:160", "1:900"), "1", c(160L, 900L), trait = "height", method = "BLINK")

    out <- compute_all_overlap_regions(pairs, NULL, NULL, strategy = "union",
                                       gea_sig_snps = gea_snps, gwas_sig_snps = gwas_snps)

    # 1:900 is outside the bounds; 1:160 is present in both sources but counts once.
    expect_equal(out$snp_count, 2L)
    expect_setequal(strsplit(out$snp_ids, ",")[[1]], c("1:150", "1:160"))
    expect_setequal(strsplit(out$traits,  ",")[[1]], c("bio_1", "height"))
    expect_setequal(strsplit(out$methods, ",")[[1]], c("BLINK", "LFMM"))
})

test_that("compute_all_overlap_regions gives min_pvalue NA when no SNP falls in bounds", {
    pairs <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                     regs("w1", "1", 150L, 300L))
    out <- compute_all_overlap_regions(pairs, NULL, NULL, strategy = "union")

    expect_true(is.na(out$min_pvalue))
    expect_equal(out$snp_count, 0L)
    expect_equal(out$snp_ids, "")
})

test_that("compute_all_overlap_regions deduplicates identical bounds from different pairs", {
    # Two GWAS regions inside one GEA region: under gea_only both pairs collapse
    # to the same bounds, hence the same region_id.
    pairs <- compute_region_overlaps(regs("g1", "1", 100L, 500L),
                                     regs(c("w1", "w2"), c("1", "1"),
                                          c(150L, 300L), c(200L, 350L)))
    expect_equal(nrow(pairs), 2L)

    out <- compute_all_overlap_regions(pairs, NULL, NULL, strategy = "gea_only")
    expect_equal(nrow(out), 1L)
    expect_equal(out$region_id, "1_100-500")
})

test_that("compute_all_overlap_regions returns the typed empty skeleton on empty input", {
    for (out in list(compute_all_overlap_regions(.empty_overlap_pairs(), NULL, NULL),
                     compute_all_overlap_regions(NULL, NULL, NULL))) {
        expect_s3_class(out, "data.table")
        expect_equal(nrow(out), 0L)
        expect_equal(names(out), names(.empty_overlap_regions()))
    }
})

test_that("compute_all_overlap_regions ignores its gea_regions/gwas_regions arguments", {
    # Neither is read in the body (fct_overlap.R:142-205) — everything comes from
    # overlap_pairs and the two sig-SNP tables. Pinned so that a future change
    # which starts using them cannot land silently.
    pairs <- compute_region_overlaps(regs("g1", "1", 100L, 200L),
                                     regs("w1", "1", 150L, 300L))
    with_regions <- compute_all_overlap_regions(pairs,
                                                regs("zzz", "9", 1L, 2L),
                                                regs("qqq", "9", 3L, 4L))
    without      <- compute_all_overlap_regions(pairs, NULL, NULL)
    expect_equal(with_regions, without)
})

# ── assign_region_ids_from_regions() ──────────────────────────────────────────

test_that("assign_region_ids_from_regions stamps SNPs inside a region and NA outside", {
    snps <- sig(c("1:150", "1:900"), "1", c(150L, 900L))
    out  <- assign_region_ids_from_regions(data.table::copy(snps),
                                           regs("r1", "1", 100L, 200L))

    expect_equal(out$region_id, c("r1", NA_character_))
    expect_type(out$region_id, "character")
})

test_that("assign_region_ids_from_regions MUTATES its argument by reference", {
    # This is why mod_gea.R:547-552 and mod_gwas.R:502-507 wrap the call in
    # data.table::copy(). Asserted deliberately so the contract is documented
    # rather than folklore.
    snps <- sig("1:150", "1", 150L)
    expect_false("region_id" %in% names(snps))

    assign_region_ids_from_regions(snps, regs("r1", "1", 100L, 200L))

    expect_true("region_id" %in% names(snps))
    expect_equal(snps$region_id, "r1")
})

test_that("assign_region_ids_from_regions blanks region_id even on the early-return paths", {
    for (bad_regions in list(NULL,
                             regs(character(), character(), integer(), integer()),
                             data.table::data.table(chr = "1", start = 100L))) {
        snps <- sig("1:150", "1", 150L)
        snps[, region_id := "STALE"]

        out <- assign_region_ids_from_regions(snps, bad_regions)

        expect_true(all(is.na(out$region_id)))
        expect_true(all(is.na(snps$region_id)))   # same object, by reference
    }
})
