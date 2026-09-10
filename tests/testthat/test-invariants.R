# Unit tests for scripts/R/lib/invariants.R.
#
# THE RULE THIS FILE EXISTS TO ENFORCE: every checker gets a PAIR of tests — a
# clean fixture that must return zero violations, and a deliberately broken one
# that must return exactly the expected violation. A checker tested only on
# clean data is indistinguishable from a checker that returns nothing at all.
# (ADAPTOGENE dossier, 2026-08-20: "A regression diff is worthless until a
# baseline control proves the test can return 'identical' for the right reason.")
#
# These tests are NOT in test-known-bugs.R and must never be skipped. They
# assert that a checker FIRES on broken input — that is new code working, not a
# quarantined defect. Only assertions about the pipeline's own correct
# behaviour belong behind skip().

dt <- data.table::data.table

regions_ok <- function() dt(
    region_id = c("1_100-200", "1_500-600", "2_100-200"),
    chr       = c("1", "1", "2"),
    start     = c(100L, 500L, 100L),
    end       = c(200L, 600L, 200L),
    length    = c(100L, 100L, 100L),
    snp_count = c(2L, 1L, 1L),
    snp_ids   = c("1:120,1:180", "1:550", "2:150")
)

snps_ok <- function() dt(
    SNPID = c("1:120", "1:180", "1:550", "2:150"),
    chr   = c("1", "1", "1", "2"),
    pos   = c(120L, 180L, 550L, 150L)
)

expect_no_violations <- function(v) {
    testthat::expect_s3_class(v, "data.table")
    testthat::expect_equal(names(v), names(no_violations()))
    testthat::expect_equal(nrow(v), 0L)
}

# ── violation plumbing ────────────────────────────────────────────────────────

test_that("no_violations is the schema contract every checker returns", {
    expect_no_violations(no_violations())
})

test_that("violation recycles key and detail against each other", {
    v <- violation("c", "error", "t", c("k1", "k2"), "same detail")
    expect_equal(nrow(v), 2L)
    expect_equal(v$detail, c("same detail", "same detail"))
    expect_equal(nrow(violation("c", "error", "t", character(), character())), 0L)
})

test_that("combine_violations drops empties and NULLs", {
    a <- violation("a", "error", "t", "k", "d")
    expect_equal(nrow(combine_violations(a, no_violations(), NULL, a)), 2L)
    expect_no_violations(combine_violations(NULL, no_violations()))
})

# ── check_regions_table ───────────────────────────────────────────────────────

test_that("check_regions_table passes a well-formed table", {
    expect_no_violations(check_regions_table(regions_ok(), snps_ok()))
    expect_no_violations(check_regions_table(NULL))
})

test_that("check_regions_table catches start > end", {
    r <- regions_ok(); r[1L, `:=`(start = 900L, end = 200L)]
    v <- check_regions_table(r)
    expect_true("region_start_after_end" %in% v$check)
    expect_equal(v[check == "region_start_after_end", key], "1_100-200")
})

test_that("check_regions_table catches length != end - start", {
    r <- regions_ok(); r[1L, length := 101L]   # the off-by-one an end-start+1 rule would give
    v <- check_regions_table(r)
    expect_true("region_length_mismatch" %in% v$check)
})

test_that("check_regions_table catches snp_count disagreeing with snp_ids", {
    r <- regions_ok(); r[1L, snp_count := 5L]
    v <- check_regions_table(r)
    expect_true("region_snp_count_disagrees_with_snp_ids" %in% v$check)
})

test_that("check_regions_table catches an empty region", {
    r <- regions_ok(); r[2L, `:=`(snp_count = 0L, snp_ids = "")]
    v <- check_regions_table(r)
    expect_true("region_has_no_snps" %in% v$check)
})

test_that("check_regions_table catches a duplicated region_id", {
    r <- rbind(regions_ok(), regions_ok()[1L])
    v <- check_regions_table(r)
    expect_true("region_id_duplicated" %in% v$check)
})

test_that("check_regions_table catches overlapping regions on the same chromosome", {
    r <- regions_ok(); r[2L, start := 150L]   # now overlaps 1_100-200
    v <- check_regions_table(r)
    expect_true("regions_overlap_within_group" %in% v$check)
})

test_that("check_regions_table allows the same span under DIFFERENT traits", {
    r <- rbind(
        dt(region_id = "1_100-200_a", trait = "a", chr = "1", start = 100L, end = 200L,
           length = 100L, snp_count = 1L, snp_ids = "1:120"),
        dt(region_id = "1_100-200_b", trait = "b", chr = "1", start = 100L, end = 200L,
           length = 100L, snp_count = 1L, snp_ids = "1:120")
    )
    expect_no_violations(check_regions_table(r))
})

test_that("check_regions_table catches a region naming a SNP that does not exist", {
    r <- regions_ok(); r[1L, `:=`(snp_ids = "1:120,9:999", snp_count = 2L)]
    v <- check_regions_table(r, snps_ok())
    expect_true("region_names_unknown_snp" %in% v$check)
    expect_match(v[check == "region_names_unknown_snp", detail], "9:999")
})

test_that("check_regions_table catches a SNP that lies outside its own region", {
    # 1:550 is a real SNP, but at position 550 — outside 1_100-200.
    r <- regions_ok(); r[1L, `:=`(snp_ids = "1:120,1:550", snp_count = 2L)]
    v <- check_regions_table(r, snps_ok())
    expect_true("region_snp_outside_bounds" %in% v$check)
})

test_that("check_regions_table reports a missing required column instead of erroring", {
    v <- check_regions_table(dt(region_id = "r1", chr = "1"))
    expect_equal(v$check, "regions_schema")
    expect_match(v$detail, "start")
})

# ── check_pvalues_table ───────────────────────────────────────────────────────

pv_ok <- function() dt(SNPID = c("1:1", "1:2", "1:3"), chr = c("1", "1", "1"),
                       pos = 1:3, EMMAX = c(0.1, 0.5, NA), LFMM = c(1e-9, 1, 0))

test_that("check_pvalues_table passes valid p-values, NA included", {
    expect_no_violations(check_pvalues_table(pv_ok(), c("EMMAX", "LFMM")))
})

test_that("check_pvalues_table catches out-of-range, negative, NaN and Inf", {
    for (bad in list(1.5, -0.1, NaN, Inf)) {
        p <- pv_ok(); p[1L, EMMAX := bad]
        v <- check_pvalues_table(p, c("EMMAX", "LFMM"))
        expect_true("pvalue_out_of_range" %in% v$check,
                    info = paste("did not fire on", bad))
    }
})

test_that("check_pvalues_table catches a duplicated SNPID", {
    p <- rbind(pv_ok(), pv_ok()[1L])
    v <- check_pvalues_table(p, c("EMMAX", "LFMM"))
    expect_true("snpid_duplicated" %in% v$check)
})

test_that("check_pvalues_table warns when chr was read as a number", {
    p <- pv_ok(); p[, chr := as.integer(chr)]
    v <- check_pvalues_table(p, c("EMMAX", "LFMM"))
    expect_true("chr_not_character" %in% v$check)
    expect_equal(v[check == "chr_not_character", severity], "warn")
})

# ── check_min_pvalue_against_sig_snps ─────────────────────────────────────────
#
# selected_snps.tsv's per-method columns hold TRAIT NAMES, not p-values, so the
# only place min_pvalue can be validated against is the long-format
# *_sig_snps_*.tsv tables. These fixtures use the real schema.

sel_ok  <- function() dt(SNPID = c("1:1", "1:2"), chr = "1", pos = c(1L, 2L),
                         EMMAX = c("bio_2", ""), LFMM = c("bio_2,bio_3", "bio_3"),
                         min_pvalue = c(1e-8, 1e-5))
sigs_ok <- function() dt(SNPID  = c("1:1", "1:1", "1:2"),
                         trait  = c("bio_2", "bio_3", "bio_3"),
                         pvalue = c(1e-8, 1e-6, 1e-5))

test_that("check_min_pvalue_against_sig_snps passes when min_pvalue is the smallest sig p", {
    expect_no_violations(check_min_pvalue_against_sig_snps(sel_ok(), sigs_ok()))
})

test_that("check_min_pvalue_against_sig_snps catches a stale min_pvalue", {
    s <- sel_ok(); s[1L, min_pvalue := 1e-6]   # 1e-8 is present in the sig table
    v <- check_min_pvalue_against_sig_snps(s, sigs_ok())
    expect_true("min_pvalue_disagrees_with_sig_tables" %in% v$check)
    expect_equal(v[check == "min_pvalue_disagrees_with_sig_tables", key], "1:1")
})

test_that("check_min_pvalue_against_sig_snps catches a selected SNP with no sig row at all", {
    s <- rbind(sel_ok(), dt(SNPID = "9:9", chr = "9", pos = 9L,
                            EMMAX = "bio_2", LFMM = "", min_pvalue = 1e-9))
    v <- check_min_pvalue_against_sig_snps(s, sigs_ok())
    expect_true("selected_snp_absent_from_sig_tables" %in% v$check)
    expect_equal(v[check == "selected_snp_absent_from_sig_tables", key], "9:9")
})

test_that("check_min_pvalue_against_sig_snps tolerates pooling several adjust variants", {
    # The same (SNPID, trait, method) carries the same p at every threshold, and
    # a looser threshold only adds rows with LARGER p — so pooling cannot lower
    # the observed minimum.
    pooled <- rbind(sigs_ok(), dt(SNPID = "1:1", trait = "bio_2", pvalue = 1e-8),
                    dt(SNPID = "1:1", trait = "bio_9", pvalue = 3e-4))
    expect_no_violations(check_min_pvalue_against_sig_snps(sel_ok(), pooled))
})

# ── selected_snps_traits ──────────────────────────────────────────────────────

test_that("selected_snps_traits reads trait names from the VALUES, not the column names", {
    # The trap: the column names are METHODS (EMMAX, LFMM); the traits are inside.
    expect_setequal(selected_snps_traits(sel_ok()), c("bio_2", "bio_3"))
})

test_that("selected_snps_traits handles the quoted-empty marker and blank cells", {
    s <- dt(SNPID = "1:1", chr = "1", pos = 1L,
            EMMAX = '""', LFMM = "", RDA = NA_character_, min_pvalue = 1e-8)
    expect_equal(selected_snps_traits(s), character(0))
})

test_that("selected_snps_traits returns character(0) rather than erroring on empty input", {
    expect_equal(selected_snps_traits(NULL), character(0))
    expect_equal(selected_snps_traits(dt(SNPID = character(), min_pvalue = numeric())),
                 character(0))
})

# ── check_chromosome_names ────────────────────────────────────────────────────

test_that("check_chromosome_names accepts a SUBSET, not only an exact match", {
    # The genes table covers fewer chromosomes than the SNP table. Legitimate:
    # a region can contain no gene. Equality here would be a false positive.
    v <- check_chromosome_names(
        list(snps  = dt(chr = c("1", "2", "3")),
             genes = dt(chr = c("1", "3"))),
        canonical = c("1", "2", "3"))
    expect_no_violations(v)
})

test_that("check_chromosome_names catches a surviving 'chr' prefix", {
    v <- check_chromosome_names(list(genes = dt(chr = c("1", "chr2H"))))
    expect_true("chromosome_name_not_normalized" %in% v$check)
    expect_equal(v[check == "chromosome_name_not_normalized", key], "chr2H")
})

test_that("check_chromosome_names catches a chromosome absent from the canonical set", {
    v <- check_chromosome_names(list(genes = dt(chr = c("1", "7"))),
                                canonical = c("1", "2"))
    expect_true("chromosome_not_in_canonical_set" %in% v$check)
    expect_equal(v[check == "chromosome_not_in_canonical_set", key], "7")
})

# ── check_summary_accounting ──────────────────────────────────────────────────

summary_ok <- function() dt(
    step = "processing",
    metric = c("samples_total", "samples_after_filtering", "samples_removed",
               "samples_het_outliers_removed", "samples_removed_relatedness",
               "samples_with_coordinates", "samples_dropped_missing_coordinates"),
    value = c("51", "47", "1", "0", "3", "46", "1")
)

test_that("check_summary_accounting passes a summary that closes", {
    # 51 - 1 - 0 - 3 = 47, and 46 + 1 = 47. These are SIMDATA's real numbers.
    expect_no_violations(check_summary_accounting(summary_ok()))
})

test_that("check_summary_accounting catches sample arithmetic that does not close", {
    s <- summary_ok(); s[metric == "samples_after_filtering", value := "49"]
    v <- check_summary_accounting(s)
    expect_true("sample_accounting_does_not_close" %in% v$check)
})

test_that("check_summary_accounting catches coordinate arithmetic that does not close", {
    s <- summary_ok(); s[metric == "samples_with_coordinates", value := "40"]
    v <- check_summary_accounting(s)
    expect_true("coordinate_accounting_does_not_close" %in% v$check)
})

test_that("check_summary_accounting catches a duplicated (step, metric) pair", {
    s <- rbind(summary_ok(), summary_ok()[1L])
    v <- check_summary_accounting(s)
    expect_true("summary_metric_duplicated" %in% v$check)
})

test_that("check_summary_accounting catches the predictor-count discrepancy", {
    # SIMDATA's real values: 19 predictors listed, n_climate_variables = 16,
    # because write_summary.R:362 does ncol(climate_site) - 4 while that file
    # carries only `sample` plus the bio columns. Filed 2026-09-10.
    s <- rbind(summary_ok(), dt(
        step = "structure",
        metric = c("climate_predictors", "n_climate_variables"),
        value = c(paste0("bio_", 1:19, collapse = ","), "16")))
    v <- check_summary_accounting(s)
    expect_true("climate_predictor_count_disagrees" %in% v$check)

    s_fixed <- data.table::copy(s)
    s_fixed[metric == "n_climate_variables", value := "19"]
    expect_false("climate_predictor_count_disagrees" %in% check_summary_accounting(s_fixed)$check)
})

test_that("check_summary_counts compares a stated count with the real table", {
    s <- dt(step = "gea", metric = "selected_snps_total", value = "8")
    expect_no_violations(check_summary_counts(s, "gea", list(selected_snps_total = 8L)))
    v <- check_summary_counts(s, "gea", list(selected_snps_total = 12L))
    expect_equal(v$check, "summary_count_disagrees_with_table")
})

# ── check_genes_table ─────────────────────────────────────────────────────────

genes_ok <- function() dt(region_id = c("1_100-200", "1_100-200"),
                          gene_id = c("g1", "g2"),
                          gene_start = c(110L, 150L), gene_end = c(120L, 160L),
                          exon_snp_count = c(1L, 0L))

test_that("check_genes_table passes a well-formed table", {
    expect_no_violations(check_genes_table(genes_ok(), regions_ok()))
})

test_that("check_genes_table catches gene_start > gene_end", {
    g <- genes_ok(); g[1L, `:=`(gene_start = 900L, gene_end = 120L)]
    v <- check_genes_table(g, regions_ok())
    expect_true("gene_start_after_end" %in% v$check)
})

test_that("check_genes_table catches a gene pointing at a region that does not exist", {
    g <- genes_ok(); g[1L, region_id := "9_1-2"]
    v <- check_genes_table(g, regions_ok())
    expect_true("gene_references_unknown_region" %in% v$check)
})

test_that("check_genes_table catches exon_snp_count exceeding the region's SNP count", {
    # This is the shape genes_in_regions.R:174 produces on real data: it counts
    # features hit, which can exceed the SNPs actually present.
    g <- genes_ok(); g[1L, exon_snp_count := 9L]   # region 1_100-200 holds 2 SNPs
    v <- check_genes_table(g, regions_ok())
    expect_true("gene_snp_count_exceeds_region" %in% v$check)
})

# ── check_offsets_table ───────────────────────────────────────────────────────

test_that("check_offsets_table passes finite non-negative offsets", {
    expect_no_violations(check_offsets_table(
        dt(site = "s1", sample = c("a", "b"), genetic_offset = c(0, 0.13))))
})

test_that("check_offsets_table catches non-finite and negative offsets", {
    for (bad in list(NA_real_, NaN, Inf)) {
        v <- check_offsets_table(dt(sample = "a", genetic_offset = bad))
        expect_true("offset_not_finite" %in% v$check, info = paste("missed", bad))
    }
    v <- check_offsets_table(dt(sample = "a", genetic_offset = -0.5))
    expect_true("offset_negative" %in% v$check)
})

test_that("check_offsets_table reports a missing value column instead of erroring", {
    v <- check_offsets_table(dt(sample = "a", something_else = 1))
    expect_equal(v$check, "offset_column_missing")
})

# ── check_single_threshold_variant ────────────────────────────────────────────

test_that("check_single_threshold_variant passes one threshold per method", {
    expect_no_violations(check_single_threshold_variant(
        list(EMMAX = "bonf_0.05", LFMM = "bonf_0.05")))
    expect_no_violations(check_single_threshold_variant(list()))
})

test_that("check_single_threshold_variant catches a stale variant left by a config change", {
    # SIMDATA's real shape: RDA is configured at bonf 0.01, but a bonf_0.05 sig
    # table from an earlier config is still on disk — Snakemake does not delete
    # outputs a config change orphaned.
    v <- check_single_threshold_variant(
        list(EMMAX = "bonf_0.05", RDA = c("bonf_0.01", "bonf_0.05")))
    expect_equal(v$check, "multiple_threshold_variants_on_disk")
    expect_equal(v$key, "RDA")
    expect_match(v$detail, "bonf_0.01, bonf_0.05")
})

test_that("check_single_threshold_variant ignores a repeated identical token", {
    expect_no_violations(check_single_threshold_variant(
        list(EMMAX = c("bonf_0.05", "bonf_0.05"))))
})

# ── check_referential_integrity ───────────────────────────────────────────────

test_that("check_referential_integrity passes when downstream is a subset", {
    expect_no_violations(check_referential_integrity(c("bio_2", "bio_3"),
                                                     c("bio_1", "bio_2", "bio_3")))
})

test_that("check_referential_integrity catches a downstream value with no upstream source", {
    # The SIMDATA shape: pairwise_overlap_table names bio_1, but GEA was re-run
    # with a predictor set that no longer contains it.
    v <- check_referential_integrity(c("bio_1", "bio_2"), c("bio_2", "bio_3"),
                                     table_name = "pairwise_overlap_table.tsv",
                                     what = "trait")
    expect_equal(nrow(v), 1L)
    expect_equal(v$key, "bio_1")
    expect_match(v$detail, "not present upstream")
})

# ── check_column_names_unique ─────────────────────────────────────────────────

test_that("check_column_names_unique passes distinct names and catches repeats", {
    expect_no_violations(check_column_names_unique(c("a", "b", "c"), "t.tsv"))
    v <- check_column_names_unique(c("region_id", "exon_snps", "region_id", "exon_snps"), "t.tsv")
    expect_setequal(v$key, c("region_id", "exon_snps"))
    expect_equal(unique(v$check), "duplicate_column_names")
})

# ── check_filtering_monotone ──────────────────────────────────────────────────

test_that("check_filtering_monotone passes a shrinking filter chain", {
    # SIMDATA's real chain.
    expect_no_violations(check_filtering_monotone(dt(
        stage     = c("Raw VCF", "After sample missingness filter", "After LD pruning"),
        n_samples = c(51L, 50L, 47L),
        n_snps    = c(354L, 354L, 247L))))
})

test_that("check_filtering_monotone catches a stage that gains rows", {
    v <- check_filtering_monotone(dt(
        stage     = c("Raw VCF", "After MAF filter"),
        n_samples = c(51L, 51L),
        n_snps    = c(354L, 400L)))
    expect_true("filtering_stage_not_monotone" %in% v$check)
    expect_match(v$detail, "n_snps rises from 354 to 400")
})

# ── check_sig_snp_overlaps ────────────────────────────────────────────────────

test_that("check_sig_snp_overlaps passes a correct cross-trait annotation", {
    s <- dt(SNPID = "1:100", trait = "bio_2",
            overlap_traits = "bio_3", overlap_snps = "1:150")
    expect_no_violations(check_sig_snp_overlaps(s, known_snpids = c("1:100", "1:150")))
})

test_that("check_sig_snp_overlaps catches overlap_traits naming the SNP's own trait", {
    # SIMDATA's real shape: 1:10683918 (trait bio_2) reports overlap_traits =
    # bio_2, because sig_snps.R:153 compares the wrong sides of the foverlaps
    # result. Filed at docs/pipeline_improvement_requests.md.
    s <- dt(SNPID = "1:10683918", trait = "bio_2",
            overlap_traits = "bio_2", overlap_snps = "1:9683918")
    v <- check_sig_snp_overlaps(s)
    expect_true("overlap_traits_includes_own_trait" %in% v$check)
})

test_that("check_sig_snp_overlaps catches overlap_snps naming a window start rather than a SNP", {
    # 1:9683918 is exactly pos - 1000000: the window boundary, not a real SNP.
    s <- dt(SNPID = "1:10683918", trait = "bio_2",
            overlap_traits = "bio_3", overlap_snps = "1:9683918")
    v <- check_sig_snp_overlaps(s, known_snpids = c("1:10683918", "1:10900000"))
    expect_true("overlap_snps_names_unknown_snp" %in% v$check)
})

test_that("check_sig_snp_overlaps ignores rows with no overlap recorded", {
    s <- dt(SNPID = c("1:1", "1:2"), trait = "bio_2",
            overlap_traits = c("", NA), overlap_snps = c("", NA))
    expect_no_violations(check_sig_snp_overlaps(s, known_snpids = c("1:1", "1:2")))
})
