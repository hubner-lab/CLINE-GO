# Unit tests for the PURE subset of R/fct_regions.R.
#
# That file is ~1100 lines, but most of it is processx launchers, a system()/awk
# caller and disk readers. Only the functions below are pure enough to unit-test;
# the rest belongs to a later smoke-test tier.
#
# SCOPE WARNING — what a green run here does NOT tell you. zzz.R:39-42 shares
# only regions.R and pval_threshold.R with the pipeline. find_genes_in_region()
# is the APP'S OWN reimplementation of gene finding; scripts/R/lib/genes_in_regions.R
# is a separate implementation the app never calls. These tests pin the app's
# behaviour and say nothing about whether the two agree. Same for fct_combine.R
# vs combine_sigsnps.R, and for region_distance.R, which the app does not source
# at all.

# ── .parse_region_id() ────────────────────────────────────────────────────────

test_that(".parse_region_id round-trips with format_region_id", {
    parsed <- .parse_region_id("1_100000-200000")
    expect_equal(parsed, list(chr = "1", start = 100000L, end = 200000L))
    expect_equal(format_region_id("1_100000-200000"), "1:100000-200000")
})

test_that(".parse_region_id handles a non-numeric chromosome name", {
    expect_equal(.parse_region_id("2H_100-200"),
                 list(chr = "2H", start = 100L, end = 200L))
})

test_that(".parse_region_id returns NULL rather than erroring on malformed input", {
    expect_null(.parse_region_id("no-underscore-here"))
    expect_null(.parse_region_id("1_100-200-300"))
    expect_null(.parse_region_id("1_abc-def"))
    expect_null(.parse_region_id("1_100"))
})

test_that(".parse_region_id splits on the FIRST underscore, so a per-trait id does not parse", {
    # The pipeline writes per-trait ids as "{chr}_{start}-{end}_{trait}", e.g.
    # "1_9391481-11683918_bio_2". Splitting on the first underscore leaves
    # "11683918_bio_2" as the end coordinate, which is not an integer -> NULL.
    # Current behaviour; pinned because .parse_region_id has a second consumer
    # at fct_discovery.R:382-390 (haplotype fuzzy matching).
    expect_null(.parse_region_id("1_9391481-11683918_bio_2"))
    # A chromosome name containing an underscore is likewise unrepresentable:
    # "chr_A_100-200" splits to chr="chr" and "A_100-200", whose first half is
    # not an integer, so the whole parse returns NULL rather than chr="chr_A".
    expect_null(.parse_region_id("chr_A_100-200"))
})

# ── .regions_overlap() ────────────────────────────────────────────────────────

test_that(".regions_overlap uses closed intervals on the same chromosome", {
    r <- function(chr, s, e) list(chr = chr, start = s, end = e)

    expect_true(.regions_overlap(r("1", 100, 200), r("1", 150, 250)))
    expect_true(.regions_overlap(r("1", 100, 200), r("1", 200, 300)))  # touching counts
    expect_true(.regions_overlap(r("1", 100, 200), r("1", 120, 130)))  # contained
    expect_false(.regions_overlap(r("1", 100, 200), r("1", 201, 300)))
    expect_false(.regions_overlap(r("1", 100, 200), r("2", 100, 200))) # different chr
})

# ── find_genes_in_region() ────────────────────────────────────────────────────

gff <- function(gene_id, chr, start, end) {
    data.table::data.table(
        gene_id = as.character(gene_id),
        chr     = as.character(chr),
        start   = as.integer(start),
        end     = as.integer(end)
    )
}
region <- function(chr, start, end) {
    data.table::data.table(chr = as.character(chr),
                           start = as.integer(start),
                           end = as.integer(end))
}

test_that("find_genes_in_region keeps genes inside, on the edge, and straddling a boundary", {
    genes <- gff(c("inside", "left_edge", "right_edge", "straddle_left", "outside"),
                 "1",
                 c(1200L,  1000L, 2000L, 900L,  5000L),
                 c(1300L,  1050L, 2100L, 1010L, 5100L))

    out <- find_genes_in_region(genes, region("1", 1000L, 2000L), promoter_length = 0L)

    expect_setequal(out$gene_id, c("inside", "left_edge", "right_edge", "straddle_left"))
})

test_that("find_genes_in_region excludes genes on another chromosome", {
    genes <- gff(c("same", "other"), c("1", "2"), c(1200L, 1200L), c(1300L, 1300L))
    out <- find_genes_in_region(genes, region("1", 1000L, 2000L), promoter_length = 0L)
    expect_equal(out$gene_id, "same")
})

test_that("find_genes_in_region extends the window UPSTREAM only", {
    genes <- gff(c("upstream", "downstream"), "1",
                 c(500L,  2500L),
                 c(600L,  2600L))

    none <- find_genes_in_region(genes, region("1", 1000L, 2000L), promoter_length = 0L)
    expect_equal(nrow(none), 0L)

    ext <- find_genes_in_region(genes, region("1", 1000L, 2000L), promoter_length = 1000L)
    expect_equal(ext$gene_id, "upstream")   # downstream stays excluded
})

test_that("find_genes_in_region floors the extended window at 1", {
    genes <- gff("g1", "1", 1L, 50L)
    out <- find_genes_in_region(genes, region("1", 100L, 200L), promoter_length = 10000L)
    expect_equal(out$gene_id, "g1")
})

test_that("find_genes_in_region returns the GENE coordinates, not the query window", {
    # foverlaps(query, genes_chr) puts genes on the y side (unprefixed) and the
    # query on the i.* side, which the function then drops. Pinned because two
    # pipeline call sites read the wrong side of exactly this idiom
    # (sig_snps.R:152, genes_in_regions.R:174).
    genes <- gff("g1", "1", 1200L, 1300L)
    out <- find_genes_in_region(genes, region("1", 1000L, 2000L), promoter_length = 0L)

    expect_equal(out$start, 1200L)
    expect_equal(out$end,   1300L)
    expect_false(any(grepl("^i\\.", names(out))))
})

test_that("find_genes_in_region returns an UNTYPED empty table on every miss path", {
    # CURRENT BEHAVIOUR, and an inconsistency worth knowing about: the hit path
    # returns a schema, the four miss paths (fct_regions.R:68,69,84,94) return a
    # bare data.table() with zero columns. Callers must not assume column names.
    genes <- gff("g1", "1", 1200L, 1300L)

    misses <- list(
        find_genes_in_region(NULL, region("1", 1000L, 2000L)),
        find_genes_in_region(gff(character(), character(), integer(), integer()),
                             region("1", 1000L, 2000L)),
        find_genes_in_region(genes, NULL),
        find_genes_in_region(genes, region("9", 1000L, 2000L)),          # no such chr
        find_genes_in_region(genes, region("1", 8000L, 9000L), 0L)       # no overlap
    )
    for (out in misses) {
        expect_s3_class(out, "data.table")
        expect_equal(nrow(out), 0L)
    }
    expect_equal(ncol(misses[[1]]), 0L)
})

# ── pure argument builders ────────────────────────────────────────────────────

test_that(".module_to_hap_source maps only the two modules that have one", {
    expect_equal(.module_to_hap_source(MOD_GWAS),     "gwas")
    expect_equal(.module_to_hap_source(MOD_GEAXGWAS), "gea")
    expect_null(.module_to_hap_source(MOD_GEA))
    expect_null(.module_to_hap_source("nonsense"))
})

test_that(".build_trait_method_arg groups methods under their trait", {
    snps <- data.table::data.table(
        trait  = c("bio_1", "bio_1", "bio_2"),
        method = c("LFMM",  "EMMAX", "EMMAX")
    )
    expect_equal(.build_trait_method_arg(snps), "bio_1:EMMAX|LFMM,bio_2:EMMAX")
})

test_that(".build_trait_method_arg deduplicates and sorts methods", {
    snps <- data.table::data.table(
        trait  = c("bio_1", "bio_1", "bio_1"),
        method = c("LFMM",  "LFMM",  "EMMAX")
    )
    expect_equal(.build_trait_method_arg(snps), "bio_1:EMMAX|LFMM")
})

test_that(".build_trait_method_arg returns NULL on empty input", {
    expect_null(.build_trait_method_arg(NULL))
    expect_null(.build_trait_method_arg(
        data.table::data.table(trait = character(), method = character())))
})

test_that(".build_assoc_tables_arg builds one colon-separated entry per file", {
    files <- c("/p/GEA/tables/methods/EMMAX/EMMAX_pvalues_K3.tsv",
               "/p/GEA/tables/methods/LFMM/LFMM_pvalues_K3.tsv")
    expect_equal(
        .build_assoc_tables_arg(files),
        paste0("EMMAX:bonf_0.05:", files[1], ",LFMM:bonf_0.05:", files[2])
    )
})

test_that(".build_assoc_tables_arg SILENTLY DROPS files that do not match the pattern", {
    # Current behaviour (fct_regions.R:1113-1119). A sig_snps file or a renamed
    # p-value table vanishes from the regionplot arg with no warning.
    files <- c("/p/EMMAX_pvalues_K3.tsv",
               "/p/EMMAX_pvalues_K3_sig_snps_bonf_0.05.tsv",   # dropped
               "/p/something_else.tsv")                        # dropped
    expect_equal(.build_assoc_tables_arg(files),
                 "EMMAX:bonf_0.05:/p/EMMAX_pvalues_K3.tsv")
    expect_null(.build_assoc_tables_arg(files[2:3]))
    expect_null(.build_assoc_tables_arg(character()))
})

# ── compute_all_regions() — needs the mounted shared lib ──────────────────────

test_that("compute_all_regions clusters SNPs and re-derives traits/methods", {
    skip_if_not(exists("cluster_snps_to_regions", mode = "function"),
                "cluster_snps_to_regions not sourced (regions.R not on path)")

    # Two SNPs 100 bp apart on chr1 (one region), one far away on chr2.
    snps <- data.table::data.table(
        SNPID  = c("1:1000", "1:1100", "2:5000"),
        chr    = c("1", "1", "2"),
        pos    = c(1000L, 1100L, 5000L),
        pvalue = c(1e-8, 1e-6, 1e-7),
        method = c("EMMAX", "LFMM", "EMMAX"),
        trait  = c("bio_1", "bio_2", "bio_1")
    )

    out <- compute_all_regions(snps, distance = 1000L)

    expect_equal(nrow(out), 2L)
    expect_equal(names(out),
                 c("region_id", "chr", "start", "end", "length",
                   "snp_count", "snp_ids", "traits", "methods", "min_pvalue"))
    chr1 <- out[chr == "1"]
    expect_equal(chr1$snp_count, 2L)
    expect_equal(chr1$traits,  "bio_1,bio_2")
    expect_equal(chr1$methods, "EMMAX,LFMM")
    expect_equal(chr1$length, chr1$end - chr1$start)
})

test_that("compute_all_regions keeps the MINIMUM p-value when a SNP appears twice", {
    skip_if_not(exists("cluster_snps_to_regions", mode = "function"),
                "cluster_snps_to_regions not sourced (regions.R not on path)")

    snps <- data.table::data.table(
        SNPID  = c("1:1000", "1:1000"),
        chr    = c("1", "1"),
        pos    = c(1000L, 1000L),
        pvalue = c(1e-4, 1e-9),
        method = c("EMMAX", "LFMM"),
        trait  = c("bio_1", "bio_1")
    )

    out <- compute_all_regions(snps, distance = 1000L)

    expect_equal(nrow(out), 1L)
    expect_equal(out$snp_count, 1L)
    expect_equal(out$min_pvalue, 1e-9)
})

test_that("compute_all_regions returns the typed empty skeleton on empty input", {
    skip_if_not(exists("cluster_snps_to_regions", mode = "function"),
                "cluster_snps_to_regions not sourced (regions.R not on path)")

    out <- compute_all_regions(NULL)
    expect_s3_class(out, "data.table")
    expect_equal(nrow(out), 0L)
    expect_equal(names(out), names(.empty_regions()))
})
