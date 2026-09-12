# Correct-behaviour assertions for defects that are KNOWN and DELIBERATELY NOT
# FIXED yet.
#
# Every block here starts with skip(), so the suite stays green and usable as a
# merge gate. They are written as the behaviour the code SHOULD have — never as
# an assertion of the buggy output — so that fixing a bug means deleting one
# skip() line and getting a real regression test for free. Do not "fix" a test
# here by weakening it to match current output.
#
# Each skip() names where the defect is filed: docs/pipeline_improvement_requests.md
# and/or the ADAPTOGENE pipeline dossier's ## Findings / plans/.

# ---------------------------------------------------------------------------
# 1. sig_snps.R:152,159 — overlap_traits / overlap_snps read the wrong side of
#    foverlaps. CONFIRMED in every *_sig_snps_*.tsv the pipeline has written.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-10), Findings [sci][med].
# ---------------------------------------------------------------------------

test_that("overlap_traits names the OTHER trait, not the row's own", {
    skip("known bug: sig_snps.R:152 uses i.trait (query side) — filed 2026-09-10")

    dt <- data.table::data.table(
        SNPID = c("1:100", "1:110"), chr = "1", pos = c(100L, 110L),
        pvalue = 1e-8, pval_threshold = 1e-6, trait = c("bio_1", "bio_2")
    )
    r <- annotate_cross_trait_overlaps(dt, 1000L)
    expect_identical(r[trait == "bio_1"]$overlap_traits, "bio_2")
    expect_identical(r[trait == "bio_2"]$overlap_traits, "bio_1")
})

test_that("overlap_snps gives the NEIGHBOUR's chr:pos, not pos minus distance", {
    skip("known bug: sig_snps.R:159 uses i.s (window start) — filed 2026-09-10")

    dt <- data.table::data.table(
        SNPID = c("1:100", "1:110"), chr = "1", pos = c(100L, 110L),
        pvalue = 1e-8, pval_threshold = 1e-6, trait = c("bio_1", "bio_2")
    )
    r <- annotate_cross_trait_overlaps(dt, 1000L)
    expect_identical(r[trait == "bio_1"]$overlap_snps, "1:110")
    expect_identical(r[trait == "bio_2"]$overlap_snps, "1:100")
})

# ---------------------------------------------------------------------------
# 2. sig_snps.R:44 — diag_list[[trait]] <<- ... inside parallel::mclapply. The
#    superassignment dies with the forked worker, so diagnostics is empty at
#    cpu >= 2. Latent: find_sig_snps.R never reads the slot.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-10), Findings [bug][low].
# ---------------------------------------------------------------------------

test_that("diagnostics is populated at cpu = 2, not only at cpu = 1", {
    skip("known bug: sig_snps.R:44 superassigns inside mclapply — filed 2026-09-10")

    dt <- data.table::data.table(
        SNPID = paste0("1:", 1:3), chr = "1", pos = 1:3,
        bio_1 = c(1e-9, 0.5, 0.5), bio_2 = c(0.5, 1e-9, 0.5)
    )
    one <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 1))
    two <- quiet(find_significant_snps_per_trait(dt, "custom", 1e-6, 2))
    expect_identical(nrow(one$diagnostics), 2L)
    expect_identical(nrow(two$diagnostics), nrow(one$diagnostics))
    expect_setequal(two$diagnostics$trait, c("bio_1", "bio_2"))
})

# ---------------------------------------------------------------------------
# 3. region_distance.R:110,120,121 — ld_table[group == group & ...]. The bare
#    data.table column masks the function argument, so the filter is a tautology
#    and the group= parameter is inert.
#    Filed: plans/2026-08-18-bug-resolve-clumping-distance-group-filter-i.md
# ---------------------------------------------------------------------------

test_that("resolve_clumping_distance honours the requested group", {
    skip("known bug: region_distance.R:110 `group == group` is a tautology — pitched 2026-08-18")

    path <- withr::local_tempfile(fileext = ".tsv")
    ld <- data.table::data.table(
        # Very different C_hat, so the two groups imply very different distances.
        group = c("cluster1", "All"), scope = "genome_wide", method = "hill_weir",
        C_hat = c(5e-2, 1e-4), n_samples = 50, max_dist_bp = 5e5,
        half_decay_bp = NA_real_, r2_02_bp = NA_real_, r2_intercept = NA_real_
    )
    data.table::fwrite(ld, path, sep = "\t")

    expected_all <- as.integer(round(
        resolve_row_distance(as.list(ld[group == "All"][1]), 0.2)))
    expect_identical(
        suppressMessages(resolve_clumping_distance("auto_genome_wide", path,
                                                   r2_threshold = 0.2, group = "All")),
        expected_all
    )
})

test_that("auto_per_chromosome does not mix rows from other groups", {
    skip("known bug: region_distance.R:120-121 same tautology — pitched 2026-08-18")

    path <- withr::local_tempfile(fileext = ".tsv")
    ld <- data.table::data.table(
        group = c("All", "All", "cluster1"),
        scope = c("genome_wide", "1", "1"),
        method = "hill_weir",
        C_hat = c(1e-3, 1e-3, 5e-2),   # cluster1's chr 1 decays much faster
        n_samples = 50, max_dist_bp = 5e5,
        half_decay_bp = NA_real_, r2_02_bp = NA_real_, r2_intercept = NA_real_
    )
    data.table::fwrite(ld, path, sep = "\t")

    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path,
                                                    0.2, group = "All"))
    expect_identical(d[["1"]],
                     as.integer(round(invert_hill_weir(1e-3, 50, 0.2, 5e5))))
})

# ---------------------------------------------------------------------------
# 4. genes_in_regions.R:95 — the promoter window is always upstream of
#    gene_start, ignoring strand. For a '-' strand gene the promoter lies
#    downstream of gene_end.
#    Filed: plans/2026-08-18-sci-promoter-snp-windows-ignore-gene-strand.md
#    Note: fixing this also needs read_gff() to KEEP the strand column (GFF3
#    field 7), which it currently drops — see gff_parsing.R:35-38.
# ---------------------------------------------------------------------------

test_that("a minus-strand gene's promoter sits downstream of gene_end", {
    skip("known bug: genes_in_regions.R:95 promoter window is strand-blind — pitched 2026-08-18")

    exon_path <- withr::local_tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
                 paste("1", "t", "exon", 9000, 9100, ".", "-", ".",
                       "Parent=other.1", sep = "\t")), exon_path)

    regions <- data.table::data.table(region_id = "R1", chr = "1",
                                      start = 1L, end = 5000L)
    genes <- data.table::data.table(chr = "1", start = 2000L, end = 3000L,
                                    gene_id = "g1", strand = "-")

    # For a '-' strand gene transcribed right-to-left, the promoter is
    # [gene_end, gene_end + promoter_length] = [3000, 3500]. A SNP at 3200 is in
    # it; a SNP at 1800 (upstream of gene_start) is NOT.
    upstream_snp   <- data.table::data.table(chr = "1", pos = 1800L)
    downstream_snp <- data.table::data.table(chr = "1", pos = 3200L)

    r_down <- quiet(find_genes_for_regions(regions, genes, exon_path,
                                           promoter_length = 500L,
                                           allsnps_dt = downstream_snp))
    r_up   <- quiet(find_genes_for_regions(regions, genes, exon_path,
                                           promoter_length = 500L,
                                           allsnps_dt = upstream_snp))
    expect_identical(r_down$genes_per_region$promoter_snp_count, 1L)
    expect_identical(r_up$genes_per_region$promoter_snp_count, 0L)
})

# ---------------------------------------------------------------------------
# 5. genes_in_regions.R:174 — .count_snps_in_features() builds snp_id from the
#    plain (feature-side) `s` column, so exon_snps / promoter_snps report the
#    FEATURE start and the counts count features, not SNPs.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-10), Findings [sci][med].
#    This is CLAUDE.md's open "validate exon/promoter SNP counting" TODO.
# ---------------------------------------------------------------------------

test_that(".count_snps_in_features lists SNP positions, not feature starts", {
    skip("known bug: genes_in_regions.R:174 uses plain `s` (feature side) — filed 2026-09-10")

    snps  <- data.table::data.table(chr = "1", pos = c(120L, 130L))
    feats <- data.table::data.table(chr = "1", start = 100L, end = 200L,
                                    gene_id = "g1")
    out <- .count_snps_in_features(snps, feats, col_name = "exon")
    expect_identical(out$exon_snps, "1:120,1:130")
})

test_that("exon_snp_count counts SNPs, not exons hit", {
    skip("known bug: genes_in_regions.R:174 — filed 2026-09-10")

    exon_path <- withr::local_tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
                 paste("1", "t", "exon", 2000, 2200, ".", "+", ".",
                       "Parent=g1.1", sep = "\t")), exon_path)

    r <- quiet(find_genes_for_regions(
        data.table::data.table(region_id = "R1", chr = "1", start = 1L, end = 5000L),
        data.table::data.table(chr = "1", start = 2000L, end = 3000L, gene_id = "g1"),
        exon_path,
        allsnps_dt = data.table::data.table(chr = "1", pos = c(2050L, 2100L, 2150L))
    ))
    expect_identical(r$genes_per_region$exon_snp_count, 3L)
})

# ---------------------------------------------------------------------------
# 6. regions.R:94 — min(grp$min_pvalue, na.rm = TRUE) on an all-NA cluster
#    returns Inf (with a warning) rather than NA.
#    Filed: plans/2026-08-18-bug-cluster-snps-to-regions-reports-inf-not.md
# ---------------------------------------------------------------------------

test_that("a region with no p-values reports min_pvalue = NA, not Inf", {
    skip("known bug: regions.R:94 min(na.rm=TRUE) yields Inf — pitched 2026-08-18")

    dt <- data.table::data.table(SNPID = "1:5000", chr = "1", pos = 5000L)
    r  <- suppressWarnings(cluster_snps_to_regions(dt, 1000L))
    expect_true(is.na(r$min_pvalue))
    expect_false(is.infinite(r$min_pvalue))
})

# ---------------------------------------------------------------------------
# 7. sig_snps.R:120-124 — the guard is `is.null(sig_snps) || nrow(...) == 0`,
#    but the branch it guards then subsets sig_snps with `[, := ]`, which errors
#    on NULL. The zero-row case works; the NULL case does not.
# ---------------------------------------------------------------------------

test_that("annotate_cross_trait_overlaps handles NULL the way its own guard implies", {
    skip("known bug: sig_snps.R:121 subsets NULL after guarding for it — not yet filed separately")

    expect_no_error(annotate_cross_trait_overlaps(NULL, 1000L))
})

# ---------------------------------------------------------------------------
# 8. manhattan_utils.R:85,95,104 — three palette getters index with
#    `x[1:length(y)]`. At length 0 that is `x[1:0]` = `x[c(1, 0)]`, which returns
#    ONE element instead of none, and setNames() then labels it NA. An empty
#    trait/method/region set is ordinary (a trait with no significant SNPs), so
#    these hand a one-colour palette with an NA name to ggplot2 rather than an
#    empty one. get_chr_colors (:76) is CORRECT — rep(length.out = 0) handles 0.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-12).
#    Fix in all three: seq_len(length(y)) instead of 1:length(y).
# ---------------------------------------------------------------------------

test_that("get_trait_colors returns an empty palette for no traits", {
    skip("known bug: manhattan_utils.R:85 uses 1:length(traits) — filed 2026-09-12")

    got <- get_trait_colors(character(0))
    expect_length(got, 0L)
})

test_that("get_method_shapes returns an empty vector for no methods", {
    skip("known bug: manhattan_utils.R:104 uses 1:length(methods) — filed 2026-09-12")

    got <- get_method_shapes(character(0))
    expect_length(got, 0L)
    expect_false(any(is.na(names(got))))
})

test_that("get_region_colors returns an empty vector for zero regions", {
    skip("known bug: manhattan_utils.R:95 uses colors[1:n_regions] — filed 2026-09-12")

    expect_length(get_region_colors(0), 0L)
})

# ---------------------------------------------------------------------------
# 9. enrichment.R:66-80 — get_go_descriptions() passes its keys straight to
#    AnnotationDbi::select() with no validation. That call HARD-ERRORS when none
#    of the keys resolve ("None of the keys entered are valid keys for 'GOID'"),
#    so a GFF whose go_field holds a typo, a non-GO string or an obsolete term
#    kills build_term2gene_from_gff() — and with it the whole enrichment rule —
#    instead of degrading to GO-ID labels the way the "GO.db not available" path
#    does. Found 2026-09-12 while writing test-enrichment.R, whose fixture had to
#    switch to real GO ids to get past it.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-12).
# ---------------------------------------------------------------------------

test_that("build_term2gene_from_gff survives a GFF with unresolvable GO ids", {
    skip("known bug: enrichment.R:69-74 select() errors on all-invalid keys — filed 2026-09-12")

    d <- withr::local_tempdir()
    p <- file.path(d, "annot.gff3")
    writeLines(c("##gff-version 3",
                 paste("1", "src", "mRNA", "1000", "1500", ".", "+", ".",
                       "ID=g1;Ontology=GO:9999999", sep = "\t"),
                 paste("1", "src", "mRNA", "2000", "2500", ".", "+", ".",
                       "ID=g2;Ontology=GO:9999999", sep = "\t")),
               p)

    # The TERM2GENE half needs no GO.db at all, so an unresolvable id should cost
    # only the human-readable names.
    got <- suppressMessages(build_term2gene_from_gff(p, "mRNA", "Ontology"))
    expect_identical(nrow(got$term2gene), 2L)
    expect_setequal(got$all_genes_with_go, c("g1", "g2"))
})

# ---------------------------------------------------------------------------
# 8. invariants.R:371 — check_summary_counts()'s `step` parameter is shadowed by
#    the data.table column of the same name, so `step == step` is TRUE for every
#    row and the step filter does nothing. The sibling at :319 gets it right
#    (`step == step_`). With a metric name present under two steps the match is
#    length 2, `next` fires, and the count goes unchecked in silence.
#    Found 2026-09-12 while building the check_invariants.R wrapper fixture.
#    Filed: docs/pipeline_improvement_requests.md (2026-09-12).
# ---------------------------------------------------------------------------

test_that("check_summary_counts scopes its comparison to the step it was given", {
    skip("known bug: invariants.R:371 `step == step` self-compares the column — filed 2026-09-12")

    # The same metric name under two steps: exactly what write_summary.R writes.
    s <- data.table::data.table(
        step   = c("gea", "gwas"),
        metric = c("selected_snps_total", "selected_snps_total"),
        value  = c("8", "3")
    )
    # gea really has 8, so scoped to gea this is clean...
    expect_identical(nrow(check_summary_counts(s, "gea", list(selected_snps_total = 8L))), 0L)
    # ...and scoped to gwas, 8 is wrong: gwas states 3.
    v <- check_summary_counts(s, "gwas", list(selected_snps_total = 8L))
    expect_identical(nrow(v), 1L)
    expect_identical(v$check, "summary_count_disagrees_with_table")
})
