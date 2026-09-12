# scripts/R/lib/enrichment.R — the GFF -> TERM2GENE shaping and the refusal
# guards. NOT the hypergeometric test.
#
# SCOPE, deliberately narrow. clusterProfiler::enricher() (:112-121), the
# all-NA-qvalue -> p.adjust substitution (:128-131) and the AnnotationDbi GO.db
# lookup (:70-74) are the statistics; they need the full annotation stack and are
# out of scope here. What IS testable, and what actually breaks, is the string
# work in front of them: the go_field lookahead regex, comma-splitting into
# TERM2GENE pairs, the empty/"NA" filtering, and the three early returns.
#
# THE MASKING TRAP, stated so nobody "fixes" it later: the "GO.db not available
# -> term2name NULL" branch (:83-84) is UNREACHABLE in cline-go:latest, because
# requireNamespace("GO.db") is TRUE there (GO.db 3.22.0 is installed). A test
# asserting that branch would pass only by accident of the harness, which is the
# same class of trap the app suite documents for qvalue
# (test-equivalence-app-pipeline.R:570-586). It is therefore not tested.
#
# SEARCH-PATH HYGIENE: build_term2gene_from_gff() reaches get_go_descriptions(),
# which runs library(GO.db) at CALL time (:69) — a global attach that also pulls
# AnnotationDbi, whose select() would mask dplyr::select() for every
# alphabetically-later file in the same test_dir() run. Every test here that
# calls the full function restores the search path afterwards.

# Attach-and-restore guard. Records what was attached before the call and detaches
# anything new, so one test cannot change name resolution for the rest of the run.
with_search_path_restored <- function(code) {
    before <- search()
    on.exit({
        added <- setdiff(search(), before)
        for (pkg in added) {
            try(detach(pkg, character.only = TRUE, unload = FALSE), silent = TRUE)
        }
    }, add = TRUE)
    force(code)
}

# A GFF3 whose chosen feature rows carry a GO field. `go_field` is configurable
# (GFF.go_field), so the fixture parameterises it rather than hardcoding "Ontology".
write_go_gff <- function(d, feature = "mRNA", go_field = "Ontology",
                         entries = NULL, name = "annot.gff3") {
    if (is.null(entries)) {
        # REAL GO ids on purpose. get_go_descriptions() (:66-80) passes them
        # straight to AnnotationDbi::select() with no validation, and that call
        # HARD-ERRORS when none of the keys are valid ("None of the keys entered
        # are valid keys for 'GOID'"). Invented ids like "GO:0001" therefore make
        # this fixture kill the function under test rather than exercise it —
        # which is itself a filed defect, quarantined in test-known-bugs.R.
        entries <- list(
            list(id = "g1", go = "GO:0008150,GO:0003674"),
            list(id = "g2", go = "GO:0003674"),
            list(id = "g3", go = "GO:0005575,GO:0008150,GO:0003674")
        )
    }
    p <- file.path(d, name)
    rows <- vapply(seq_along(entries), function(i) {
        e <- entries[[i]]
        attrs <- paste0("ID=", e$id, ";Name=", toupper(e$id))
        if (!is.null(e$go)) attrs <- paste0(attrs, ";", go_field, "=", e$go)
        paste("1", "src", feature, i * 1000L, i * 1000L + 500L, ".", "+", ".",
              attrs, sep = "\t")
    }, character(1))
    writeLines(c("##gff-version 3", rows), p)
    p
}

test_that("build_term2gene_from_gff refuses an unset go_field", {
    # Three spellings of "not configured" — all must return NULL rather than
    # building a bogus regex like "(?<==)[^;]+".
    for (gf in list(NULL, "NULL", "")) {
        expect_null(quiet(build_term2gene_from_gff("/nonexistent.gff3", "mRNA", gf)))
    }
})

test_that("build_term2gene_from_gff returns NULL when the feature type is absent", {
    d <- withr::local_tempdir()
    p <- write_go_gff(d, feature = "mRNA")
    # GFF.feature is a config value; asking for one the file does not contain must
    # yield empty, not an error (the Tier 1 contract for gff_parsing.R too).
    expect_null(quiet(build_term2gene_from_gff(p, "gene", "Ontology")))
})

test_that("build_term2gene_from_gff returns NULL when no feature carries a GO term", {
    d <- withr::local_tempdir()
    p <- write_go_gff(d, entries = list(list(id = "g1", go = NULL),
                                        list(id = "g2", go = NULL)))
    expect_null(quiet(build_term2gene_from_gff(p, "mRNA", "Ontology")))
})

test_that("build_term2gene_from_gff splits comma-separated GO terms into pairs", {
    d <- withr::local_tempdir()
    p <- write_go_gff(d)
    got <- with_search_path_restored(
        quiet(build_term2gene_from_gff(p, "mRNA", "Ontology")))

    expect_named(got, c("term2gene", "term2name", "all_genes_with_go"))
    t2g <- got$term2gene
    expect_identical(names(t2g), c("term", "gene"))
    # 2 + 1 + 3 annotations = 6 pairs, one row per (term, gene).
    expect_identical(nrow(t2g), 6L)
    expect_setequal(unique(t2g$term),
                    c("GO:0008150", "GO:0003674", "GO:0005575"))
    expect_setequal(t2g$gene[t2g$term == "GO:0003674"], c("g1", "g2", "g3"))
    expect_setequal(got$all_genes_with_go, c("g1", "g2", "g3"))
})

test_that("build_term2gene_from_gff reads the go_field NAME it is given", {
    # The regex is a lookahead on the field name, so a project using a different
    # attribute key must work and the wrong key must find nothing.
    d <- withr::local_tempdir()
    p <- write_go_gff(d, go_field = "go_terms")

    got <- with_search_path_restored(
        quiet(build_term2gene_from_gff(p, "mRNA", "go_terms")))
    expect_identical(nrow(got$term2gene), 6L)

    expect_null(quiet(build_term2gene_from_gff(p, "mRNA", "Ontology")))
})

test_that("build_term2gene_from_gff drops the literal string NA and empties", {
    d <- withr::local_tempdir()
    p <- write_go_gff(d, entries = list(
        list(id = "g1", go = "GO:0008150"),
        list(id = "g2", go = "NA"),          # literal 'NA', not a missing value
        list(id = "g3", go = "")))
    got <- with_search_path_restored(
        quiet(build_term2gene_from_gff(p, "mRNA", "Ontology")))

    expect_identical(nrow(got$term2gene), 1L)
    expect_identical(got$all_genes_with_go, "g1")
})

test_that("run_enrichment_for_region refuses a region with fewer than 2 GO genes", {
    # The refusal happens BEFORE clusterProfiler is touched (:100-106), which is
    # why it is testable here at all.
    t2g <- data.frame(term = c("GO:0008150", "GO:0008150"), gene = c("g1", "g2"),
                      stringsAsFactors = FALSE)
    bg <- c("g1", "g2")

    # One overlapping gene.
    expect_null(quiet(run_enrichment_for_region(
        "1_100-200", "bio_1", c("g1", "gX"), t2g, NULL, bg)))
    # No overlapping gene, even though the region has plenty.
    expect_null(quiet(run_enrichment_for_region(
        "1_100-200", "bio_1", c("gX", "gY", "gZ"), t2g, NULL, bg)))
    # No genes at all.
    expect_null(quiet(run_enrichment_for_region(
        "1_100-200", "bio_1", character(0), t2g, NULL, bg)))
})

test_that("save_enrichment_result writes the documented column order", {
    # The TSV schema is a downstream contract (the Shiny GO table reads it), and
    # it is built with setcolorder (:174-177) — easy to break silently.
    d <- withr::local_tempdir()
    res_df <- data.frame(
        ID = c("GO:0001", "GO:0002"),
        Description = c("alpha", "beta"),
        GeneRatio = c("2/3", "1/3"), BgRatio = c("3/9", "2/9"),
        pvalue = c(0.001, 0.04), p.adjust = c(0.002, 0.08),
        qvalue = c(0.002, 0.08), geneID = c("g1/g2", "g3"),
        Count = c(2L, 1L), stringsAsFactors = FALSE)
    obj <- new_enrich_result(res_df, gene = c("g1", "g2", "g3"))

    out <- quiet(save_enrichment_result(
        list(region_id = "1_100-200", trait = "bio_1", enrich_obj = obj),
        file.path(d, "tables"), file.path(d, "inter")))

    expect_identical(names(out),
        c("region_id", "trait", "GO_id", "description", "gene_ratio", "bg_ratio",
          "pvalue", "p_adjust", "qvalue", "gene_ids", "gene_count"))

    # Per-trait subdirectories, and both artefacts written.
    tsv <- file.path(d, "tables", "bio_1", "Region_1_100-200_enrichment.tsv")
    qsf <- file.path(d, "inter",  "bio_1", "Region_1_100-200_enrichment.qs")
    expect_true(file.exists(tsv))
    expect_true(file.exists(qsf))
    expect_gt(file.size(tsv), 0)

    back <- data.table::fread(tsv)
    expect_identical(nrow(back), 2L)
    expect_identical(unique(back$region_id), "1_100-200")
    expect_identical(unique(back$trait), "bio_1")
})
