# scripts/R/lib/enrichment_plots.R — the pure data shaping, not the plots.
#
# remap_enrich_genes() is the one target here that is fully pure and fully
# testable with no GO.db and no clusterProfiler call: it needs only an
# enrichResult S4 shell (DOSE is in the image), and it rewrites the "/"-joined
# geneID strings the GO tables and cnetplot labels are drawn from. A silent
# failure there relabels genes in a published figure.
#
# create_emapplot/create_cnetplot/create_dotplot/save_enrichment_plot and
# process_all_enrichment are NOT tested: they render, and the repo forbids reading
# plot content (rule 1). process_all_enrichment additionally counts INTENDED files
# rather than written ones (:197-206 vs the swallowing tryCatch at :112,126,138) —
# filed, not asserted here.

mk_enrich <- function(gene_ids = c("g1/g2", "g3"),
                      genes = c("g1", "g2", "g3")) {
    res_df <- data.frame(
        ID = paste0("GO:000", seq_along(gene_ids)),
        Description = paste0("term", seq_along(gene_ids)),
        geneID = gene_ids,
        Count = lengths(strsplit(gene_ids, "/")),
        stringsAsFactors = FALSE)
    new_enrich_result(res_df, gene = genes)
}

test_that("remap_enrich_genes is a no-op when the label map is NULL", {
    obj <- mk_enrich()
    got <- remap_enrich_genes(obj, NULL)
    expect_identical(got@result$geneID, c("g1/g2", "g3"))
    expect_identical(got@gene, c("g1", "g2", "g3"))
})

test_that("remap_enrich_genes relabels every gene in a slash-joined string", {
    obj <- mk_enrich()
    map <- c(g1 = "ALPHA", g2 = "BETA", g3 = "GAMMA")
    got <- remap_enrich_genes(obj, map)
    expect_identical(unname(got@result$geneID), c("ALPHA/BETA", "GAMMA"))
    expect_identical(unname(got@gene), c("ALPHA", "BETA", "GAMMA"))
})

test_that("remap_enrich_genes leaves an unmapped gene as its own id", {
    # A partial label map is the normal case: only some genes carry a Name.
    obj <- mk_enrich()
    got <- remap_enrich_genes(obj, c(g1 = "ALPHA"))
    expect_identical(unname(got@result$geneID), c("ALPHA/g2", "g3"))
    expect_identical(unname(got@gene), c("ALPHA", "g2", "g3"))
})

test_that("remap_enrich_genes preserves the slash separator and gene count", {
    obj <- mk_enrich(gene_ids = "g1/g2/g3")
    got <- remap_enrich_genes(obj, c(g2 = "BETA"))
    expect_identical(unname(got@result$geneID), "g1/BETA/g3")
    expect_identical(length(strsplit(got@result$geneID, "/")[[1]]), 3L)
})

test_that("build_label_map returns NULL for every documented no-op case", {
    d <- withr::local_tempdir()
    f <- file.path(d, "genes.tsv")
    data.table::fwrite(data.table::data.table(gene_id = c("g1", "g2"),
                                              Name = c("ALPHA", "BETA")),
                       f, sep = "\t")

    expect_null(build_label_map(NULL, "Name"))
    expect_null(build_label_map("NULL", "Name"))
    expect_null(build_label_map(file.path(d, "missing.tsv"), "Name"))
    # label_field == "gene_id" means "no relabelling wanted" — the identity map is
    # deliberately NULL rather than a wasted lookup table.
    expect_null(build_label_map(f, "gene_id"))
})

test_that("build_label_map falls back to NULL when the label column is absent", {
    d <- withr::local_tempdir()
    f <- file.path(d, "genes.tsv")
    data.table::fwrite(data.table::data.table(gene_id = c("g1"), Name = c("ALPHA")),
                       f, sep = "\t")
    expect_null(quiet(build_label_map(f, "description")))
})

test_that("build_label_map maps gene_id to the requested column, dropping blanks", {
    d <- withr::local_tempdir()
    f <- file.path(d, "genes.tsv")
    data.table::fwrite(data.table::data.table(
        gene_id = c("g1", "g2", "g3", "g1"),
        Name    = c("ALPHA", "", NA, "ALPHA")),   # blank, NA and a duplicate row
        f, sep = "\t")

    map <- quiet(build_label_map(f, "Name"))
    expect_identical(map, c(g1 = "ALPHA"))
})

test_that("top_regions_for_plotting returns NULL for every no-op case", {
    d <- withr::local_tempdir()
    f <- file.path(d, "regions.tsv")
    data.table::fwrite(data.table::data.table(
        region_id = "1_1-2", trait = "bio_1", snp_count = 5L, min_pvalue = 0.01),
        f, sep = "\t")

    expect_null(top_regions_for_plotting(NULL, 5))
    expect_null(top_regions_for_plotting("NULL", 5))
    expect_null(top_regions_for_plotting(file.path(d, "missing.tsv"), 5))
    # top_n <= 0 means "keep all", expressed as NULL.
    expect_null(top_regions_for_plotting(f, 0))
    expect_null(top_regions_for_plotting(f, -1))
})

test_that("top_regions_for_plotting takes the top N per trait by snp_count", {
    d <- withr::local_tempdir()
    f <- file.path(d, "regions.tsv")
    data.table::fwrite(data.table::data.table(
        region_id = c("1_a", "1_b", "1_c", "2_a", "2_b"),
        trait     = c("bio_1", "bio_1", "bio_1", "bio_2", "bio_2"),
        snp_count = c(3L, 9L, 5L, 1L, 7L),
        min_pvalue = c(0.01, 0.02, 0.03, 0.04, 0.05)),
        f, sep = "\t")

    got <- quiet(top_regions_for_plotting(f, 2))
    # Per trait, descending snp_count: bio_1 -> b(9), c(5); bio_2 -> b(7), a(1).
    expect_identical(got, c("1_b", "1_c", "2_b", "2_a"))
})

test_that("top_regions_for_plotting returns NULL for an empty regions table", {
    d <- withr::local_tempdir()
    f <- file.path(d, "regions.tsv")
    data.table::fwrite(data.table::data.table(
        region_id = character(), trait = character(),
        snp_count = integer(), min_pvalue = numeric()), f, sep = "\t")
    expect_null(top_regions_for_plotting(f, 5))
})
