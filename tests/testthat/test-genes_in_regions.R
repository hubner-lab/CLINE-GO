# scripts/R/lib/genes_in_regions.R
#
# Gene-region overlap plus exon/promoter SNP counting. Regions arriving here are
# ALREADY extended by the clumping distance (create_regions.R does that), so the
# overlap test is a plain interval intersection — no further padding.

# --- fixtures --------------------------------------------------------------

regions <- function(region_id, chr, start, end) {
    data.table::data.table(region_id = region_id, chr = as.character(chr),
                           start = as.integer(start), end = as.integer(end))
}

gff <- function(gene_id, chr, start, end, Name = gene_id) {
    data.table::data.table(chr = as.character(chr), start = as.integer(start),
                           end = as.integer(end), gene_id = gene_id, Name = Name)
}

# find_genes_for_regions() only reads gff_path when it needs exon features, i.e.
# when allsnps_dt is supplied.
write_gff_exons <- function(lines, envir = parent.frame()) {
    path <- withr::local_tempfile(fileext = ".gff3", .local_envir = envir)
    writeLines(c("##gff-version 3", lines), path)
    path
}

# --- gene-region overlap ---------------------------------------------------

test_that("a gene fully inside a region is found", {
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "1", 2000, 3000), "unused"))
    expect_identical(r$genes_per_region$gene_id, "g1")
    expect_identical(r$genes_per_region$region_id, "R1")
})

test_that("a gene overlapping only one edge of a region is found", {
    left  <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                          gff("g1", "1", 500, 1500), "unused"))
    right <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                          gff("g1", "1", 4500, 6000), "unused"))
    expect_identical(left$genes_per_region$gene_id, "g1")
    expect_identical(right$genes_per_region$gene_id, "g1")
})

test_that("a gene touching a region at exactly one base is found", {
    # Boundary is inclusive on both sides.
    lo <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                       gff("g1", "1", 500, 1000), "unused"))
    hi <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                       gff("g1", "1", 5000, 6000), "unused"))
    expect_identical(nrow(lo$genes_per_region), 1L)
    expect_identical(nrow(hi$genes_per_region), 1L)
})

test_that("a gene one base outside the region is NOT found", {
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "1", 5001, 6000), "unused"))
    expect_identical(nrow(r$genes_per_region), 0L)
})

test_that("a gene on a different chromosome is excluded", {
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "2", 2000, 3000), "unused"))
    expect_identical(nrow(r$genes_per_region), 0L)
})

test_that("one gene spanning two regions yields two rows, collapsed to one", {
    regs <- rbind(regions("R1", "1", 1000, 2000), regions("R2", "1", 2500, 4000))
    r <- quiet(find_genes_for_regions(regs, gff("g1", "1", 1500, 3000), "unused"))
    expect_identical(nrow(r$genes_per_region), 2L)
    expect_identical(nrow(r$genes_collapsed), 1L)
    expect_identical(r$genes_collapsed$region_id, "R1,R2")
})

# --- output shape ----------------------------------------------------------

test_that("region coordinates are dropped and gene coordinates renamed", {
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "1", 2000, 3000), "unused"))
    g <- r$genes_per_region
    expect_true(all(c("gene_start", "gene_end") %in% colnames(g)))
    expect_false(any(c("i.start", "i.end", "i.chr", "start", "end") %in% colnames(g)))
    expect_identical(g$gene_start, 2000L)
    expect_identical(g$gene_end, 3000L)
})

test_that("the key columns lead and the GFF attribute columns follow", {
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "1", 2000, 3000), "unused"))
    expect_identical(head(colnames(r$genes_per_region), 5),
                     c("region_id", "gene_id", "chr", "gene_start", "gene_end"))
    expect_true("Name" %in% colnames(r$genes_per_region))
    expect_identical(tail(colnames(r$genes_per_region), 4),
                     c("exon_snps", "promoter_snps",
                       "exon_snp_count", "promoter_snp_count"))
})

test_that("genes_collapsed is sorted by chr then gene_start", {
    regs <- regions(c("R1", "R2"), c("2", "1"), c(1000, 1000), c(9000, 9000))
    gfft <- rbind(gff("gB", "2", 2000, 3000), gff("gA", "1", 5000, 6000))
    r <- quiet(find_genes_for_regions(regs, gfft, "unused"))
    expect_identical(r$genes_collapsed$chr, c("1", "2"))
})

# --- exon / promoter counting ----------------------------------------------

test_that("without allsnps_dt the snp columns are empty strings and zero counts", {
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "1", 2000, 3000), "unused"))
    g <- r$genes_per_region
    expect_identical(g$exon_snps, "")
    expect_identical(g$promoter_snps, "")
    expect_identical(g$exon_snp_count, 0L)
    expect_identical(g$promoter_snp_count, 0L)
})

test_that("the promoter window extends UPSTREAM of gene_start and is floored at 1", {
    # Gene at 2000-3000, promoter_length 500 -> promoter is [1500, 2000]. A SNP at
    # 1600 is in it; one at 2500 (inside the gene body) and one at 1400 are not.
    path <- write_gff_exons(paste("1", "test", "exon", 9000, 9100, ".", "+", ".",
                                  "Parent=other.1", sep = "\t"))
    snps <- data.table::data.table(chr = "1", pos = c(1400L, 1600L, 2500L))
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1, 5000),
                                      gff("g1", "1", 2000, 3000), path,
                                      promoter_length = 500L, allsnps_dt = snps))
    expect_gt(r$genes_per_region$promoter_snp_count, 0L)

    # Nothing upstream at all -> no promoter SNPs, proving the window is not
    # symmetric around gene_start.
    snps_down <- data.table::data.table(chr = "1", pos = c(2500L, 3500L))
    r2 <- quiet(find_genes_for_regions(regions("R1", "1", 1, 5000),
                                       gff("g1", "1", 2000, 3000), path,
                                       promoter_length = 500L, allsnps_dt = snps_down))
    expect_identical(r2$genes_per_region$promoter_snp_count, 0L)
})

test_that("a gene starting below promoter_length still gets a valid window", {
    path <- write_gff_exons(paste("1", "test", "exon", 9000, 9100, ".", "+", ".",
                                  "Parent=other.1", sep = "\t"))
    snps <- data.table::data.table(chr = "1", pos = c(50L))
    # Gene at 100-500 with promoter_length 10000 would give a negative start; it
    # must be floored at 1, and the SNP at 50 must still be caught.
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1, 5000),
                                      gff("g1", "1", 100, 500), path,
                                      promoter_length = 10000L, allsnps_dt = snps))
    expect_gt(r$genes_per_region$promoter_snp_count, 0L)
})

test_that("exon counting uses the GFF's exon features, not the gene body", {
    path <- write_gff_exons(c(
        paste("1", "test", "exon", 2000, 2100, ".", "+", ".", "Parent=g1.1", sep = "\t")
    ))
    # 2050 is inside the exon; 2500 is inside the gene but intronic.
    inside <- data.table::data.table(chr = "1", pos = 2050L)
    intron <- data.table::data.table(chr = "1", pos = 2500L)
    r_in  <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                          gff("g1", "1", 2000, 3000), path,
                                          allsnps_dt = inside))
    r_out <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                          gff("g1", "1", 2000, 3000), path,
                                          allsnps_dt = intron))
    expect_gt(r_in$genes_per_region$exon_snp_count, 0L)
    expect_identical(r_out$genes_per_region$exon_snp_count, 0L)
})

test_that("an unreadable GFF degrades to zero exon SNPs instead of erroring", {
    snps <- data.table::data.table(chr = "1", pos = 2050L)
    r <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 5000),
                                      gff("g1", "1", 2000, 3000),
                                      "/nonexistent/path.gff3", allsnps_dt = snps))
    expect_identical(r$genes_per_region$exon_snp_count, 0L)
})

# --- .count_snps_in_features ------------------------------------------------

test_that(".count_snps_in_features returns the empty schema for an empty feature set", {
    snps <- data.table::data.table(chr = "1", pos = 100L)
    out  <- .count_snps_in_features(snps, data.table::data.table(
        chr = character(), start = integer(), end = integer(), gene_id = character()),
        col_name = "exon")
    expect_identical(nrow(out), 0L)
    expect_identical(colnames(out), c("gene_id", "exon_snps"))
})

test_that(".count_snps_in_features requires SNPs to be strictly within a feature", {
    feats <- data.table::data.table(chr = "1", start = 100L, end = 200L, gene_id = "g1")
    expect_identical(nrow(.count_snps_in_features(
        data.table::data.table(chr = "1", pos = 250L), feats, "exon")), 0L)
    expect_identical(nrow(.count_snps_in_features(
        data.table::data.table(chr = "1", pos = 150L), feats, "exon")), 1L)
})

# --- empty paths -----------------------------------------------------------

test_that("all three early exits return the 5-column empty table in BOTH slots", {
    # Schema mismatch worth pinning: the empty frame has no exon/promoter columns
    # at all, so it cannot be rbind()ed with a populated result.
    expected <- c("region_id", "gene_id", "chr", "gene_start", "gene_end")

    no_regions <- quiet(find_genes_for_regions(regions(character(), character(),
                                                       integer(), integer()),
                                               gff("g1", "1", 1, 2), "unused"))
    no_genes   <- quiet(find_genes_for_regions(regions("R1", "1", 1, 2),
                                               gff(character(), character(),
                                                   integer(), integer()), "unused"))
    no_overlap <- quiet(find_genes_for_regions(regions("R1", "1", 1000, 2000),
                                               gff("g1", "1", 8000, 9000), "unused"))

    for (r in list(no_regions, no_genes, no_overlap)) {
        expect_identical(nrow(r$genes_per_region), 0L)
        expect_identical(nrow(r$genes_collapsed), 0L)
        expect_identical(colnames(r$genes_per_region), expected)
    }
})

test_that("the caller's regions and GFF tables are not mutated", {
    regs <- regions("R1", 1L, 1000, 5000)   # numeric chr, coerced internally
    gfft <- gff("g1", 1L, 2000, 3000)
    quiet(find_genes_for_regions(regs, gfft, "unused"))
    expect_true(is.character(regs$chr))     # regions() already made it character
    expect_identical(colnames(gfft), c("chr", "start", "end", "gene_id", "Name"))
})
