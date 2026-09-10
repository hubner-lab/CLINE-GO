# scripts/R/lib/gff_parsing.R
#
# Both readers shell out to `grep -v '#'`, so every fixture has to be a real file
# on disk — an in-memory data.table will not exercise them.

# --- fixtures --------------------------------------------------------------

write_gff <- function(lines, envir = parent.frame()) {
    path <- withr::local_tempfile(fileext = ".gff3", .local_envir = envir)
    writeLines(c("##gff-version 3", lines), path)
    path
}

gff_line <- function(chr, feature, start, end, attrs, strand = "+") {
    paste(chr, "test", feature, start, end, ".", strand, ".", attrs, sep = "\t")
}

# --- extract_gene_id -------------------------------------------------------

test_that("extract_gene_id prefers Parent= over ID=", {
    expect_identical(extract_gene_id("ID=gene1.1_exon_2;Parent=gene1.1"), "gene1")
})

test_that("extract_gene_id falls back to ID= when there is no Parent=", {
    expect_identical(extract_gene_id("ID=gene1;Name=foo"), "gene1")
})

test_that("extract_gene_id strips isoform and exon suffixes", {
    expect_identical(extract_gene_id("ID=gene1.1"), "gene1")
    expect_identical(extract_gene_id("ID=gene1.12_exon_3"), "gene1")
    # Only a trailing .N is an isoform — an internal dot is part of the name.
    expect_identical(extract_gene_id("ID=HORVU.MOREX.r3.1HG0000010"),
                     "HORVU.MOREX.r3.1HG0000010")
})

test_that("extract_gene_id is vectorised", {
    out <- extract_gene_id(c("ID=a.1", "Parent=b.2", "ID=c"))
    expect_identical(out, c("a", "b", "c"))
})

test_that("extract_gene_id returns NA when neither key is present", {
    expect_true(is.na(extract_gene_id("Name=nothing;Note=here")))
})

# --- clean_attr_value ------------------------------------------------------

test_that("clean_attr_value strips the leading key=", {
    expect_identical(clean_attr_value("ID=gene1"), "gene1")
    expect_identical(clean_attr_value(c("Name=foo", "Note=bar")), c("foo", "bar"))
})

test_that("clean_attr_value strips only the FIRST key= it finds", {
    expect_identical(clean_attr_value("Note=a=b"), "a=b")
})

# --- read_gff --------------------------------------------------------------

test_that("read_gff returns chr as character, start/end as integer", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 200, "ID=gene1;Name=alpha"),
        gff_line("2", "gene", 300, 400, "ID=gene2;Name=beta")
    ))
    g <- quiet(read_gff(path, "gene"))
    expect_true(is.character(g$chr))
    expect_true(is.integer(g$start))
    expect_true(is.integer(g$end))
    expect_identical(g$chr, c("1", "2"))
    expect_identical(g$start, c(100L, 300L))
})

test_that("read_gff keeps only the requested feature type", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 200, "ID=gene1"),
        gff_line("1", "mRNA", 100, 200, "ID=gene1.1;Parent=gene1"),
        gff_line("1", "exon", 100, 150, "ID=gene1.1_exon_1;Parent=gene1.1")
    ))
    expect_identical(nrow(quiet(read_gff(path, "gene"))), 1L)
    expect_identical(nrow(quiet(read_gff(path, "mRNA"))), 1L)
})

test_that("read_gff splits attributes into one column per key", {
    path <- write_gff(gff_line("1", "gene", 100, 200,
                               "ID=gene1;Name=alpha;Ontology_term=GO:0006355"))
    g <- quiet(read_gff(path, "gene"))
    expect_true(all(c("ID", "Name", "Ontology_term") %in% colnames(g)))
    expect_identical(g$Name, "alpha")
    expect_identical(g$Ontology_term, "GO:0006355")
    expect_identical(g$gene_id, "gene1")
})

test_that("read_gff turns literal 'NA' and empty attribute values into NA", {
    path <- write_gff(gff_line("1", "gene", 100, 200, "ID=gene1;Name=NA;Note="))
    g <- quiet(read_gff(path, "gene"))
    expect_true(is.na(g$Name))
    expect_true(is.na(g$Note))
})

test_that("read_gff infers attribute names from the FIRST matching row only", {
    # Documented behaviour, not a bug to fix here: a GFF with heterogeneous
    # attribute keys is split against row 1's key list, so row 2's extra key is
    # never given a column of its own.
    path <- write_gff(c(
        gff_line("1", "gene", 100, 200, "ID=gene1;Name=alpha"),
        gff_line("1", "gene", 300, 400, "ID=gene2;Name=beta;Note=extra")
    ))
    g <- quiet(read_gff(path, "gene"))
    # separate() expands `description` in place, so the attribute columns land
    # where it was and gene_id (added before the split) trails them.
    expect_identical(colnames(g),
                     c("chr", "start", "end", "ID", "Name", "gene_id"))
    expect_false("Note" %in% colnames(g))
})

test_that("read_gff on a missing feature type warns and returns the 4-column empty frame", {
    path <- write_gff(gff_line("1", "gene", 100, 200, "ID=gene1;Name=alpha"))
    g <- expect_message(read_gff(path, "tRNA"), "No features of type")
    expect_identical(nrow(g), 0L)
    # Schema differs from the populated case: no attribute columns at all.
    expect_identical(colnames(g), c("chr", "start", "end", "gene_id"))
})

test_that("read_gff preserves a non-numeric chromosome name verbatim", {
    # The pipeline strips the 'chr' prefix with sed in processing.smk, NOT here —
    # read_gff must not second-guess whatever names it is handed.
    path <- write_gff(gff_line("chr2H", "gene", 100, 200, "ID=gene1"))
    expect_identical(quiet(read_gff(path, "gene"))$chr, "chr2H")
})

# --- read_gff_exons --------------------------------------------------------

test_that("read_gff_exons loads exon features with their parent gene id", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 500, "ID=gene1"),
        gff_line("1", "exon", 100, 150, "ID=gene1.1_exon_1;Parent=gene1.1"),
        gff_line("1", "exon", 400, 500, "ID=gene1.1_exon_2;Parent=gene1.1")
    ))
    e <- quiet(read_gff_exons(path))
    expect_identical(colnames(e), c("chr", "start", "end", "gene_id"))
    expect_identical(nrow(e), 2L)
    expect_identical(unique(e$gene_id), "gene1")
})

test_that("read_gff_exons falls back to CDS when there is no exon feature", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 500, "ID=gene1"),
        gff_line("1", "CDS", 120, 180, "ID=cds1;Parent=gene1.1")
    ))
    e <- expect_message(read_gff_exons(path), "falling back to CDS")
    expect_identical(nrow(e), 1L)
    expect_identical(e$start, 120L)
})

test_that("read_gff_exons returns an empty table when neither exon nor CDS exists", {
    path <- write_gff(gff_line("1", "gene", 100, 500, "ID=gene1"))
    e <- quiet(read_gff_exons(path))
    expect_identical(nrow(e), 0L)
})
