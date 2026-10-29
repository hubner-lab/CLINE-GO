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

# 4de2a0: attributes are read by KEY. The names used to come from row 1 and every row was
# split against them by position.
test_that("read_gff gives a key that only a later row carries its own column", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 200, "ID=gene1;Name=alpha"),
        gff_line("1", "gene", 300, 400, "ID=gene2;Name=beta;Note=extra")
    ))
    g <- quiet(read_gff(path, "gene"))
    # The attribute columns take description's place, in order of first appearance, and
    # gene_id trails them (the column order the positional split had).
    expect_identical(colnames(g),
                     c("chr", "start", "end", "strand", "ID", "Name", "Note", "gene_id"))
    expect_identical(g$Note, c(NA, "extra"))
})

test_that("read_gff matches each value to its own key when a row omits one", {
    # The d150 transcript table: a row without em_EC had its em_PFAMs value read as em_EC,
    # and a key row 1 lacked (em_GOs) was dropped.
    path <- write_gff(c(
        gff_line("1", "transcript", 100, 200, "ID=g1.t1;em_desc=kinase;em_EC=2.7.11.1;em_PFAMs=Pkinase"),
        gff_line("1", "transcript", 300, 400, "ID=g2.t1;em_desc=BURP;em_PFAMs=BURP;em_GOs=GO:0005575"),
        gff_line("1", "transcript", 500, 600, "ID=g3.t1;em_GOs=GO:0008150;em_desc=none")
    ))
    g <- quiet(read_gff(path, "transcript"))
    expect_identical(colnames(g), c("chr", "start", "end", "strand", "ID", "em_desc", "em_EC",
                                    "em_PFAMs", "em_GOs", "gene_id"))
    expect_identical(g$em_EC, c("2.7.11.1", NA, NA))
    expect_identical(g$em_PFAMs, c("Pkinase", "BURP", NA))
    expect_identical(g$em_GOs, c(NA, "GO:0005575", "GO:0008150"))
    expect_identical(g$em_desc, c("kinase", "BURP", "none"))
    expect_identical(g$gene_id, c("g1.t1", "g2.t1", "g3.t1"))
})

test_that("read_gff keeps a key's first value and skips pieces that are no key=value", {
    path <- write_gff(gff_line("1", "gene", 100, 200, "ID=gene1;;Note=a=b;Note=second;Flag;=x"))
    g <- quiet(read_gff(path, "gene"))
    expect_identical(colnames(g), c("chr", "start", "end", "strand", "ID", "Note", "gene_id"))
    expect_identical(g$Note, "a=b")
})

test_that("an attribute named like a core column never overwrites it", {
    # An Ensembl gene row carries gene_id=; it replaced the id every gene join keys on.
    path <- write_gff(gff_line("1", "gene", 100, 200,
                               "ID=gene:ENSG01;gene_id=ENSG01;chr=x;Name=a"))
    g <- quiet(read_gff(path, "gene"))
    expect_identical(g$gene_id, "gene:ENSG01")
    expect_identical(g$chr, "1")
    expect_identical(g$gene_id_attr, "ENSG01")
    expect_identical(g$chr_attr, "x")
    expect_identical(colnames(g), c("chr", "start", "end", "strand", "ID", "gene_id_attr",
                                    "chr_attr", "Name", "gene_id"))
})

test_that("read_gff on a missing feature type warns and returns the empty frame", {
    path <- write_gff(gff_line("1", "gene", 100, 200, "ID=gene1;Name=alpha"))
    g <- expect_message(read_gff(path, "tRNA"), "No features of type")
    expect_identical(nrow(g), 0L)
    # Schema differs from the populated case: no attribute columns at all.
    expect_identical(colnames(g), c("chr", "start", "end", "strand", "gene_id"))
})

# GFF3 field 7 is kept because the promoter window depends on it: on the - strand
# the transcription start is gene_end, not gene_start (genes_in_regions.R).
test_that("read_gff keeps GFF3 field 7 verbatim", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 200, "ID=gene1", strand = "+"),
        gff_line("1", "gene", 300, 400, "ID=gene2", strand = "-"),
        gff_line("1", "gene", 500, 600, "ID=gene3", strand = ".")
    ))
    g <- quiet(read_gff(path, "gene"))
    expect_identical(g$strand, c("+", "-", "."))
})

test_that("an attribute literally named 'strand' does not shadow field 7", {
    path <- write_gff(gff_line("1", "gene", 100, 200, "ID=gene1;strand=bogus",
                               strand = "-"))
    g <- quiet(read_gff(path, "gene"))
    expect_identical(g$strand, "-")
    expect_identical(g$strand_attr, "bogus")
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

# Finding 7d89be: exon_snp_count joins exon gene ids to the genes table's gene ids, so every
# exon id must be one read_gff() produces, for GFF.feature = gene AND for the transcript level.
# Stripping ".N" off the exon's Parent= only met the gene side when transcript and gene shared
# a namespace; these are the conventions that do not.
exon_ids_join <- function(path, features) {
    e <- quiet(read_gff_exons(path))
    for (f in features) {
        g <- quiet(read_gff(path, f))
        expect_true(all(e$gene_id %in% g$gene_id),
                    info = paste0("feature = ", f, "; exon ids: ", paste(unique(e$gene_id), collapse = ", ")))
    }
    e
}

test_that("read_gff_exons resolves Ensembl transcript: parents to the gene: id", {
    path <- write_gff(c(
        gff_line("1H", "gene", 100, 900, "ID=gene:HORVU.MOREX.r3.1HG0000010;biotype=protein_coding"),
        gff_line("1H", "mRNA", 100, 900,
                 "ID=transcript:HORVU.MOREX.r3.1HG0000010.1;Parent=gene:HORVU.MOREX.r3.1HG0000010"),
        gff_line("1H", "exon", 100, 300,
                 "Parent=transcript:HORVU.MOREX.r3.1HG0000010.1;Name=HORVU.MOREX.r3.1HG0000010.1.exon1"),
        gff_line("1H", "exon", 700, 900,
                 "Parent=transcript:HORVU.MOREX.r3.1HG0000010.1;Name=HORVU.MOREX.r3.1HG0000010.1.exon2")
    ))
    e <- exon_ids_join(path, c("gene", "mRNA"))
    expect_identical(e$gene_id, rep("gene:HORVU.MOREX.r3.1HG0000010", 2L))
})

test_that("read_gff_exons resolves RefSeq rna- parents to the gene- id", {
    path <- write_gff(c(
        gff_line("NC_1", "gene", 100, 900, "ID=gene-LOC123456;Name=LOC123456"),
        gff_line("NC_1", "mRNA", 100, 900, "ID=rna-XM_045123.1;Parent=gene-LOC123456"),
        gff_line("NC_1", "exon", 100, 300, "ID=exon-XM_045123.1-1;Parent=rna-XM_045123.1"),
        gff_line("NC_1", "CDS", 150, 300, "ID=cds-XP_044123.1;Parent=rna-XM_045123.1")
    ))
    e <- exon_ids_join(path, c("gene", "mRNA"))
    expect_identical(e$gene_id, "gene-LOC123456")
})

test_that("read_gff_exons gives a multi-parent exon one row per parent", {
    path <- write_gff(c(
        gff_line("1", "gene", 100, 900, "ID=A"),
        gff_line("1", "mRNA", 100, 900, "ID=A.1;Parent=A"),
        gff_line("1", "mRNA", 100, 700, "ID=A.2;Parent=A"),
        gff_line("1", "exon", 100, 300, "Parent=A.1,A.2")
    ))
    e <- exon_ids_join(path, c("gene", "mRNA"))
    # The ".N" strip on the whole Parent= value used to give "A.1,A", which joins nothing.
    expect_identical(e$gene_id, c("A", "A"))
    expect_identical(e$start, c(100L, 100L))
})

test_that("read_gff_exons keeps the ids of the TAIR, SIMDATA and Trifolium conventions", {
    tair <- write_gff(c(
        gff_line("1", "gene", 100, 900, "ID=AT1G01010"),
        gff_line("1", "mRNA", 100, 900, "ID=AT1G01010.1;Parent=AT1G01010"),
        gff_line("1", "exon", 100, 300, "Parent=AT1G01010.1")
    ))
    expect_identical(exon_ids_join(tair, c("gene", "mRNA"))$gene_id, "AT1G01010")

    simdata <- write_gff(c(   # SIMDATA.gff3: gene > mRNA > CDS, no exon
        gff_line("1", "gene", 100, 900, "ID=gene001;Name=DRY1"),
        gff_line("1", "mRNA", 100, 900, "ID=gene001.1;Parent=gene001;Name=DRY1-mRNA"),
        gff_line("1", "CDS", 120, 880, "ID=gene001.1.cds;Parent=gene001.1")
    ))
    expect_identical(exon_ids_join(simdata, c("gene", "mRNA"))$gene_id, "gene001")

    trifolium <- write_gff(c(   # Trifolium_genes_chr1-6.gff3: transcript (no Parent=) > exon
        gff_line("1", "transcript", 100, 900, "ID=g1.t1;em_desc=x"),
        gff_line("1", "exon", 100, 300, "Parent=g1.t1")
    ))
    expect_identical(exon_ids_join(trifolium, "transcript")$gene_id, "g1.t1")
})

# --- normalize_chr ---------------------------------------------------------

test_that("normalize_chr matches the pipeline's chromosome-name contract", {
    # The R twin of scripts/normalize_gff.py canonical(): what filter_vcf +
    # plink --output-chr MT emit for the same raw names (measured, plink
    # v1.90b7.2). The app's raw-GFF loader and the raw-VCF density table are
    # joined against pipeline tables carrying exactly these strings.
    expect_identical(normalize_chr(c("chr1", "Chr1", "CHR1", "1")), rep("1", 4))
    # plink parses its five codes case-insensitively and writes them back as
    # X / Y / XY / MT — a TAIR-style `Mt` therefore has to fold to `MT`.
    expect_identical(normalize_chr(c("X", "x", "chrX", "Y", "XY", "M", "MT", "Mt", "chrM")),
                     c("X", "X", "X", "Y", "XY", "MT", "MT", "MT", "MT"))
    # Everything else is kept verbatim, prefix aside (plink --allow-extra-chr).
    expect_identical(normalize_chr(c("2H", "chr2H", "Chr2H", "Pt", "Un", "scaffold_1")),
                     c("2H", "2H", "2H", "Pt", "Un", "scaffold_1"))
    # Idempotent — safe on an already-normalised name.
    x <- c("chrX", "Chr1", "Mt", "2H")
    expect_identical(normalize_chr(normalize_chr(x)), normalize_chr(x))
})
