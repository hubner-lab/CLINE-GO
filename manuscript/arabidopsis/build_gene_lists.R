#!/usr/bin/env Rscript
# build_gene_lists.R: the known-gene lists the Arabidopsis1135 regions are scored against.
#
# The rules below were fixed, and committed, before any Arabidopsis1135 region table existed, so no
# list can be fitted to the results. Each list is independent of our scans.
#
# Reads data/arabidopsis1135/sources/ and data/Arabidopsis_1001G.gff3. Writes
# data/arabidopsis1135/lists/, one table per list with the columns item, label, chr, start, end:
#   flor_id.tsv             FLOR-ID's flowering-time genes (archived site, 306 genes)
#   grenenet_blocks.tsv     GrENE-net parallel-selection blocks in at least 15 gardens
#   grenenet_block_genes.tsv  the genes in those blocks
#   aragwas_ft.tsv          AraGWAS FT10 and FT16 hits, label = trait
#   classic_genes.tsv       FRI, FLC, FT, SVP, DOG1
#   abiotic.tsv             experimental GO annotations to four abiotic-stress responses
#   lists_summary.tsv, provenance.tsv
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -v $PWD:/data -w /data \
#     cline-go:6cdedcc Rscript manuscript/arabidopsis/build_gene_lists.R

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(jsonlite)
  library(readxl)
  library(rvest)
})

# The rules ------------------------------------------------------------------------------------
flor_id_genes <- 306L                  # the count Bouche et al. 2016 (NAR) state
blocks_sheet <- "num_significant_garden_15"
min_gardens <- 15L                     # blocks significant in at least this many gardens
n_blocks_expected <- 22L
n_block_genes_expected <- 837L         # the design's 836 lost the first gene to a header row
classic_genes <- c(FRI = "AT4G00650", FLC = "AT5G10140", FT = "AT1G65480",
                   SVP = "AT2G22540", DOG1 = "AT5G45830")
# FT10: 1001G panel (261) and the earlier panel of about 190 lines (31). FT16: 262 and 30.
aragwas_studies <- c("30" = "FT16", "31" = "FT10", "261" = "FT10", "262" = "FT16")
aragwas_min_maf <- 0.05                # reproduces every study's own hit counts, and it is our MAF
abiotic_terms <- c("GO:0009409" = "response to cold",
                   "GO:0009408" = "response to heat",
                   "GO:0009414" = "response to water deprivation",
                   "GO:0009651" = "response to salt stress")
experimental_evidence <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")

source_dir <- here("data", "arabidopsis1135", "sources")
list_dir <- here("data", "arabidopsis1135", "lists")
gff_path <- here("data", "Arabidopsis_1001G.gff3")
dir.create(list_dir, showWarnings = FALSE)

write_tsv <- function(x, name) fwrite(x, file.path(list_dir, name), sep = "\t", na = "")

# Genes from the project's own GFF ---------------------------------------------------------------
gff <- fread(cmd = sprintf("awk -F'\\t' '$3 == \"gene\" || $3 == \"ncRNA_gene\"' %s", shQuote(gff_path)),
             sep = "\t", header = FALSE, quote = "",
             col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "phase",
                           "attributes"),
             colClasses = c(V1 = "character"))
gff[, gene_id := sub(".*gene_id=([^;]+).*", "\\1", attributes)]
gff[, symbol := fifelse(grepl("(^|;)Name=", attributes),
                        sub(".*(^|;)Name=([^;]+).*", "\\2", attributes), "")]
genes <- gff[, .(gene_id, symbol, feature, chr, start, end)]
stopifnot(!anyDuplicated(genes$gene_id))

gene_rows <- function(ids, labels = NULL) {
  rows <- data.table(item = ids, label = if (is.null(labels)) "" else labels)
  rows <- genes[rows, on = c(gene_id = "item")]
  missing <- rows[is.na(chr), gene_id]
  if (length(missing)) message("Not in the GFF: ", paste(missing, collapse = ", "))
  rows[!is.na(chr), .(item = gene_id, label, chr, start, end)]
}

# FLOR-ID ---------------------------------------------------------------------------------------
flor_rows <- html_elements(read_html(file.path(source_dir, "florid_archive",
                                               "gene_list_flowering.html")), "tbody tr")
flor_id <- rbindlist(lapply(flor_rows, function(row) {
  agi <- unique(toupper(regmatches(as.character(row),
                                   gregexpr("A[Tt][1-5CMcm][Gg][0-9]{5}", as.character(row)))[[1]]))
  cells <- html_text2(html_elements(row, xpath = "./td"))
  data.table(agi = agi, name = cells[2])
}))
flor_id <- unique(flor_id, by = "agi")
message(sprintf("FLOR-ID archived list: %d rows, %d genes", length(flor_rows), nrow(flor_id)))
stopifnot(nrow(flor_id) == flor_id_genes)
flor_table <- gene_rows(flor_id$agi, flor_id$name)

# GrENE-net parallel blocks --------------------------------------------------------------------
# Each table has one sheet per threshold. Table S14 has no header row: 837 genes at 15 gardens.
tables_dir <- file.path(source_dir, "dl_grene", "Supplemental_tables")
blocks <- as.data.table(read_excel(
  file.path(tables_dir, "Table_S13_significant_parallel_blocks_aross_gardens.xlsx"),
  sheet = blocks_sheet
))
setnames(blocks, c("chr", "start", "end", "n_gardens"))
blocks[, chr := as.character(chr)]
stopifnot(nrow(blocks) == n_blocks_expected, all(blocks$n_gardens == min_gardens))
block_table <- blocks[, .(item = sprintf("block_%s_%.0f", chr, start), label = n_gardens,
                          chr, start, end)]

s14_genes <- as.data.table(read_excel(
  file.path(tables_dir, "Table_S14_genes_overlapping_significant_parallel_blocks.xlsx"),
  sheet = blocks_sheet, col_names = FALSE
))[[1]]
stopifnot(length(s14_genes) == n_block_genes_expected, !anyDuplicated(s14_genes))
block_gene_table <- gene_rows(s14_genes)

setkey(blocks, chr, start, end)
genes_in_blocks <- unique(foverlaps(genes, blocks, nomatch = NULL,
                                    by.x = c("chr", "start", "end")), by = "gene_id")$gene_id
message(sprintf("GrENE-net: %d blocks in >= %d gardens, %d genes in Table S14, %d of them inside the blocks on our GFF, %d GFF genes inside the blocks",
                nrow(blocks), min_gardens, length(s14_genes),
                sum(s14_genes %in% genes_in_blocks), length(genes_in_blocks)))

# AraGWAS ---------------------------------------------------------------------------------------
read_study <- function(pk, trait) {
  meta <- fromJSON(file.path(source_dir, "aragwas_ft", sprintf("study_%s_meta.json", pk)))
  scan <- fromJSON(file.path(source_dir, "aragwas_ft", sprintf("study_%s.json", pk)))
  hits <- rbindlist(lapply(1:5, function(i) {
    part <- scan[[paste0("chr", i)]]
    data.table(chr = as.character(i), pos = part$positions, score = part$scores, maf = part$mafs)
  }))
  hits <- hits[maf >= aragwas_min_maf &
                 score >= pmin(meta$permutationThreshold, meta$bonferroniThreshold)]
  hits[, `:=`(study = pk, trait = trait, n_samples = meta$numberSamples,
              above_permutation = score >= meta$permutationThreshold,
              above_bonferroni = score >= meta$bonferroniThreshold)]
  message(sprintf("AraGWAS %s (%s, n = %d): %d above permutation, %d above Bonferroni",
                  pk, trait, meta$numberSamples, sum(hits$above_permutation),
                  sum(hits$above_bonferroni)))
  stopifnot(sum(hits$above_permutation) == meta$nHitsPerm,
            sum(hits$above_bonferroni) == meta$nHitsBonf)
  hits
}
aragwas_hits <- rbindlist(Map(read_study, names(aragwas_studies), aragwas_studies))
aragwas_table <- aragwas_hits[, .(
  item = sprintf("%s_%s_%.0f", study, chr, pos), label = trait, chr, start = pos, end = pos,
  study, n_samples, score, maf, above_permutation, above_bonferroni
)]

# Classic genes ---------------------------------------------------------------------------------
classic_table <- gene_rows(unname(classic_genes), names(classic_genes))
stopifnot(nrow(classic_table) == length(classic_genes))

# Abiotic stress --------------------------------------------------------------------------------
gaf_path <- file.path(source_dir, "go", "tair.gaf.gz")
gaf <- fread(cmd = sprintf("zcat %s | awk -F'\\t' '!/^!/'", shQuote(gaf_path)), sep = "\t",
             header = FALSE, quote = "", fill = TRUE, select = c(2, 3, 4, 5, 7, 9),
             col.names = c("agi", "symbol", "qualifier", "term", "evidence", "aspect"))
abiotic <- gaf[term %in% names(abiotic_terms) & evidence %in% experimental_evidence &
                 aspect == "P" & !grepl("NOT", qualifier, fixed = TRUE) &
                 grepl("^AT[1-5CM]G[0-9]{5}$", agi)]
abiotic <- abiotic[, .(terms = paste(sort(unique(abiotic_terms[term])), collapse = "; ")), by = agi]
message(sprintf("Abiotic stress: %d genes from %d annotations of %d terms", nrow(abiotic),
                gaf[term %in% names(abiotic_terms) & evidence %in% experimental_evidence, .N],
                length(abiotic_terms)))
abiotic_table <- gene_rows(abiotic$agi, abiotic$terms)

# Write -----------------------------------------------------------------------------------------
tables <- list(flor_id = flor_table, grenenet_blocks = block_table,
               grenenet_block_genes = block_gene_table, aragwas_ft = aragwas_table,
               classic_genes = classic_table, abiotic = abiotic_table)
for (name in names(tables)) write_tsv(tables[[name]], paste0(name, ".tsv"))

write_tsv(data.table(list = names(tables), n_items = vapply(tables, nrow, integer(1)),
                     n_chromosomes = vapply(tables, function(x) uniqueN(x$chr), integer(1))),
          "lists_summary.tsv")

source_files <- c(file.path(source_dir, "florid_archive", "gene_list_flowering.html"),
                  file.path(tables_dir, "Table_S13_significant_parallel_blocks_aross_gardens.xlsx"),
                  file.path(tables_dir, "Table_S14_genes_overlapping_significant_parallel_blocks.xlsx"),
                  gaf_path,
                  file.path(source_dir, "aragwas_ft", sprintf("study_%s.json", names(aragwas_studies))),
                  gff_path)
write_tsv(data.table(file = sub(paste0(here(), "/"), "", source_files, fixed = TRUE),
                     sha256 = unname(tools::sha256sum(source_files)),
                     modified = format(file.mtime(source_files), "%Y-%m-%d")),
          "provenance.tsv")

print(fread(file.path(list_dir, "lists_summary.tsv")))
