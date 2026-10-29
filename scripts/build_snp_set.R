#!/usr/bin/env Rscript
# build_snp_set.R — build one SNP set declared in config (GEA.snp_sets /
# GWAS.snp_sets) into the curated SNP set store mode=maladaptation reads.
#
# The selection is the Shiny GEA tab's own (scripts/R/lib/snp_set_select.R):
# each declared trait::method cell is thresholded with its resolved rule
# (inclusive <=), the cells are combined with the declared strategy, and the
# result is collapsed to SNPID / chr / pos / min_pvalue. A config block that
# reproduces a set saved from the app therefore reproduces its TSV.
#
# Usage:
#   build_snp_set.R DECL OUT_TSV OUT_RECORD K_BEST PVAL_FILE [PVAL_FILE ...]
#
#   DECL        declaration.json written by the snp_set_declare rule — or a
#               snp_set.yaml record (a superset), which is how a set saved from
#               Shiny is rebuilt for a round-trip check
#   OUT_TSV     <store>/<set>/selected_snps.tsv
#   OUT_RECORD  <store>/<set>/snp_set.yaml (declaration + k_best, n_snps, created)
#   K_BEST      the K the p-value tables were computed under (recorded only)
#   PVAL_FILE   one wide p-value table per method; the method is the name of the
#               file's parent directory (tables/methods/<method>/...), as in
#               combine_selected_snps.R
#
# Zero selected SNPs is NOT an error: a header-only TSV and a record with
# n_snps: 0 are written and the script exits 0, so a strict declaration cannot
# fail mode=gea after every association job has run. mode=maladaptation skips a
# set of fewer than 3 SNPs (MIN_OFFSET_SNPS, workflow/rules/common.smk) under
# snp_sets: all and refuses it when named.

suppressPackageStartupMessages(library(data.table))

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/lib/snp_set_select.R")
source("/pipeline/scripts/R/lib/region_distance.R")   # an auto_genome_wide gap

args <- commandArgs(trailingOnly = TRUE)
################
DECL_PATH  <- args[1]
OUT_TSV    <- args[2]
OUT_RECORD <- args[3]
K_BEST     <- suppressWarnings(as.integer(args[4]))
PVAL_FILES <- args[-(1:4)]
################
if (length(args) < 5L || any(is.na(args[1:4])) || any(!nzchar(args[1:4])) || is.na(K_BEST)) {
    stop("Usage: build_snp_set.R <declaration.json|snp_set.yaml> <out_tsv> <out_record> ",
         "<k_best> <pvalue_file> [<pvalue_file> ...]\n  got ", length(args),
         " argument(s): ", paste(args, collapse = " | "))
}

clash <- intersect(normalizePath(c(DECL_PATH, PVAL_FILES), mustWork = FALSE),
                   normalizePath(c(OUT_TSV, OUT_RECORD), mustWork = FALSE))
if (length(clash) > 0) stop("Refusing to overwrite an input file: ", paste(clash, collapse = ", "))

decl <- read_snp_set_decl(DECL_PATH)
if (!identical(decl$regime, "snp")) {
    stop("SNP set '", decl$name, "': regime '", decl$regime,
         "' cannot be built from config (only 'snp')")
}
message("INFO: Building SNP set ", decl$name, " (", decl$source_module, ", origin ",
        decl$origin, ")")
gap_spec <- decl$clumping_distance
if (is.character(gap_spec)) {
    decl$clumping_distance <- resolve_clumping_distance_strict(gap_spec, decl$ld_decay_table %||% "NULL",
                                                               decl$clumping_ld_rule %||% "critical",
                                                               decl$ld_decay_group %||% "All")
    decl$clumping_distance_spec <- gap_spec
}
message("INFO: Strategy: ", decl$strategy, "; clumping distance: ", decl$clumping_distance,
        " bp", if (is.character(gap_spec)) paste0(" (", gap_spec, ", ", decl$clumping_ld_rule, ")") else "",
        "; ", length(decl$pairs), " trait::method pair(s)")

methods_of_files <- basename(dirname(PVAL_FILES))
if (anyDuplicated(methods_of_files)) {
    stop("Two p-value files for one method: ",
         paste(PVAL_FILES[methods_of_files %in% methods_of_files[duplicated(methods_of_files)]],
               collapse = ", "))
}
needed <- unique(split_snp_set_pairs(decl$pairs)$method)
absent <- setdiff(needed, methods_of_files)
if (length(absent) > 0) {
    stop("SNP set '", decl$name, "' names method(s) with no p-value file: ",
         paste(absent, collapse = ", "))
}

pvals <- lapply(PVAL_FILES[methods_of_files %in% needed], function(f) {
    message("INFO: Reading ", f)
    fread(f, sep = "\t", header = TRUE,
          colClasses = c(chr = "character", SNPID = "character"))
})
names(pvals) <- methods_of_files[methods_of_files %in% needed]

sel <- select_snp_set(pvals, decl$cells, decl$pairs,
                      strategy = decl$strategy, gap = decl$clumping_distance)

message("INFO: Per-cell rules and hits:")
for (i in seq_len(nrow(sel$cells))) {
    r <- sel$cells[i, ]
    message(sprintf("  %-40s %-6s %-10s cutoff=%-12s status=%-14s n=%d",
                    r$pair, r$type, format(r$value), format(r$threshold),
                    r$status, r$n_sig))
}

n <- nrow(sel$snps)
if (n == 0L) {
    message("WARNING: SNP set '", decl$name, "' selected 0 SNPs — writing an empty set. ",
            "mode=maladaptation skips it under snp_sets: all.")
} else if (n < 3L) {
    message("WARNING: SNP set '", decl$name, "' selected ", n, " SNP(s). ",
            "mode=maladaptation skips a set of fewer than 3 SNPs under snp_sets: all.")
}

dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)
tmp_tsv <- tempfile(tmpdir = dirname(OUT_TSV), fileext = ".tsv")
fwrite(sel$snps, tmp_tsv, sep = "\t", quote = FALSE)
if (!file.rename(tmp_tsv, OUT_TSV)) {
    unlink(tmp_tsv)
    stop("could not move the SNP set into place: ", OUT_TSV)
}

record <- decl
record$k_best  <- K_BEST
record$n_snps  <- n
record$created <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
write_snp_set_record(record, OUT_RECORD)

message("INFO: Wrote ", n, " SNPs to ", OUT_TSV, " and the record to ", OUT_RECORD)
