#!/usr/bin/env Rscript
# promote_snp_set.R — Promote pipeline GEA/GWAS significant SNPs into a maladaptation SNP set.
#
# Run BETWEEN mode=gea (or mode=gwas) and mode=maladaptation, because
# resolve_active_snp_sets() in common.smk globs _intermediate/snp_sets/* at Snakemake
# parse time — the set directory must exist before `snakemake --config mode=maladaptation` runs.
#
# Usage:
#   Rscript promote_snp_set.R SOURCE_SELECTED_SNPS SET_NAME SNP_SETS_DIR [SOURCE_MODULE]
#
# Args:
#   1  SOURCE_SELECTED_SNPS  Path to GEA/GWAS combined selected_snps.tsv
#                            (e.g. {PROJECT}_results/GEA/tables/selected_snps.tsv)
#                            Schema: SNPID, chr, pos, <method cols...>, min_pvalue
#   2  SET_NAME              Name for the promoted set (alphanumeric + _ + .)
#                            Used as the snp_set wildcard in maladaptation rules.
#                            Must be listed in Maladaptation.snp_sets in the config.
#   3  SNP_SETS_DIR          Path to _intermediate/snp_sets/ directory
#                            (e.g. {PROJECT}_results/_intermediate/snp_sets/)
#   4  SOURCE_MODULE         (optional) "GEA" | "GWAS" — recorded in manifest. Default "pipeline"
#
# Output:
#   {SNP_SETS_DIR}/{SET_NAME}/selected_snps.tsv  — canonical SNP set (SNPID, chr, pos, min_pvalue)
#   {SNP_SETS_DIR}/manifest.json                  — upserted with this set's entry (Shiny-compatible)
#
# Manifest contract (shared with scripts/clinego.app/R/fct_snp_sets.R):
#   {SNP_SETS_DIR}/manifest.json  — JSON array of objects:
#     { name, n_snps, created, source_module }
#   Written by this script OR by Shiny save_snp_set(); pipeline resolves sets by glob only.

suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: promote_snp_set.R SOURCE_SELECTED_SNPS SET_NAME SNP_SETS_DIR [SOURCE_MODULE]")
}

SOURCE_FILE    <- args[1]
SET_NAME       <- args[2]
SNP_SETS_DIR   <- args[3]
SOURCE_MODULE  <- if (length(args) >= 4) args[4] else "pipeline"

#=============================================================================
# Validate set name
#=============================================================================
if (!grepl("^[A-Za-z0-9_.]+$", SET_NAME)) {
    stop("ERROR: SET_NAME must match [A-Za-z0-9_.]+, got: ", SET_NAME)
}
message("INFO: Promoting SNP set '", SET_NAME, "' from ", SOURCE_FILE)

#=============================================================================
# Read source selected_snps.tsv
#=============================================================================
if (!file.exists(SOURCE_FILE)) {
    stop("ERROR: Source file not found: ", SOURCE_FILE)
}

# Read with correct column types (mirrors io_selected_snps.R::read_selected_snps)
dt <- fread(SOURCE_FILE, sep = '\t', header = TRUE,
            colClasses = c("SNPID" = "character", "chr" = "character"))
message("INFO: Source table loaded — ", nrow(dt), " rows, columns: ",
        paste(names(dt), collapse = ", "))

if (nrow(dt) == 0) {
    stop("ERROR: Source selected_snps.tsv is empty — no significant SNPs to promote")
}

# Validate SNPID format (expected chr:pos)
bad_ids <- !grepl("^[^:]+:[0-9]+$", dt$SNPID)
if (any(bad_ids)) {
    message("WARNING: ", sum(bad_ids), " SNPIDs don't match expected chr:pos format. First: ",
            dt$SNPID[which(bad_ids)[1]])
}

#=============================================================================
# Build canonical SNP set: unique by SNPID, keep SNPID + chr + pos + min_pvalue
#=============================================================================
required_cols <- c("SNPID", "chr", "pos", "min_pvalue")
missing_req <- setdiff(required_cols, names(dt))
if (length(missing_req) > 0) {
    stop("ERROR: Source file missing required columns: ", paste(missing_req, collapse = ", "))
}

# Collapse to unique SNPIDs, keeping minimum pvalue across any method/trait columns
uniq <- dt[, .(chr      = chr[1L],
               pos      = pos[1L],
               min_pvalue = min(min_pvalue, na.rm = TRUE)),
           by = SNPID]
data.table::setcolorder(uniq, c("SNPID", "chr", "pos", "min_pvalue"))
data.table::setorder(uniq, chr, pos)
message("INFO: Unique SNPs: ", nrow(uniq))

#=============================================================================
# Write SNP set TSV (atomic write: tempfile → rename)
#=============================================================================
set_dir <- file.path(SNP_SETS_DIR, SET_NAME)
dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)

tsv_path <- file.path(set_dir, "selected_snps.tsv")
tmp_tsv  <- tempfile(tmpdir = set_dir, fileext = ".tsv")
fwrite(uniq, tmp_tsv, sep = '\t', quote = FALSE)
file.rename(tmp_tsv, tsv_path)
message("INFO: SNP set written → ", tsv_path)

#=============================================================================
# Upsert manifest.json (atomic write: tempfile → rename)
# Manifest schema: JSON array of { name, n_snps, created, source_module, ... }
# Mirrors .write_manifest_atomic() in scripts/clinego.app/R/fct_snp_sets.R
#=============================================================================
manifest_path <- file.path(SNP_SETS_DIR, "manifest.json")
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)

# Read existing manifest (empty list if absent/corrupt)
man <- tryCatch(
    jsonlite::read_json(manifest_path, simplifyVector = FALSE),
    error = function(e) list()
)

# Remove existing entry for this set name (upsert semantics)
man <- Filter(function(x) !identical(x$name, SET_NAME), man)

# Build new entry
new_entry <- list(
    name          = SET_NAME,
    n_snps        = nrow(uniq),
    created       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    source_module = SOURCE_MODULE
)
man <- c(man, list(new_entry))

# Atomic write
tmp_json <- tempfile(tmpdir = dirname(manifest_path), fileext = ".json")
tryCatch({
    jsonlite::write_json(man, tmp_json, auto_unbox = TRUE, pretty = TRUE)
    file.rename(tmp_json, manifest_path)
    message("INFO: Manifest upserted → ", manifest_path,
            " (", length(man), " set(s) total)")
}, error = function(e) {
    if (file.exists(tmp_json)) unlink(tmp_json)
    message("WARNING: Could not write manifest: ", e$message)
})

message("INFO: promote_snp_set.R complete — set '", SET_NAME, "' ready for mode=maladaptation")
message("INFO: Add 'Maladaptation.snp_sets: [", SET_NAME, "]' (or 'all') to your config")
