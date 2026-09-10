# fct_snp_sets.R
# Functions for managing curated SNP sets (saved from GEA tab for Maladaptation).
#
# Shared contract with the Snakemake pipeline:
#   Set store:  {PROJECT}_results/_intermediate/snp_sets/{name}/selected_snps.tsv
#               header'd TSV, columns: SNPID, chr, pos, min_pvalue; unique by SNPID
#   Manifest:   {PROJECT}_results/_intermediate/snp_sets/manifest.json
#               array of {name, n_snps, source_module, threshold_type, threshold_value,
#               regime, strategy, clumping_distance, traits, methods, created}
#               Written by Shiny ONLY — pipeline resolves sets by glob.
#   GF suffixes emitted by pipeline: {name}_spatial | {name}_nospatial

# ─── Manifest read ────────────────────────────────────────────────────────────

#' Read the SNP sets manifest (returns list of objects; empty list on missing/corrupt)
#' @noRd
read_snp_sets_manifest <- function(project) {
    if (is.null(project)) return(list())
    path <- snp_sets_manifest_path(project)
    if (!file_ok(path)) return(list())
    tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(e) list()
    )
}

# ─── List for picker ──────────────────────────────────────────────────────────

#' Return data.frame(name, n_snps) in manifest order; 0-row if no sets
#' @noRd
list_snp_sets <- function(project) {
    m <- read_snp_sets_manifest(project)
    if (length(m) == 0) return(data.frame(name = character(0), n_snps = integer(0),
                                           stringsAsFactors = FALSE))
    data.frame(
        name   = vapply(m, function(x) x$name   %||% "", character(1)),
        n_snps = vapply(m, function(x) as.integer(x$n_snps %||% 0L), integer(1)),
        stringsAsFactors = FALSE
    )
}

#' Return TRUE if a set with this name already exists
#' @noRd
set_exists <- function(project, name) {
    isTRUE(name %in% list_snp_sets(project)$name)
}

# ─── Save ─────────────────────────────────────────────────────────────────────

#' Save a curated SNP set (writes TSV + atomically upserts manifest).
#'
#' @param project  Project name string.
#' @param name     Set name — must match [A-Za-z0-9_.]+.
#' @param sigsnps_dt  interactive_sigsnps() result (8-col long data.table).
#' @param params_list List of manifest metadata fields (see contract above).
#' @return Invisible integer: number of unique SNPs written.
#' @noRd
save_snp_set <- function(project, name, sigsnps_dt, params_list) {
    dir <- file.path(snp_sets_dir(project), name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    # Collapse long dt → UNIQUE by SNPID, keep chr/pos, min pvalue across methods/traits
    dt <- data.table::as.data.table(sigsnps_dt)
    uniq <- dt[, .(chr = chr[1L], pos = pos[1L],
                   min_pvalue = min(min_pvalue, na.rm = TRUE)), by = SNPID]
    data.table::setcolorder(uniq, c("SNPID", "chr", "pos", "min_pvalue"))
    data.table::setorder(uniq, chr, pos)

    # Atomic TSV write
    tsv     <- snp_set_path(project, name)
    tmp_tsv <- tempfile(tmpdir = dir, fileext = ".tsv")
    data.table::fwrite(uniq, tmp_tsv, sep = "\t", quote = FALSE)
    file.rename(tmp_tsv, tsv)

    # Atomic manifest upsert (replace same-name entry or append)
    n     <- nrow(uniq)
    entry <- c(list(
        name          = name,
        n_snps        = n,
        created       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    ), params_list)

    man  <- read_snp_sets_manifest(project)
    keep <- Filter(function(x) !identical(x$name, name), man)
    .write_manifest_atomic(project, c(keep, list(entry)))

    invisible(n)
}

# ─── Delete ───────────────────────────────────────────────────────────────────

#' Delete a curated SNP set and (optionally) its GF pipeline results.
#'
#' GF result directories use suffixes {name}_spatial and {name}_nospatial
#' (both variants from spatial_correction: both). We match EXACTLY these two
#' suffixes rather than a raw glob prefix to avoid over-deleting a set named
#' "foo" when "foo_bar" also exists.
#' @noRd
delete_snp_set <- function(project, name, remove_gf_results = TRUE) {
    unlink(file.path(snp_sets_dir(project), name), recursive = TRUE, force = TRUE)

    man  <- read_snp_sets_manifest(project)
    man  <- Filter(function(x) !identical(x$name, name), man)
    .write_manifest_atomic(project, man)

    if (isTRUE(remove_gf_results)) {
        # Known suffixes the pipeline emits for this set name
        suffixes <- c(name,
                      paste0(name, "_spatial"),
                      paste0(name, "_nospatial"))
        # Remove results for all registered maladaptation methods
        all_methods <- c("gradient_forest", "geometric_offset", "rda_offset")
        for (method in all_methods) {
            for (base in c(
                mod_path(project, MOD_MALAD, "plots",  method),
                mod_path(project, MOD_MALAD, "tables", method),
                mod_path(project, MOD_INTER, method)
            )) {
                for (suf in suffixes) {
                    target <- file.path(base, suf)
                    if (file.exists(target) || dir.exists(target))
                        unlink(target, recursive = TRUE, force = TRUE)
                }
            }
        }
    }

    invisible(TRUE)
}

# ─── Private helpers ──────────────────────────────────────────────────────────

#' Atomically write the manifest JSON (tempfile + rename idiom from fct_region_params.R)
#' @noRd
.write_manifest_atomic <- function(project, man) {
    path <- snp_sets_manifest_path(project)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".json")
    tryCatch({
        jsonlite::write_json(man, tmp, auto_unbox = TRUE, pretty = TRUE)
        file.rename(tmp, path)
    }, error = function(e) {
        if (file.exists(tmp)) unlink(tmp)
        message("WARNING: Could not write snp_sets manifest: ", e$message)
    })
}
