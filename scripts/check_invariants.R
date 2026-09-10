#!/usr/bin/env Rscript
# Validate a pipeline results tree against the invariants in
# scripts/R/lib/invariants.R.
#
# Usage:
#   Rscript /pipeline/scripts/check_invariants.R <results_dir> [--modules GEA,GWAS] [--warn-only]
#
# Exits 1 when any violation has severity "error" (0 with --warn-only, or when
# only "warn" rows are present), so it can gate a merge the same way
# tests/run_tests.R does.
#
# Every actual CHECK lives in scripts/R/lib/invariants.R so it can be unit-tested
# against hand-built fixtures without a results tree — see
# tests/testthat/test-invariants.R. What stays here is path resolution, reading,
# and the data SHAPING that turns files into the checkers' arguments.
#
# That shaping is not nothing, and it is NOT covered by any test. Both schema
# misreads found while building this file lived here, not in the lib:
#   * selected_snps.tsv's per-method columns hold TRAIT NAMES, not p-values
#   * pairwise_overlap_table.tsv's source_a/source_b are "climate"/"phenotype",
#     not "GEA"/"GWAS"
# The reusable half has since moved to the lib (selected_snps_traits()); what
# remains — SOURCE_TO_MODULE, the glob exclusions, the numeric coercions — is
# where the next schema change will bite first. Check it against a real tree
# after any output-format change.
#
# NOTE ON EXPECTED FAILURES: this validator is EXPECTED to report violations on
# the current SIMDATA tree. They correspond to defects already filed in
# docs/pipeline_improvement_requests.md. A clean run on SIMDATA would itself be
# a signal that something stopped working. See the "predicted red" list in that
# file.

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
    stop("usage: check_invariants.R <results_dir> [--modules GEA,GWAS] [--warn-only]")
}
RESULTS_DIR <- normalizePath(args[[1L]], mustWork = TRUE)
WARN_ONLY   <- "--warn-only" %in% args
MODULES     <- {
    i <- match("--modules", args)
    if (is.na(i) || i >= length(args)) c("GEA", "GWAS")
    else strsplit(args[[i + 1L]], ",", fixed = TRUE)[[1]]
}

# The lib lives beside this script.
.here <- local({
    fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(fa) == 0) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})
source(file.path(.here, "R", "lib", "invariants.R"))

# ── reading ───────────────────────────────────────────────────────────────────

# Every id-bearing column is forced to character. A numeric-looking chromosome
# or SNPID inferred as integer breaks every %in% and join downstream — the trap
# CLAUDE.md warns about and tests/testthat/test-io.R pins.
ID_COLS <- c("chr", "SNPID", "site", "sample", "region_id", "gene_id", "trait",
             "snp_ids", "methods", "traits", "other_traits",
             "overlap_traits", "overlap_snps", "step", "metric", "value")

read_tsv <- function(path) {
    if (!file.exists(path) || file.size(path) == 0) return(NULL)
    hdr <- names(fread(path, nrows = 0L))
    cc  <- intersect(ID_COLS, hdr)
    fread(path, sep = "\t", colClasses = setNames(rep("character", length(cc)), cc),
          showProgress = FALSE)
}

mod_tbl <- function(module, ...) file.path(RESULTS_DIR, module, "tables", ...)

# ── run ───────────────────────────────────────────────────────────────────────

V <- list()
add <- function(x) if (!is.null(x) && nrow(x) > 0) V[[length(V) + 1L]] <<- x

message("INFO: checking ", RESULTS_DIR, " (modules: ", paste(MODULES, collapse = ", "), ")")

summary_dt <- read_tsv(file.path(RESULTS_DIR, "pipeline_summary.tsv"))
add(check_summary_accounting(summary_dt))
add(check_filtering_monotone(
    read_tsv(file.path(RESULTS_DIR, "Processing", "tables", "filtering_summary.tsv"))))

chr_tables    <- list()
sel_by_module <- list()
canonical_chr <- character(0)   # union of chr over every TESTED SNP (p-value tables)

for (mod in MODULES) {
    sel      <- read_tsv(mod_tbl(mod, "selected_snps.tsv"))
    per      <- read_tsv(mod_tbl(mod, "regions_per_trait.tsv"))
    comb     <- read_tsv(mod_tbl(mod, "regions_combined.tsv"))
    genes    <- read_tsv(mod_tbl(mod, "genes_per_region.tsv"))
    if (is.null(sel) && is.null(per) && is.null(comb)) next
    sel_by_module[[mod]] <- sel

    if (!is.null(sel)) {
        # SCHEMA NOTE, easy to get wrong: in selected_snps.tsv the per-method
        # columns hold TRAIT NAMES, not p-values —
        #     SNPID  chr  pos  EMMAX  LFMM        RDA  min_pvalue
        #     1:...  1    ...  bio_2  bio_2,bio_3 ""   6.2e-08
        # so min_pvalue is the minimum over p-values that live in
        # methods/{m}/{m}_pvalues_K*.tsv, not over anything in this table.
        # check_min_pvalue_column() is therefore NOT applicable here.
        # Empty string (written as the two-character token "") is the pipeline's
        # "this method had no hit" marker.
        sel_num <- copy(sel)
        sel_num[, min_pvalue := suppressWarnings(as.numeric(min_pvalue))]
        add(check_pvalues_table(sel_num, "min_pvalue", paste0(mod, "/selected_snps.tsv")))
        chr_tables[[paste0(mod, "/selected_snps.tsv")]] <- sel

        if (!is.null(summary_dt)) {
            step <- tolower(mod)
            exp  <- list(selected_snps_total = nrow(sel))
            if (!is.null(per))  exp$regions_per_trait <- nrow(per)
            if (!is.null(comb)) exp$regions_combined  <- nrow(comb)
            add(check_summary_counts(summary_dt, step, exp))
        }
    }

    for (nm in c("regions_per_trait.tsv", "regions_combined.tsv")) {
        tb <- if (nm == "regions_per_trait.tsv") per else comb
        if (is.null(tb)) next
        add(check_regions_table(tb, sel, paste0(mod, "/", nm)))
        add(check_column_names_unique(names(tb), paste0(mod, "/", nm)))
        chr_tables[[paste0(mod, "/", nm)]] <- tb
    }

    for (nm in c("genes_per_region.tsv", "genes_per_region_collapsed.tsv", "genes_combined.tsv")) {
        tb <- read_tsv(mod_tbl(mod, nm))
        if (is.null(tb)) next
        add(check_column_names_unique(names(tb), paste0(mod, "/", nm)))
        if (nm == "genes_per_region.tsv") {
            add(check_genes_table(tb, per, paste0(mod, "/", nm)))
            chr_tables[[paste0(mod, "/", nm)]] <- tb
        }
    }

    # The real p-value tables: one column per trait, every cell a p-value.
    for (f in Sys.glob(mod_tbl(mod, "methods", "*", "*_pvalues_K*.tsv"))) {
        if (grepl("_sig_snps_|_wza_", basename(f))) next
        tb <- read_tsv(f)
        if (is.null(tb)) next
        tcols <- setdiff(names(tb), c("SNPID", "chr", "pos", "n_snps", "mean_maf"))
        for (cn in tcols) tb[, (cn) := suppressWarnings(as.numeric(get(cn)))]
        add(check_pvalues_table(tb, tcols, file.path(mod, "methods", basename(f))))
        add(check_column_names_unique(names(tb), file.path(mod, "methods", basename(f))))
        # Every SNP that was TESTED. This is the canonical chromosome set: a
        # region or gene table naming a chromosome absent from it is a real
        # violation. Without it, check_chromosome_names() derives "canonical"
        # from the union of what it was given, which makes every table a subset
        # by construction and leaves only the ^chr-prefix half of the check live.
        if ("chr" %in% names(tb)) {
            canonical_chr <- union(canonical_chr, unique(as.character(tb$chr)))
        }
    }

    # Per-method significant-SNP tables: the overlap_traits / overlap_snps
    # annotation must name other traits and real SNPs.
    sig_pool     <- list()
    sig_variants <- list()
    for (f in Sys.glob(mod_tbl(mod, "methods", "*", "*_sig_snps_*.tsv"))) {
        tb <- read_tsv(f)
        if (is.null(tb)) next
        add(check_sig_snp_overlaps(tb, known_snpids = if (!is.null(sel)) sel$SNPID else NULL,
                                   table_name = file.path(mod, "methods", basename(f))))
        if (all(c("SNPID", "pvalue") %in% names(tb))) {
            sig_pool[[f]] <- tb[, .(SNPID = as.character(SNPID),
                                    pvalue = suppressWarnings(as.numeric(pvalue)))]
        }
        m <- basename(dirname(f))
        # "EMMAX_pvalues_K3_sig_snps_bonf_0.05.tsv" -> "bonf_0.05"
        tok <- sub("\\.tsv$", "", sub("^.*_sig_snps_", "", basename(f)))
        sig_variants[[m]] <- c(sig_variants[[m]], tok)
    }
    add(check_single_threshold_variant(sig_variants, paste0(mod, "/methods")))
    if (!is.null(sel) && length(sig_pool) > 0) {
        add(check_min_pvalue_against_sig_snps(sel, rbindlist(sig_pool, use.names = TRUE),
                                              paste0(mod, "/selected_snps.tsv")))
    }
    # DELIBERATELY NOT CHECKED: summary's sig_snps_{method} against the sig
    # tables. The summary counts the ONE threshold the config selected, while
    # the tables on disk may include stale variants from earlier configs — so
    # the comparison is ambiguous by construction and would report a config
    # artefact as an arithmetic error. check_single_threshold_variant() above
    # catches the underlying staleness directly instead.
    if (!is.null(summary_dt)) {
        genes_combined <- read_tsv(mod_tbl(mod, "genes_combined.tsv"))
        if (!is.null(genes_combined)) {
            add(check_summary_counts(summary_dt, tolower(mod),
                                     list(genes_found = nrow(genes_combined))))
        }
    }
}

# Cross-module referential integrity: GEAxGWAS is built from the GEA and GWAS
# selected-SNP tables and goes stale silently when either is re-run alone.
pairwise <- read_tsv(file.path(RESULTS_DIR, "GEAxGWAS", "tables", "pairwise_overlap_table.tsv"))
if (!is.null(pairwise) && all(c("trait_a", "source_a", "trait_b", "source_b") %in% names(pairwise))) {
    # A trait of this module is a value INSIDE the per-method columns of
    # selected_snps.tsv, not a column name — the column names are methods.
    traits_of <- function(mod) {
        sel <- sel_by_module[[mod]]
        if (is.null(sel)) return(character(0))
        mcols <- setdiff(names(sel), c("SNPID", "chr", "pos", "min_pvalue"))
        vals  <- unlist(lapply(mcols, function(cn) as.character(sel[[cn]])), use.names = FALSE)
        vals  <- gsub('"', "", vals)
        vals  <- unlist(strsplit(vals[nzchar(vals)], ",", fixed = TRUE), use.names = FALSE)
        unique(trimws(vals))
    }
    # source_a / source_b are "climate" (GEA) and "phenotype" (GWAS), not the
    # module names.
    SOURCE_TO_MODULE <- c(climate = "GEA", phenotype = "GWAS")
    long <- rbindlist(list(
        pairwise[, .(source = tolower(source_a), trait = trait_a)],
        pairwise[, .(source = tolower(source_b), trait = trait_b)]))
    long[, module := SOURCE_TO_MODULE[source]]
    for (mod in intersect(MODULES, c("GEA", "GWAS"))) {
        ds <- long[!is.na(module) & module == mod, unique(trait)]
        us <- traits_of(mod)
        if (length(ds) == 0 || length(us) == 0) next
        add(check_referential_integrity(
            ds, us,
            check_name = "pairwise_table_references_unknown_trait",
            table_name = "GEAxGWAS/pairwise_overlap_table.tsv",
            what = paste0(mod, " trait")))
    }
}

# Genetic offsets, one table per (method, run) cell.
meta <- read_tsv(file.path(RESULTS_DIR, "Processing", "tables", "metadata.tsv"))
for (f in Sys.glob(file.path(RESULTS_DIR, "Maladaptation", "tables", "*", "*",
                             "genetic_offset_site.tsv"))) {
    tb <- read_tsv(f)
    if (is.null(tb)) next
    label <- file.path("Maladaptation", basename(dirname(dirname(f))),
                       basename(dirname(f)), "genetic_offset_site.tsv")
    add(check_offsets_table(tb, "genetic_offset", label))
    # No offset may be reported for a site the run does not have.
    if (!is.null(meta) && "site" %in% names(meta) && "site" %in% names(tb)) {
        add(check_referential_integrity(
            tb$site, meta$site,
            check_name = "offset_reported_for_unknown_site",
            table_name = label, what = "site"))
    }
}

# The "normalize early" promise, across every chr-bearing table of the run.
add(check_chromosome_names(chr_tables,
                           canonical = if (length(canonical_chr) > 0) canonical_chr else NULL))

# ── report ────────────────────────────────────────────────────────────────────

viol <- if (length(V) == 0) no_violations() else rbindlist(V, use.names = TRUE)
n_err  <- sum(viol$severity == "error")
n_warn <- sum(viol$severity == "warn")

if (nrow(viol) == 0) {
    message("INFO: no violations.")
} else {
    setorder(viol, severity, check, table)
    cat("\n")
    for (i in seq_len(nrow(viol))) {
        cat(sprintf("[%-5s] %-45s %s\n           %s\n           %s\n",
                    viol$severity[i], viol$check[i], viol$table[i],
                    paste0("key: ", viol$key[i]), viol$detail[i]))
    }
    cat("\n")
}
message(sprintf("INFO: %d violation(s): %d error, %d warn", nrow(viol), n_err, n_warn))

if (n_err > 0 && !WARN_ONLY) quit(status = 1)
