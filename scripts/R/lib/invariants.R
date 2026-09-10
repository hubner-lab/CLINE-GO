# Data-independent invariants over a pipeline results tree.
#
# These are the checks that must hold for ANY dataset, as opposed to golden
# files, which pin one dataset's numbers. They catch a class golden files
# structurally cannot: two output files that are each internally consistent but
# disagree with each other, because one mode was re-run and the other was not.
#
# CONTRACT, uniform across every checker here:
#   * Input is data.tables, never paths. File reading belongs to the CLI wrapper
#     (scripts/check_invariants.R) so every checker stays unit-testable.
#   * Output is ALWAYS a violations data.table with the columns
#         check, severity, table, key, detail
#     — zero rows when the invariant holds. Never NULL, never a bare logical.
#   * A checker NEVER stop()s on a violation. A checker that throws can only
#     report its first finding, and the point is to report all of them.
#     stop() is reserved for a caller passing something structurally unusable.
#   * severity is "error" (fails the run) or "warn" (reported, does not fail).
#
# This file has no dependency beyond data.table.

# ── violation plumbing ────────────────────────────────────────────────────────

# The empty violations table. Also the schema contract every checker returns.
no_violations <- function() {
    # as.data.table(list(...)), NOT data.table(...): `key` is a reserved argument
    # of data.table(), so a column of that name is silently consumed as the key
    # spec and the result errors with "some columns are not in the data.table".
    data.table::as.data.table(list(
        check    = character(),
        severity = character(),
        table    = character(),
        key      = character(),
        detail   = character()
    ))
}

# Build violation rows. `key` and `detail` recycle against each other so a
# checker can emit one row per offending record in a single call.
violation <- function(check, severity, table, key, detail) {
    # `key` alone decides how many violations there are. detail is recycled to
    # match. Deliberately NOT max(length(key), length(detail)): callers build
    # detail with paste0(), and paste0("x ", character(0)) returns "x " — length
    # ONE, not zero — so a max() would invent a bogus violation every time a
    # checker found nothing.
    n <- length(key)
    if (n == 0L) return(no_violations())
    if (length(detail) == 0L) detail <- ""
    data.table::as.data.table(list(   # see no_violations() on why not data.table()
        check    = rep_len(as.character(check),    n),
        severity = rep_len(as.character(severity), n),
        table    = rep_len(as.character(table),    n),
        key      = rep_len(as.character(key),      n),
        detail   = rep_len(as.character(detail),   n)
    ))
}

# Concatenate checker outputs. Accepts NULLs so callers can build lists freely.
combine_violations <- function(...) {
    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, list(...))
    if (length(parts) == 0) return(no_violations())
    data.table::rbindlist(parts, use.names = TRUE)
}

# ── structural: regions ───────────────────────────────────────────────────────

# regions: region_id, chr, start, end, length, snp_count, snp_ids [, trait]
# snps:    optional selected-SNPs table (SNPID, chr, pos) to cross-check membership
check_regions_table <- function(regions, snps = NULL, table_name = "regions") {
    if (is.null(regions) || nrow(regions) == 0) return(no_violations())
    required <- c("region_id", "chr", "start", "end")
    missing  <- setdiff(required, names(regions))
    if (length(missing) > 0) {
        return(violation("regions_schema", "error", table_name, table_name,
                         paste0("missing required column(s): ", paste(missing, collapse = ", "))))
    }

    r <- data.table::copy(regions)
    r[, start := as.integer(start)][, end := as.integer(end)]
    out <- list()

    bad <- r[start > end]
    out$bounds <- violation("region_start_after_end", "error", table_name,
                            bad$region_id, paste0("start ", bad$start, " > end ", bad$end))

    if ("length" %in% names(r)) {
        # end - start, NOT end - start + 1: regions.R:90 and fct_overlap.R:185 agree.
        bad <- r[as.integer(length) != (end - start)]
        out$length <- violation("region_length_mismatch", "error", table_name,
                                bad$region_id,
                                paste0("length ", bad$length, " != end - start = ", bad$end - bad$start))
    }

    if ("snp_count" %in% names(r)) {
        bad <- r[as.integer(snp_count) < 1L]
        out$empty <- violation("region_has_no_snps", "error", table_name,
                               bad$region_id, paste0("snp_count = ", bad$snp_count))
    }

    if (all(c("snp_count", "snp_ids") %in% names(r))) {
        listed <- vapply(r$snp_ids, function(s) {
            if (is.na(s) || !nzchar(s)) return(0L)
            length(unique(strsplit(s, ",", fixed = TRUE)[[1]]))
        }, integer(1L), USE.NAMES = FALSE)
        bad_i <- which(listed != as.integer(r$snp_count))
        out$count <- violation("region_snp_count_disagrees_with_snp_ids", "error", table_name,
                               r$region_id[bad_i],
                               paste0("snp_count = ", r$snp_count[bad_i],
                                      " but snp_ids lists ", listed[bad_i]))
    }

    dup <- r[, .N, by = "region_id"][N > 1L]
    out$dup <- violation("region_id_duplicated", "error", table_name,
                         dup$region_id, paste0("appears ", dup$N, " times"))

    # Regions must not overlap within a group. Per-trait tables are grouped by
    # trait; a combined table is one group.
    grp_col <- if ("trait" %in% names(r)) "trait" else NULL
    r[, .grp := if (is.null(grp_col)) "" else as.character(get(grp_col))]
    data.table::setorder(r, .grp, chr, start)
    r[, .prev_end := data.table::shift(end),   by = c(".grp", "chr")]
    r[, .prev_id  := data.table::shift(region_id), by = c(".grp", "chr")]
    bad <- r[!is.na(.prev_end) & start <= .prev_end]
    out$overlap <- violation("regions_overlap_within_group", "error", table_name,
                             bad$region_id,
                             paste0("starts at ", bad$start, ", but ", bad$.prev_id,
                                    " on chr ", bad$chr, " ends at ", bad$.prev_end,
                                    if (!is.null(grp_col)) paste0(" (", grp_col, " = ", bad$.grp, ")") else ""))

    if (!is.null(snps) && nrow(snps) > 0 && "snp_ids" %in% names(r)) {
        known <- unique(as.character(snps$SNPID))
        pos_by_id <- stats::setNames(as.integer(snps$pos), as.character(snps$SNPID))
        chr_by_id <- stats::setNames(as.character(snps$chr), as.character(snps$SNPID))

        miss_key <- character(0); miss_det <- character(0)
        oob_key  <- character(0); oob_det  <- character(0)
        for (i in seq_len(nrow(r))) {
            s <- r$snp_ids[i]
            if (is.na(s) || !nzchar(s)) next
            ids <- unique(strsplit(s, ",", fixed = TRUE)[[1]])
            unknown <- setdiff(ids, known)
            if (length(unknown) > 0) {
                miss_key <- c(miss_key, rep(r$region_id[i], length(unknown)))
                miss_det <- c(miss_det, paste0("SNP ", unknown, " is not in the SNP table"))
            }
            for (id in intersect(ids, known)) {
                p <- pos_by_id[[id]]; c_ <- chr_by_id[[id]]
                if (!identical(c_, as.character(r$chr[i])) || p < r$start[i] || p > r$end[i]) {
                    oob_key <- c(oob_key, r$region_id[i])
                    oob_det <- c(oob_det, paste0("SNP ", id, " at ", c_, ":", p,
                                                " lies outside ", r$chr[i], ":",
                                                r$start[i], "-", r$end[i]))
                }
            }
        }
        out$missing_snp <- violation("region_names_unknown_snp", "error", table_name, miss_key, miss_det)
        out$oob_snp     <- violation("region_snp_outside_bounds", "error", table_name, oob_key, oob_det)

        if (nrow(r) > length(known)) {
            out$more_regions <- violation("more_regions_than_snps", "error", table_name, table_name,
                                          paste0(nrow(r), " regions but only ", length(known), " SNPs"))
        }
    }

    do.call(combine_violations, unname(out))
}

# ── numeric: p-values ─────────────────────────────────────────────────────────

# pcols: names of p-value columns. Empty-string cells (the pipeline's "no hit"
# marker in selected_snps.tsv) are treated as absent, not as a bad p-value.
check_pvalues_table <- function(dt, pcols, table_name = "pvalues") {
    if (is.null(dt) || nrow(dt) == 0) return(no_violations())
    pcols <- intersect(pcols, names(dt))
    if (length(pcols) == 0) return(no_violations())

    key_col <- if ("SNPID" %in% names(dt)) as.character(dt$SNPID) else as.character(seq_len(nrow(dt)))
    out <- list()

    for (p in pcols) {
        v <- suppressWarnings(as.numeric(dt[[p]]))
        # is.na(NaN) is TRUE, so a !is.na() guard would swallow NaN. Test for it
        # explicitly and only then fall through to the plain-NA exclusion.
        bad_i <- which(is.nan(v) | (!is.na(v) & (is.infinite(v) | v < 0 | v > 1)))
        out[[paste0("range_", p)]] <- violation(
            "pvalue_out_of_range", "error", table_name,
            key_col[bad_i], paste0(p, " = ", dt[[p]][bad_i]))
    }

    if ("SNPID" %in% names(dt)) {
        dup <- data.table::data.table(id = key_col)[, .N, by = "id"][N > 1L]
        out$dup <- violation("snpid_duplicated", "error", table_name,
                             dup$id, paste0("appears ", dup$N, " times"))
    }

    if ("chr" %in% names(dt) && !is.character(dt$chr)) {
        out$chrtype <- violation("chr_not_character", "warn", table_name, table_name,
                                 paste0("chr column read as ", class(dt$chr)[1],
                                        " — a numeric-looking chromosome will not join against VCF headers"))
    }

    do.call(combine_violations, unname(out))
}

# selected_snps.tsv's min_pvalue is the minimum p over that SNP's SIGNIFICANT
# (trait, method) rows — combine_sigsnps.R:98-100 — NOT a row-wise minimum over
# anything in selected_snps.tsv itself, whose per-method columns hold TRAIT
# NAMES, not p-values. The validation source is therefore the long-format
# per-method *_sig_snps_*.tsv tables.
#
# sig_rows: long table with at least SNPID and pvalue, pooled across methods.
#   Pooling across several `adjust` variants of the same method is safe: a
#   given (SNPID, trait, method) carries the same p-value at every threshold,
#   and a looser threshold can only ADD rows with LARGER p, never a smaller one.
check_min_pvalue_against_sig_snps <- function(selected, sig_rows,
                                              table_name = "selected_snps.tsv",
                                              tol = 1e-9) {
    if (is.null(selected) || nrow(selected) == 0) return(no_violations())
    if (!all(c("SNPID", "min_pvalue") %in% names(selected))) return(no_violations())
    if (is.null(sig_rows) || nrow(sig_rows) == 0) return(no_violations())
    if (!all(c("SNPID", "pvalue") %in% names(sig_rows))) return(no_violations())

    obs <- data.table::data.table(
        SNPID  = as.character(sig_rows$SNPID),
        pvalue = suppressWarnings(as.numeric(sig_rows$pvalue))
    )[!is.na(pvalue), .(observed = min(pvalue)), by = "SNPID"]

    sel <- data.table::data.table(
        SNPID  = as.character(selected$SNPID),
        stated = suppressWarnings(as.numeric(selected$min_pvalue))
    )
    m <- obs[sel, on = "SNPID"]

    absent <- m[is.na(observed)]
    a <- violation("selected_snp_absent_from_sig_tables", "error", table_name,
                   absent$SNPID,
                   "SNP appears in selected_snps.tsv but in no *_sig_snps_*.tsv")

    cmp <- m[!is.na(observed)]
    bad <- cmp[is.na(stated) | abs(stated - observed) > tol * pmax(1, abs(observed))]
    b <- violation("min_pvalue_disagrees_with_sig_tables", "error", table_name,
                   bad$SNPID,
                   paste0("min_pvalue = ", bad$stated,
                          " but the smallest p across the sig tables is ", bad$observed))
    combine_violations(a, b)
}

# The trait names a module actually produced. In selected_snps.tsv the per-method
# columns are NAMED for methods and hold comma-separated TRAIT names as values,
# with "" (sometimes written with literal quotes) meaning "this method had no hit
# for this SNP". Extracting them is easy to get backwards — the column names are
# not traits — so it lives here, with fixtures, rather than in the CLI wrapper.
selected_snps_traits <- function(selected,
                                 fixed_cols = c("SNPID", "chr", "pos", "min_pvalue")) {
    if (is.null(selected) || nrow(selected) == 0) return(character(0))
    mcols <- setdiff(names(selected), fixed_cols)
    if (length(mcols) == 0) return(character(0))
    vals <- unlist(lapply(mcols, function(cn) as.character(selected[[cn]])), use.names = FALSE)
    vals <- gsub('"', "", vals)
    vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
    if (length(vals) == 0) return(character(0))
    unique(trimws(unlist(strsplit(vals, ",", fixed = TRUE), use.names = FALSE)))
}

# ── the "normalize early" promise ─────────────────────────────────────────────

# tables: named list of data.tables, each with a chr column.
# canonical: the run's chromosome set. When NULL it is the union of all tables,
#   so the check degrades to "no chr prefix" plus internal agreement.
#
# SUBSET, not equality: a genes table legitimately covers fewer chromosomes than
# its SNP table when a region contains no gene.
check_chromosome_names <- function(tables, canonical = NULL) {
    tables <- Filter(function(x) !is.null(x) && nrow(x) > 0 && "chr" %in% names(x), tables)
    if (length(tables) == 0) return(no_violations())

    seen <- lapply(tables, function(x) unique(as.character(x$chr)))
    if (is.null(canonical)) canonical <- sort(unique(unlist(seen, use.names = FALSE)))
    canonical <- as.character(canonical)
    out <- list()

    for (nm in names(seen)) {
        prefixed <- grep("^chr", seen[[nm]], value = TRUE, ignore.case = TRUE)
        out[[paste0("prefix_", nm)]] <- violation(
            "chromosome_name_not_normalized", "error", nm, prefixed,
            paste0("chromosome '", prefixed,
                   "' still carries a 'chr' prefix; the pipeline strips it in processing mode"))

        unknown <- setdiff(seen[[nm]], canonical)
        out[[paste0("unknown_", nm)]] <- violation(
            "chromosome_not_in_canonical_set", "error", nm, unknown,
            paste0("chromosome '", unknown, "' is absent from the run's chromosome set (",
                   paste(canonical, collapse = ", "), ")"))
    }
    do.call(combine_violations, unname(out))
}

# ── pipeline_summary.tsv ──────────────────────────────────────────────────────

# summary_dt: long format, columns step / metric / value.
check_summary_accounting <- function(summary_dt) {
    if (is.null(summary_dt) || nrow(summary_dt) == 0) return(no_violations())
    if (!all(c("step", "metric", "value") %in% names(summary_dt))) {
        return(violation("summary_schema", "error", "pipeline_summary.tsv", "pipeline_summary.tsv",
                         "expected long format with columns step, metric, value"))
    }
    s <- data.table::copy(summary_dt)
    out <- list()

    dup <- s[, .N, by = c("step", "metric")][N > 1L]
    # Guarded: paste0("x", character(0)) returns "x" (length ONE), so building a
    # key from an empty result would invent a violation.
    if (nrow(dup) > 0) {
        out$dup <- violation("summary_metric_duplicated", "error", "pipeline_summary.tsv",
                             paste0(dup$step, "/", dup$metric),
                             paste0("appears ", dup$N, " times; write_summary.R:30-36 should replace, not append"))
    }

    num <- function(step_, metric_) {
        v <- s[step == step_ & metric == metric_, value]
        if (length(v) != 1L) return(NA_real_)
        suppressWarnings(as.numeric(v[[1L]]))
    }

    total  <- num("processing", "samples_total")
    after  <- num("processing", "samples_after_filtering")
    rm_    <- num("processing", "samples_removed")
    het    <- num("processing", "samples_het_outliers_removed")
    rel    <- num("processing", "samples_removed_relatedness")
    if (!any(is.na(c(total, after, rm_, het, rel)))) {
        expected <- total - rm_ - het - rel
        if (!isTRUE(all.equal(expected, after))) {
            out$acct <- violation("sample_accounting_does_not_close", "error",
                                  "pipeline_summary.tsv", "processing/samples_after_filtering",
                                  paste0(total, " total - ", rm_, " removed - ", het,
                                         " het outliers - ", rel, " related = ", expected,
                                         ", but samples_after_filtering = ", after))
        }
    }

    with_c <- num("processing", "samples_with_coordinates")
    no_c   <- num("processing", "samples_dropped_missing_coordinates")
    if (!any(is.na(c(with_c, no_c, after))) && !isTRUE(all.equal(with_c + no_c, after))) {
        out$coord <- violation("coordinate_accounting_does_not_close", "error",
                               "pipeline_summary.tsv", "processing/samples_with_coordinates",
                               paste0(with_c, " with coordinates + ", no_c, " dropped = ",
                                      with_c + no_c, ", but samples_after_filtering = ", after))
    }

    preds <- s[step == "structure" & metric == "climate_predictors", value]
    n_pred <- num("structure", "n_climate_variables")
    if (length(preds) == 1L && !is.na(n_pred)) {
        listed <- length(strsplit(preds[[1L]], ",", fixed = TRUE)[[1]])
        if (listed != n_pred) {
            out$clim <- violation("climate_predictor_count_disagrees", "error",
                                  "pipeline_summary.tsv", "structure/n_climate_variables",
                                  paste0("climate_predictors lists ", listed,
                                         " predictors but n_climate_variables = ", n_pred))
        }
    }
    do.call(combine_violations, unname(out))
}

# A summary metric that states a row count must match the table it summarises.
check_summary_counts <- function(summary_dt, step, expected) {
    if (is.null(summary_dt) || nrow(summary_dt) == 0 || length(expected) == 0) {
        return(no_violations())
    }
    s <- data.table::copy(summary_dt)
    key <- character(0); det <- character(0)
    for (metric_ in names(expected)) {
        v <- s[step == step & metric == metric_, value]
        if (length(v) != 1L) next
        stated <- suppressWarnings(as.numeric(v[[1L]]))
        if (is.na(stated) || stated == expected[[metric_]]) next
        key <- c(key, paste0(step, "/", metric_))
        det <- c(det, paste0("summary says ", stated, ", the table has ", expected[[metric_]]))
    }
    violation("summary_count_disagrees_with_table", "error", "pipeline_summary.tsv", key, det)
}

# ── genes ─────────────────────────────────────────────────────────────────────

check_genes_table <- function(genes, regions = NULL, table_name = "genes_per_region") {
    if (is.null(genes) || nrow(genes) == 0) return(no_violations())
    out <- list()
    g <- data.table::copy(genes)

    if (all(c("gene_start", "gene_end") %in% names(g))) {
        g[, gene_start := as.integer(gene_start)][, gene_end := as.integer(gene_end)]
        bad <- g[gene_start > gene_end]
        out$bounds <- violation("gene_start_after_end", "error", table_name,
                                bad$gene_id, paste0("gene_start ", bad$gene_start,
                                                    " > gene_end ", bad$gene_end))
    }

    if (!is.null(regions) && nrow(regions) > 0 &&
        all(c("region_id") %in% names(g)) && "region_id" %in% names(regions)) {
        unknown <- setdiff(unique(as.character(g$region_id)),
                           unique(as.character(regions$region_id)))
        out$orphan <- violation("gene_references_unknown_region", "error", table_name,
                                unknown, "region_id is absent from the region table")
    }

    # exon_snp_count can never exceed the number of SNPs in the region. This is
    # what genes_in_regions.R:174 violates on real data (it counts features hit,
    # not SNPs); every cell is empty on SIMDATA, so it is a real-data check.
    if (!is.null(regions) && "snp_count" %in% names(regions) &&
        "region_id" %in% names(regions)) {
        snp_by_region <- stats::setNames(as.integer(regions$snp_count),
                                         as.character(regions$region_id))
        for (cn in intersect(c("exon_snp_count", "promoter_snp_count"), names(g))) {
            cnt <- suppressWarnings(as.integer(g[[cn]]))
            cap <- snp_by_region[as.character(g$region_id)]
            bad_i <- which(!is.na(cnt) & !is.na(cap) & cnt > cap)
            out[[cn]] <- violation("gene_snp_count_exceeds_region", "error", table_name,
                                   g$region_id[bad_i],
                                   paste0(cn, " = ", cnt[bad_i], " for gene ", g$gene_id[bad_i],
                                          " but the region holds only ", cap[bad_i], " SNPs"))
        }
    }
    do.call(combine_violations, unname(out))
}

# ── genetic offsets ───────────────────────────────────────────────────────────

check_offsets_table <- function(offsets, value_col = "genetic_offset",
                                table_name = "genetic_offset_site") {
    if (is.null(offsets) || nrow(offsets) == 0) return(no_violations())
    if (!(value_col %in% names(offsets))) {
        return(violation("offset_column_missing", "error", table_name, table_name,
                         paste0("no column named ", value_col)))
    }
    v   <- suppressWarnings(as.numeric(offsets[[value_col]]))
    key <- if ("sample" %in% names(offsets)) as.character(offsets$sample) else as.character(seq_along(v))

    bad_i <- which(is.na(v) | is.nan(v) | is.infinite(v))
    a <- violation("offset_not_finite", "error", table_name, key[bad_i],
                   paste0(value_col, " = ", offsets[[value_col]][bad_i]))
    neg_i <- which(!is.na(v) & is.finite(v) & v < 0)
    b <- violation("offset_negative", "error", table_name, key[neg_i],
                   paste0(value_col, " = ", v[neg_i]))
    combine_violations(a, b)
}

# ── stale outputs from a previous config ──────────────────────────────────────

# Snakemake never deletes an output that a config change orphaned, so a results
# tree accumulates files from earlier parameter choices that look exactly like
# current ones. A method with more than one `_sig_snps_{adjust}_{value}` variant
# on disk has at most one that matches the running config; the rest are stale
# and will be read by anything that globs.
#
# variants_by_method: named list, method -> character vector of adjust tokens
#   (e.g. list(RDA = c("bonf_0.01", "bonf_0.05"))).
check_single_threshold_variant <- function(variants_by_method, table_name = "sig_snps") {
    if (length(variants_by_method) == 0) return(no_violations())
    bad <- Filter(function(v) length(unique(v)) > 1L, variants_by_method)
    if (length(bad) == 0) return(no_violations())
    violation("multiple_threshold_variants_on_disk", "error", table_name,
              names(bad),
              vapply(bad, function(v) paste0(
                  "significant-SNP tables exist at ", length(unique(v)),
                  " different thresholds (", paste(sort(unique(v)), collapse = ", "),
                  "); at most one matches the running config, the rest are stale"),
                  character(1L), USE.NAMES = FALSE))
}

# ── cross-module referential integrity ────────────────────────────────────────

# Every value a downstream table references must exist upstream. This is the
# check that catches a mode re-run without its dependants being re-run —
# something Snakemake does not track and a golden file cannot see, because both
# files are individually self-consistent.
check_referential_integrity <- function(downstream_values, upstream_values,
                                        check_name = "downstream_references_unknown_value",
                                        table_name = "downstream", what = "value") {
    downstream_values <- unique(as.character(downstream_values[!is.na(downstream_values)]))
    upstream_values   <- unique(as.character(upstream_values[!is.na(upstream_values)]))
    if (length(downstream_values) == 0) return(no_violations())
    unknown <- setdiff(downstream_values, upstream_values)
    violation(check_name, "error", table_name, unknown,
              paste0(what, " '", unknown, "' is not present upstream (upstream has: ",
                     paste(upstream_values, collapse = ", "), ")"))
}

# ── table hygiene ─────────────────────────────────────────────────────────────

check_column_names_unique <- function(col_names, table_name) {
    col_names <- as.character(col_names)
    dup <- unique(col_names[duplicated(col_names)])
    violation("duplicate_column_names", "error", table_name, dup,
              paste0("column '", dup, "' appears ",
                     vapply(dup, function(d) sum(col_names == d), integer(1L)), " times"))
}

# filtering_summary.tsv stages must never gain samples or SNPs.
check_filtering_monotone <- function(filtering, table_name = "filtering_summary.tsv") {
    if (is.null(filtering) || nrow(filtering) < 2) return(no_violations())
    out <- list()
    for (cn in intersect(c("n_samples", "n_snps"), names(filtering))) {
        v <- suppressWarnings(as.numeric(filtering[[cn]]))
        bad_i <- which(diff(v) > 0) + 1L
        stage <- if ("stage" %in% names(filtering)) as.character(filtering$stage) else as.character(seq_along(v))
        out[[cn]] <- violation("filtering_stage_not_monotone", "error", table_name,
                               stage[bad_i],
                               paste0(cn, " rises from ", v[bad_i - 1L], " to ", v[bad_i],
                                      " at stage '", stage[bad_i], "'"))
    }
    do.call(combine_violations, unname(out))
}

# ── sig-SNP overlap columns ───────────────────────────────────────────────────

# overlap_traits must never name the SNP's OWN trait (sig_snps.R:153 filters for
# other-trait matches), and every id in overlap_snps must be a real SNP.
check_sig_snp_overlaps <- function(sig_snps, known_snpids = NULL,
                                   table_name = "sig_snps") {
    if (is.null(sig_snps) || nrow(sig_snps) == 0) return(no_violations())
    if (!all(c("overlap_traits", "overlap_snps") %in% names(sig_snps))) return(no_violations())
    out <- list()
    s <- sig_snps

    if ("trait" %in% names(s)) {
        self_key <- character(0); self_det <- character(0)
        for (i in seq_len(nrow(s))) {
            ot <- s$overlap_traits[i]
            if (is.na(ot) || !nzchar(ot)) next
            traits <- strsplit(ot, ",", fixed = TRUE)[[1]]
            if (as.character(s$trait[i]) %in% traits) {
                self_key <- c(self_key, as.character(s$SNPID[i]))
                self_det <- c(self_det, paste0("overlap_traits = '", ot,
                                               "' includes the SNP's own trait '", s$trait[i], "'"))
            }
        }
        out$self <- violation("overlap_traits_includes_own_trait", "error", table_name,
                              self_key, self_det)
    }

    if (!is.null(known_snpids)) {
        known <- unique(as.character(known_snpids))
        bad_key <- character(0); bad_det <- character(0)
        for (i in seq_len(nrow(s))) {
            os <- s$overlap_snps[i]
            if (is.na(os) || !nzchar(os)) next
            ids <- unique(strsplit(os, ",", fixed = TRUE)[[1]])
            unknown <- setdiff(ids, known)
            if (length(unknown) > 0) {
                bad_key <- c(bad_key, rep(as.character(s$SNPID[i]), length(unknown)))
                bad_det <- c(bad_det, paste0("overlap_snps names '", unknown,
                                             "', which is not a SNP in this run"))
            }
        }
        out$unknown <- violation("overlap_snps_names_unknown_snp", "error", table_name,
                                 bad_key, bad_det)
    }
    do.call(combine_violations, unname(out))
}
