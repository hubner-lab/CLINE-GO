# snp_set_select.R — the ONE implementation of a curated SNP set's selection.
#
# Shared by the Shiny app (sys.source()d into the clinego.app namespace by
# zzz.R / dev.R via CLINEGO_SHARED_LIBS) and the pipeline (scripts/build_snp_set.R,
# which builds the GEA.snp_sets / GWAS.snp_sets declared in config). A set saved
# from the GEA tab and a set declared in config go through the same pair filter
# (filter_to_pairs) and combine (combine_sigsnps_long). PENDING (plan
# staged-plotting-newt, Phase 3): the app's save_snp_set() still has its own
# collapse, ordered by chr, pos only; until it calls collapse_snp_set() the two
# stores can order rows differently when two SNP IDs share a position.
#
# This is NOT the pipeline's scripts/R/lib/combine_sigsnps.R. That one feeds the
# static selected_snps.tsv tables, takes a predictor whitelist, and diverges from
# this one on Cross-method min_pvalue (tests/testthat/test-equivalence-app-pipeline.R,
# "KNOWN DIVERGENCE: app Cross-method keeps rows that never matched"). The two
# deliberately coexist: different names, different callers.
#
# Requires: data.table ATTACHED (bare := / .() / column names inside [ ]);
# compute_pval_threshold() from scripts/R/utils/pval_threshold.R sourced first.
# yaml and jsonlite are only called namespace-qualified.

# ── Strategy names ────────────────────────────────────────────────────────────

# Normalize legacy and alias strategy names to current canonical names.
.normalize_strategy <- function(s) {
    switch(s,
        # Legacy pipeline names
        All          = "Union",
        Sum          = "Union",
        Overlap      = "Cross-method",
        MethodOverlap = "Cross-method per-trait",
        PairOverlap  = "Cross-method per-trait",
        s  # pass through "Union", "Cross-method", "Cross-method per-trait" unchanged
    )
}

# ── Combine ───────────────────────────────────────────────────────────────────

# Combine per-method sig SNPs using a specified strategy.
#
# @param sigsnps_list Named list of data.tables, one per method.
#   Each must have columns: SNPID, chr, pos, pvalue, method, trait.
# @param strategy "Union", "Cross-method", "Cross-method per-trait", a legacy
#   alias (see .normalize_strategy), or a method name (single-method passthrough).
# @param gap Integer. Max bp distance for spatial overlap (Cross-method variants).
# @return Long-format data.table with SNPID, chr, pos, pvalue, method, trait,
#   min_pvalue. Empty skeleton (no min_pvalue) if nothing passes.
combine_sigsnps_long <- function(sigsnps_list, strategy = "All", gap = 200000L) {
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    strategy <- .normalize_strategy(strategy)

    # Drop empty methods
    sigsnps_list <- sigsnps_list[vapply(sigsnps_list, nrow, integer(1)) > 0]
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    methods_vec <- names(sigsnps_list)
    gap <- as.integer(gap)

    result <- if (strategy %in% methods_vec) {
        # ── Single-method ───────────────────────────────────────────────────────
        data.table::copy(sigsnps_list[[strategy]])

    } else if (strategy == "Union") {
        # ── Union: all sig SNPs from all methods ────────────────────────────────
        data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)

    } else if (strategy == "Cross-method") {
        # ── Cross-method spatial overlap (trait-agnostic) ───────────────────────
        # A SNP passes if any partner exists in another method within gap bp.
        selected_ids <- character(0)
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                ids <- .spatial_overlap_ids(sigsnps_list[[n1]], sigsnps_list[[n2]], gap)
                selected_ids <- c(selected_ids, ids)
            }
        }
        selected_ids <- unique(selected_ids)
        if (length(selected_ids) == 0) return(.empty_sigsnps())
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_snps[SNPID %in% selected_ids]

    } else if (strategy == "Cross-method per-trait") {
        # ── Cross-method spatial overlap, same-trait only ───────────────────────
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_traits <- unique(all_snps$trait)
        selected_rows <- list()
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                for (bio in all_traits) {
                    m1_bio <- sigsnps_list[[n1]][trait == bio]
                    m2_bio <- sigsnps_list[[n2]][trait == bio]
                    if (nrow(m1_bio) == 0 || nrow(m2_bio) == 0) next
                    ids1 <- .spatial_overlap_ids(m1_bio, m2_bio, gap)
                    ids2 <- .spatial_overlap_ids(m2_bio, m1_bio, gap)
                    if (length(ids1) > 0)
                        selected_rows <- c(selected_rows, list(m1_bio[SNPID %in% ids1]))
                    if (length(ids2) > 0)
                        selected_rows <- c(selected_rows, list(m2_bio[SNPID %in% ids2]))
                }
            }
        }
        if (length(selected_rows) == 0) return(.empty_sigsnps())
        unique(data.table::rbindlist(selected_rows, use.names = TRUE, fill = TRUE))

    } else {
        stop(paste0("Unknown combine strategy: '", strategy,
                    "'. Valid values: Union, Cross-method, Cross-method per-trait, or a method name."))
    }

    if (is.null(result) || nrow(result) == 0) return(.empty_sigsnps())

    # Deduplicate — same SNP can appear in multiple methods/traits
    result <- unique(result, by = c("SNPID", "method", "trait"))

    # Add min_pvalue column for combined Manhattan y-placement
    # (background PNG y-axis spans the minimum p-value per SNP across all methods)
    result[, min_pvalue := min(pvalue, na.rm = TRUE), by = "SNPID"]

    result
}

# Find SNP IDs from dt1 that have a partner in dt2 within gap bp on same chr.
.spatial_overlap_ids <- function(dt1, dt2, gap) {
    if (nrow(dt1) == 0 || nrow(dt2) == 0) return(character(0))

    # Expand dt1 to windows [pos-gap, pos+gap]; dt2 stays as points [pos, pos]
    win <- data.table::data.table(
        chr   = as.character(dt1$chr),
        start = pmax(1L, as.integer(dt1$pos) - gap),
        end   = as.integer(dt1$pos) + gap,
        SNPID = dt1$SNPID
    )
    pts <- data.table::data.table(
        chr   = as.character(dt2$chr),
        start = as.integer(dt2$pos),
        end   = as.integer(dt2$pos),
        SNPID = dt2$SNPID
    )
    data.table::setkey(win, chr, start, end)
    data.table::setkey(pts, chr, start, end)
    ov <- data.table::foverlaps(pts, win, type = "within", nomatch = NULL)
    if (nrow(ov) == 0) return(character(0))
    # Return IDs from dt1 (window side)
    unique(ov$SNPID)
}

# Empty sig SNPs skeleton.
.empty_sigsnps <- function() {
    data.table::data.table(
        SNPID  = character(),
        chr    = character(),
        pos    = integer(),
        pvalue = numeric(),
        method = character(),
        trait  = character()
    )
}

# Empty sig SNPs skeleton with region_id column (for interactive combine output).
.empty_sigsnps_assoc <- function() {
    data.table::data.table(
        SNPID     = character(),
        chr       = character(),
        pos       = integer(),
        pvalue    = numeric(),
        method    = character(),
        trait     = character(),
        region_id = character()
    )
}

# ── Pair selection ────────────────────────────────────────────────────────────

# Split "trait::method" keys at the FIRST "::".
# @return data.frame(pair, trait, method), one row per key, input order.
split_snp_set_pairs <- function(pairs) {
    pairs <- as.character(pairs)
    trait <- sub("::.*$", "", pairs)
    ok <- grepl("::", pairs, fixed = TRUE) & nzchar(trait)
    if (!all(ok)) {
        stop("malformed trait::method pair(s): ",
             paste(pairs[!ok], collapse = ", "))
    }
    data.frame(pair = pairs, trait = trait,
               method = substring(pairs, nchar(trait) + 3L),
               stringsAsFactors = FALSE)
}

# Keep, per method, only the rows whose trait is paired with that method.
#
# @param sigsnps_by_method Named list method -> data.table (needs a trait column).
# @param pairs Character vector of "trait::method" keys.
# A malformed key (no "trait::" part) matches nothing rather than stopping:
# this runs inside a Shiny reactive on whatever the matrix JS sent.
# @return The filtered list; methods left with no rows are dropped.
filter_to_pairs <- function(sigsnps_by_method, pairs) {
    pairs <- as.character(pairs)
    sp <- split_snp_set_pairs(pairs[grepl("^.+::", pairs)])
    filtered <- lapply(names(sigsnps_by_method), function(m) {
        paired_traits <- sp$trait[sp$method == m]
        if (length(paired_traits) == 0) return(NULL)
        dt <- sigsnps_by_method[[m]]
        if (nrow(dt) == 0) return(NULL)
        dt[trait %in% paired_traits]
    })
    names(filtered) <- names(sigsnps_by_method)
    Filter(function(x) !is.null(x) && nrow(x) > 0, filtered)
}

# Collapse a long combined table to the SNP set store's four columns.
#
# UNIQUE by SNPID, min p across methods/traits, ordered chr, pos, SNPID (the
# SNPID tiebreak makes the order independent of the method list's order, so
# two runs that combine the same cells write byte-identical files).
# @return data.table(SNPID, chr, pos, min_pvalue); 0 rows for an empty input.
collapse_snp_set <- function(long_dt) {
    dt <- data.table::as.data.table(long_dt)
    if (nrow(dt) == 0) {
        return(data.table::data.table(SNPID = character(), chr = character(),
                                      pos = integer(), min_pvalue = numeric()))
    }
    uniq <- dt[, .(chr = chr[1L], pos = pos[1L],
                   min_pvalue = min(min_pvalue, na.rm = TRUE)), by = SNPID]
    data.table::setcolorder(uniq, c("SNPID", "chr", "pos", "min_pvalue"))
    data.table::setorder(uniq, chr, pos, SNPID)
    uniq
}

# ── One-call selection ────────────────────────────────────────────────────────

# Threshold every declared cell, combine the selected pairs, collapse.
#
# @param pvals_by_method Named list method -> wide data.table (SNPID, chr, pos,
#   one p-value column per trait), as load_all_method_pvalues() returns.
# @param cells Named list "trait::method" -> list(type =, value =), the RESOLVED
#   rule of every cell in `pairs` (precedence already applied by the caller).
# @param pairs Character vector of "trait::method" keys to include.
# @param strategy,gap As combine_sigsnps_long().
# @return list(snps  = collapse_snp_set() table,
#              long  = combine_sigsnps_long() table,
#              cells = data.frame(pair, type, value, threshold, status, n_sig))
#   A cell whose rule cannot run (qval under 10 tests, top N above the test
#   count) contributes no SNPs, exactly as in the app; its status says why.
select_snp_set <- function(pvals_by_method, cells, pairs,
                           strategy = "Union", gap = 100000L) {
    if (length(pairs) == 0) stop("select_snp_set: no trait::method pairs to select")
    if (anyDuplicated(pairs)) stop("select_snp_set: duplicate pairs: ",
                                   paste(unique(pairs[duplicated(pairs)]), collapse = ", "))
    sp <- split_snp_set_pairs(pairs)

    no_rule <- setdiff(sp$pair, names(cells))
    if (length(no_rule)) stop("select_snp_set: no rule for pair(s): ",
                              paste(no_rule, collapse = ", "))
    no_method <- setdiff(unique(sp$method), names(pvals_by_method))
    if (length(no_method)) stop("select_snp_set: no p-value table for method(s): ",
                                paste(no_method, collapse = ", "))
    no_trait <- sp$pair[!mapply(function(t, m) t %in% names(pvals_by_method[[m]]),
                                sp$trait, sp$method)]
    if (length(no_trait)) stop("select_snp_set: trait column missing for pair(s): ",
                               paste(no_trait, collapse = ", "))

    rows <- vector("list", nrow(sp))
    summary <- data.frame(pair = sp$pair, type = NA_character_, value = NA_real_,
                          threshold = NA_real_, status = NA_character_, n_sig = 0L,
                          stringsAsFactors = FALSE)
    for (i in seq_len(nrow(sp))) {
        tr <- sp$trait[i]; m <- sp$method[i]
        rule <- cells[[sp$pair[i]]]
        pv_dt <- pvals_by_method[[m]]
        pvec <- pv_dt[[tr]]
        res <- tryCatch(
            compute_pval_threshold(pvec, rule$type, rule$value),
            error = function(e) list(status = "error", threshold = NA_real_,
                                     message = conditionMessage(e))
        )
        summary$type[i]   <- rule$type
        summary$value[i]  <- as.numeric(rule$value)
        summary$status[i] <- res$status
        if (!identical(res$status, "ok") || is.na(res$threshold)) next
        summary$threshold[i] <- res$threshold
        # Inclusive '<=': compute_pval_threshold() returns a cutoff that is a
        # member of the call set under 'qval'/'top' (pval_threshold.R:58).
        mask <- pvec <= res$threshold & !is.na(pvec)
        summary$n_sig[i] <- sum(mask)
        if (!any(mask)) next
        rows[[i]] <- pv_dt[mask, .(SNPID = SNPID, chr = chr, pos = pos,
                                   pvalue = .SD[[tr]], method = m, trait = tr),
                           .SDcols = tr]
    }

    flat <- data.table::rbindlist(rows)
    long <- if (nrow(flat) == 0) .empty_sigsnps() else
        combine_sigsnps_long(split(flat, by = "method", keep.by = TRUE),
                             strategy = strategy, gap = as.integer(gap))
    list(snps = collapse_snp_set(long), long = long, cells = summary)
}

# ── Declaration / record I/O ──────────────────────────────────────────────────
#
# declaration.json  the resolved inputs of one set, deterministic (no timestamp).
#                   Written by the pipeline's snp_set_declare rule.
# snp_set.yaml      the record: the declaration plus k_best, n_snps, created
#                   (and config_snippet for a set saved from Shiny). Written next
#                   to selected_snps.tsv by build_snp_set.R and by save_snp_set().
# A record is a superset of a declaration, so read_snp_set_decl() reads either.

SNP_SET_DECL_FIELDS <- c("name", "origin", "source_module", "regime",
                         "strategy", "clumping_distance", "pairs", "cells")

# Coerce a parsed declaration/record to canonical R types and check it.
# JSON/YAML both collapse a length-1 array to a scalar on some writers, so
# `pairs` is unlisted rather than trusted to arrive as a list.
normalize_snp_set_decl <- function(x) {
    missing <- setdiff(SNP_SET_DECL_FIELDS, names(x))
    if (length(missing)) stop("SNP set declaration lacks field(s): ",
                              paste(missing, collapse = ", "))
    pairs <- as.character(unlist(x$pairs, use.names = FALSE))
    if (length(pairs) == 0) stop("SNP set declaration '", x$name, "' has no pairs")
    cells <- x$cells
    if (!is.list(cells) || is.null(names(cells)) || any(!nzchar(names(cells)))) {
        stop("SNP set declaration '", x$name, "': cells must be a trait::method-keyed map")
    }
    cells <- lapply(names(cells), function(k) {
        r <- cells[[k]]
        ty <- as.character(r$type %||% NA_character_)
        va <- suppressWarnings(as.numeric(r$value %||% NA_real_))
        if (length(ty) != 1L || is.na(ty) || length(va) != 1L || is.na(va)) {
            stop("SNP set declaration '", x$name, "': cell ", k,
                 " needs a scalar type and a numeric value")
        }
        list(type = ty, value = va)
    })
    names(cells) <- names(x$cells)
    split_snp_set_pairs(pairs)   # stops on a malformed key

    out <- x
    out$name              <- as.character(x$name)
    out$origin            <- as.character(x$origin)
    out$source_module     <- as.character(x$source_module)
    out$regime            <- as.character(x$regime)
    out$strategy          <- as.character(x$strategy)
    out$clumping_distance <- as.integer(x$clumping_distance)
    out$pairs             <- pairs
    out$cells             <- cells
    out
}

# Read a declaration (.json) or a record (.yaml/.yml) into normalized form.
read_snp_set_decl <- function(path) {
    if (!file.exists(path)) stop("SNP set declaration not found: ", path)
    raw <- if (grepl("\\.json$", path, ignore.case = TRUE)) {
        jsonlite::read_json(path, simplifyVector = FALSE)
    } else {
        yaml::read_yaml(path)
    }
    normalize_snp_set_decl(raw)
}

# Write a record atomically (tempfile in the target dir + rename).
#
# Fields are written in a fixed order: the declaration fields, then the rest in
# the order given. `pairs` is written as a sequence even when it holds one key,
# and numbers keep 15 significant digits (yaml's default of 7 would round a
# threshold value on its way to the file).
write_snp_set_record <- function(record, path) {
    rec <- normalize_snp_set_decl(record)
    rec <- rec[c(SNP_SET_DECL_FIELDS, setdiff(names(rec), SNP_SET_DECL_FIELDS))]
    rec$pairs <- as.list(rec$pairs)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".yaml")
    ok <- FALSE
    on.exit(if (!ok && file.exists(tmp)) unlink(tmp), add = TRUE)
    writeLines(yaml::as.yaml(rec, precision = 15L), tmp)
    if (!file.rename(tmp, path)) stop("could not move SNP set record into place: ", path)
    ok <- TRUE
    invisible(path)
}
