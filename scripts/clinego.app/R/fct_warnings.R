# fct_warnings.R — warnings about plain data problems (release A2).
#
# Each rule is a pure function: data.tables in, NULL or a zero-row table when there is
# nothing to say, otherwise what the warning names. The modules render the result with
# warning_note(). A rule fires only on an exact equality or on the pipeline's own rule,
# never on a cut chosen by eye. A label is at most WARNING_LABEL_MAX_WORDS English words.
# A warning informs, it does not teach (A3): its note lists what was found, item by item
# (which scan, predictor, trait or sample, which problem, its value), then at most one short
# clause and the Config line, both from HELP_NOTES (fct_help_notes.R) under its
# WARNING_NOTE_IDS id. Both are tested: warning_note() returns NULL for an unknown id, so a
# typo would drop a warning without a word.

#' HELP_NOTES ids the warnings render through.
#' @noRd
WARNING_NOTE_IDS <- c(
    "warn_manhattan_stripe", "warn_ld_decay", "warn_climate_predictors",
    "warn_climate_dropped_samples", "warn_predictor_chips", "warn_trait_chips",
    "warn_cross_entropy", "warn_dbmem_fallback"
)

#' Longest label a warning badge may carry, in words.
#' @noRd
WARNING_LABEL_MAX_WORDS <- 8L

#' An amber badge for one warning. Its tooltip: what was found, one item per line, then
#' the entry's short clause (if any) and its Config line. No "what it means" block.
#' @param id    One of WARNING_NOTE_IDS.
#' @param label The badge text, at most WARNING_LABEL_MAX_WORDS words.
#' @param items What this project's tables say, one string per item, or NULL.
#' @param class Badge class; neutral ("bg-secondary") for a state that is unknown
#'   rather than wrong.
#' @noRd
warning_note <- function(id, label, items = NULL, class = "bg-warning text-dark") {
    entry <- HELP_NOTES[[id]]
    if (is.null(entry)) return(NULL)
    body <- htmltools::tagList(
        if (length(items) > 0)
            htmltools::tags$ul(class = "warn-items", lapply(items, htmltools::tags$li)),
        if (!is.null(entry$desc)) htmltools::p(class = "warn-why", entry$desc),
        help_note_config(entry$config)
    )
    filter_note(label = label, body = body, class = class, icon = "exclamation-triangle-fill")
}

#' At most `n` items, then one "and N more <noun>s" item.
#' @noRd
capped_items <- function(items, n, noun) {
    if (length(items) <= n) return(items)
    c(items[seq_len(n)], paste0("and ", n_noun(length(items) - n, noun), "."))
}

#' Names joined by ", ", the first `n` only, then "and N more".
#' @noRd
name_list <- function(x, n = 12L) {
    if (length(x) <= n) return(paste(x, collapse = ", "))
    paste0(paste(x[seq_len(n)], collapse = ", "), " and ", length(x) - n, " more")
}

#' Number of words in a label (what WARNING_LABEL_MAX_WORDS counts).
#' @noRd
label_words <- function(label) length(strsplit(trimws(label), "\\s+")[[1]])

#' "1 sample" / "3 samples".
#' @noRd
n_noun <- function(n, singular, plural = paste0(singular, "s")) {
    paste(n, if (n == 1) singular else plural)
}

# ── Manhattan stripes ─────────────────────────────────────────────────────────

#' Stripes on a Manhattan: one method x trait's significant SNPs that tie exactly on p
#' across chromosomes.
#'
#' A single-marker test gives two markers the same statistic, and so the same p, when
#' their genotype vectors are the same or allele-flipped. A p-value shared by significant
#' SNPs on two or more chromosomes is therefore one genotype pattern counted as many hits
#' (a lineage or site indicator, or copies the imputation made alike); on a Manhattan it is
#' a horizontal stripe. Ties are read at signif(p, 10). p >= 1, and p <= 0 (an underflow,
#' where unrelated markers tie), are not read. Trait columns of one method whose hits and
#' p-values are identical (twin predictors) count once, named together. The caller skips
#' WZA windows: a window's p is not a single-marker test.
#'
#' @param sig_list Named list, method -> data.table(SNPID, chr, pvalue, trait): the
#'   significant SNPs at the live threshold (compute_method_sigsnps_cached()).
#' @return data.table, one row per method and trait group with a stripe, most stripe hits
#'   first: method, traits (" = "-joined), n_hits, n_stripes, n_stripe_hits, largest
#'   (hits in the largest tie), largest_p (its p), n_chr (its chromosomes). Zero rows when
#'   there is none.
#' @noRd
manhattan_stripes <- function(sig_list) {
    rows <- list()
    for (m in names(sig_list)) {
        dt <- sig_list[[m]]
        if (is.null(dt) || nrow(dt) == 0) next
        dt <- data.table::as.data.table(dt)
        dt <- dt[is.finite(dt$pvalue) & dt$pvalue > 0 & dt$pvalue < 1, ]
        if (nrow(dt) == 0) next
        dt <- data.table::data.table(snp = as.character(dt$SNPID), chr = as.character(dt$chr),
                                     trait = as.character(dt$trait),
                                     p = signif(dt$pvalue, 10))
        data.table::setorderv(dt, c("trait", "snp"))
        # one signature per trait, its hits with their p: identical columns share it
        sig <- dt[, list(key = paste(snp, p, collapse = "|")), by = "trait"]
        for (trs in split(sig$trait, sig$key)) {
            ties <- dt[dt$trait == trs[1], list(n = .N, n_chr = data.table::uniqueN(chr)),
                       by = "p"]
            st <- ties[ties$n_chr >= 2, ]
            if (nrow(st) == 0) next
            big <- which.max(st$n)
            rows[[length(rows) + 1]] <- data.table::data.table(
                method = m, traits = paste(sort(trs), collapse = " = "),
                n_hits = sum(ties$n), n_stripes = nrow(st), n_stripe_hits = sum(st$n),
                largest = st$n[big], largest_p = st$p[big], n_chr = st$n_chr[big])
        }
    }
    if (length(rows) == 0) {
        return(data.table::data.table(method = character(), traits = character(),
                                      n_hits = integer(), n_stripes = integer(),
                                      n_stripe_hits = integer(), largest = integer(),
                                      largest_p = numeric(), n_chr = integer()))
    }
    out <- data.table::rbindlist(rows)
    out[order(-out$n_stripe_hits, out$method, out$traits), ]
}

#' The stripe rows of one method and trait (the trait may sit in a twin group).
#' @noRd
stripes_for <- function(stripes, method, trait) {
    if (is.null(stripes) || nrow(stripes) == 0 || is.null(method) || is.null(trait))
        return(NULL)
    in_group <- vapply(strsplit(stripes$traits, " = ", fixed = TRUE),
                       function(t) trait %in% t, logical(1))
    # computed outside [ ]: inside it a bare `method` is the column, not the argument
    keep <- which(stripes$method == method & in_group)
    stripes[keep, ]
}

#' Badge label for one stripe row.
#' @noRd
stripe_label <- function(row) {
    if (row$n_stripes == 1) paste0("Stripe: ", row$n_stripe_hits, " hits share one genotype")
    else paste0("Stripes: ", row$n_stripe_hits, " hits share ", row$n_stripes, " genotypes")
}

#' The stripe warning for a Manhattan card: NULL in the WZA regime or without a stripe.
#' @param stripes manhattan_stripes().
#' @param wza     TRUE in the WZA regime.
#' @noRd
stripe_note <- function(stripes, wza = FALSE) {
    if (isTRUE(wza) || is.null(stripes) || nrow(stripes) == 0) return(NULL)
    # one item per scan: the scan, its ties, the largest tie's p and hits
    items <- vapply(seq_len(nrow(stripes)), function(i) {
        s <- stripes[i, ]
        paste0(s$method, ", ", s$traits, ": ", s$n_stripe_hits, " of ", s$n_hits,
               " hits in ", n_noun(s$n_stripes, "tie"), "; largest ",
               if (!is.null(s$largest_p)) paste0("p = ", formatC(s$largest_p, format = "g",
                                                                 digits = 3), ", "),
               s$largest, " hits on ", s$n_chr, " chromosomes")
    }, character(1))
    warning_note("warn_manhattan_stripe", stripe_label(stripes[1, ]),
                 items = capped_items(items, 6L, "scan"))
}

# ── LD decay ──────────────────────────────────────────────────────────────────

#' The LD-decay warning: the genome-wide curve of 'All' is unreliable or not assessed.
#' @param gw ld_decay_genome_wide_status().
#' @noRd
ld_decay_note <- function(gw) {
    if (is.null(gw) || identical(gw$status, "reliable")) return(NULL)
    # short: it sits in the LD-decay cards' narrow headers, which name the curve
    label <- if (identical(gw$status, "unreliable")) "unreliable" else "not assessed"
    warning_note("warn_ld_decay", label, items = paste0("'All', genome-wide: ", gw$reason))
}

#' Theme of the Structure tab's LD-decay KPI (8bbdee): green only for a reliable curve.
#' @param gw ld_decay_genome_wide_status(), NULL without a genome-wide 'All' row.
#' @noRd
ld_decay_kpi_theme <- function(gw) {
    if (is.null(gw)) "secondary"
    else if (identical(gw$status, "reliable")) "success"
    else "warning"
}

# ── Climate predictors: identical and two-valued ──────────────────────────────

#' Identical and two-valued climate predictors.
#'
#' Read over the samples climate was extracted for (climate_present_site.tsv, one row per
#' sample). Identical: 1 - |r| < 1e-9 between two predictors. Two-valued ("binary"):
#' exactly two distinct values, i.e. the predictor splits the sites into two groups.
#' Neither is read below three sites, where every predictor is two-valued and every pair
#' perfectly correlated. A predictor with one value is check_climate_variance's
#' (climate_invariant_predictors.tsv) and is left to it.
#'
#' @param clim data.table: sample + one numeric column per predictor.
#' @param meta data.table with sample and site (Processing/tables/metadata.tsv).
#' @return data.table(predictor, flag = "identical" | "binary", detail): its twins
#'   (comma-separated) or its two values. Zero rows when nothing is flagged.
#' @noRd
predictor_flags <- function(clim, meta) {
    empty <- data.table::data.table(predictor = character(), flag = character(),
                                    detail = character())
    if (is.null(clim) || is.null(meta) || nrow(clim) == 0 ||
        !all(c("sample", "site") %in% names(meta)) || !"sample" %in% names(clim))
        return(empty)
    clim <- data.table::as.data.table(clim)
    n_sites <- data.table::uniqueN(meta$site[meta$sample %in% clim$sample])
    if (n_sites < 3) return(empty)
    preds <- setdiff(names(clim), c("sample", "site", "latitude", "longitude"))
    preds <- preds[vapply(preds, function(p) is.numeric(clim[[p]]), logical(1))]
    if (length(preds) == 0) return(empty)

    values <- lapply(preds, function(p) sort(unique(stats::na.omit(clim[[p]]))))
    names(values) <- preds
    binary <- preds[lengths(values) == 2]
    rows <- lapply(binary, function(p) data.table::data.table(
        predictor = p, flag = "binary",
        detail = paste(format(values[[p]], digits = 4, trim = TRUE), collapse = " and ")))

    varying <- preds[lengths(values) >= 2]
    if (length(varying) >= 2) {
        r <- suppressWarnings(stats::cor(as.matrix(clim[, varying, with = FALSE]),
                                         use = "pairwise.complete.obs"))
        same <- !is.na(r) & (1 - abs(r)) < 1e-9
        diag(same) <- FALSE
        for (p in varying[rowSums(same) > 0]) {
            rows[[length(rows) + 1]] <- data.table::data.table(
                predictor = p, flag = "identical",
                detail = paste(varying[same[p, ]], collapse = ", "))
        }
    }
    if (length(rows) == 0) return(empty)
    data.table::rbindlist(rows)
}

#' The predictors a project's config selects, as the chips read it: Climate.predictors,
#' or every chip offered when it is empty.
#' @noRd
selected_predictors <- function(config) {
    sel <- parse_csv_list(config_get(config, "Climate", "predictors", default = NULL))
    if (length(sel) == 0) predictor_chip_choices(config) else sel
}

#' Names in chip order: bio_<N> by N, any other name after them, alphabetically.
#' @noRd
chip_order <- function(x) {
    n <- suppressWarnings(as.integer(sub("^bio_", "", x)))
    n[!grepl("^bio_[0-9]+$", x)] <- NA_integer_
    x[order(is.na(n), n, x)]
}

#' Twin sets among the selected predictors (each set in chip order).
#' @param flags    predictor_flags().
#' @param selected Selected predictors.
#' @return List of character vectors, each of two or more selected identical predictors.
#' @noRd
selected_twin_sets <- function(flags, selected) {
    idn <- flags[flags$flag == "identical" & flags$predictor %in% selected, ]
    sets <- list()
    seen <- character(0)
    for (i in seq_len(nrow(idn))) {
        p <- idn$predictor[i]
        if (p %in% seen) next
        set <- chip_order(intersect(c(p, parse_csv_list(idn$detail[i])), selected))
        seen <- c(seen, set)
        if (length(set) >= 2) sets[[length(sets) + 1]] <- set
    }
    firsts <- vapply(sets, `[[`, "", 1L)
    sets[match(chip_order(firsts), firsts)]
}

#' Tab warning (W4): one badge for the selected predictors that are identical to another
#' selected one or take two values; its note has one item per twin set and per two-valued
#' predictor.
#' @noRd
climate_predictors_note <- function(flags, selected) {
    sets <- selected_twin_sets(flags, selected)
    bin  <- flags[flags$flag == "binary" & flags$predictor %in% selected, ]
    if (length(sets) == 0 && nrow(bin) == 0) return(NULL)
    n <- length(unique(c(unlist(sets), bin$predictor)))
    kinds <- c(if (length(sets) > 0) "identical", if (nrow(bin) > 0) "two-valued")
    warning_note("warn_climate_predictors",
        paste0(n_noun(n, "selected predictor"), ": ", paste(kinds, collapse = " or ")),
        items = c(vapply(sets, function(s) paste0(paste(s, collapse = " = "),
                                                  ": identical at every sample"), ""),
                  paste0(bin$predictor, ": two values across the sites, ", bin$detail)))
}

#' Tab counter of invariant predictors over the selected ones only (e83e54).
#' @param invariant Predictors climate_invariant_predictors.tsv names.
#' @param selected  Selected predictors.
#' @noRd
selected_invariant <- function(invariant, selected) intersect(selected, invariant)

# ── Chip marks ────────────────────────────────────────────────────────────────

#' A flagged chip's problems, e.g. "identical to bio_11 at every sample".
#' @param flags Its rows of predictor_flags() / trait_flags().
#' @noRd
chip_problems <- function(flags) {
    parts <- vapply(seq_len(nrow(flags)), function(i) {
        switch(flags$flag[i],
            identical = paste0("identical to ", flags$detail[i], " at every sample"),
            binary    = paste0("two values across the sites, ", flags$detail[i]),
            flags$detail[i])
    }, character(1))
    paste(parts, collapse = "; ")
}

#' Tooltip text for one flagged chip.
#' @param name  Chip name.
#' @param flags Its rows of predictor_flags() / trait_flags().
#' @noRd
chip_flag_title <- function(name, flags) paste0(name, ": ", chip_problems(flags))

#' Legend badge under a chip grid: each marked chip and its problem, one item per chip
#' (no corner icons: a marked chip is orange, A3). Neutral, amber once a marked chip is
#' selected (the chip script keeps that in step with clicks).
#' @param flags     The grid's predictor_flags() / trait_flags() rows.
#' @param struck    Names struck as invariant.
#' @param selected  Names selected when the grid renders.
#' @param legend_id Element id the chip script toggles.
#' @param id        "warn_predictor_chips" or "warn_trait_chips".
#' @param noun      "predictor" or "trait".
#' @noRd
chip_legend <- function(flags, struck, selected, legend_id, id, noun) {
    marked <- unique(c(flags$predictor, struck))
    if (length(marked) == 0) return(NULL)
    amber <- any(unique(flags$predictor) %in% selected)
    marked <- chip_order(marked)
    struck_reason <- if (identical(noun, "trait")) "does not vary or has no values"
                     else "no variance"
    items <- vapply(marked, function(m) {
        if (m %in% struck) paste0(m, ": ", struck_reason, ", cannot be selected")
        else paste0(m, ": ", chip_problems(flags[flags$predictor == m, ]))
    }, character(1), USE.NAMES = FALSE)
    htmltools::div(
        id = legend_id, class = "bc-legend mt-1",
        `data-amber` = "bg-warning text-dark", `data-neutral` = "bg-secondary",
        warning_note(id, n_noun(length(marked), paste("marked", noun)),
                     class = if (amber) "bg-warning text-dark" else "bg-secondary",
                     items = items)
    )
}

# ── Samples left out of the climate analyses ──────────────────────────────────

#' Samples the climate analyses leave out: no coordinates (coord_missing_summary.tsv) or
#' no climate value at their coordinates (climate_na_excluded.tsv).
#' @param coord_missing data.table(n_missing, missing_samples), or NULL.
#' @param na_excluded   data.table(sample, reason, distance_km), or NULL.
#' @return NULL, or list(no_coords = character, no_climate = data.table(sample, reason)).
#' @noRd
climate_dropped_samples <- function(coord_missing, na_excluded) {
    no_coords <- if (!is.null(coord_missing) && nrow(coord_missing) > 0 &&
                     "missing_samples" %in% names(coord_missing))
        parse_csv_list(coord_missing$missing_samples[1]) else character(0)
    no_climate <- if (!is.null(na_excluded) && nrow(na_excluded) > 0 &&
                      "sample" %in% names(na_excluded))
        data.table::data.table(
            sample = as.character(na_excluded$sample),
            reason = if ("reason" %in% names(na_excluded))
                as.character(na_excluded$reason) else "no climate value")
        else data.table::data.table(sample = character(), reason = character())
    if (length(no_coords) + nrow(no_climate) == 0) return(NULL)
    list(no_coords = no_coords, no_climate = no_climate)
}

#' Tab warning for climate_dropped_samples().
#' @noRd
climate_dropped_note <- function(dropped) {
    if (is.null(dropped)) return(NULL)
    n <- length(dropped$no_coords) + nrow(dropped$no_climate)
    # one item per reason, its samples named
    by_reason <- split(dropped$no_climate$sample, dropped$no_climate$reason)
    warning_note("warn_climate_dropped_samples",
        paste(n_noun(n, "sample"), "left out of climate analyses"),
        items = c(if (length(dropped$no_coords) > 0)
                      paste0("no coordinates: ", name_list(dropped$no_coords)),
                  vapply(names(by_reason), function(r) paste0(r, ": ", name_list(by_reason[[r]])),
                         "", USE.NAMES = FALSE)))
}

# ── Traits: identical and invariant (GWAS chips) ──────────────────────────────

#' Identical traits, over the samples that carry both.
#'
#' 1 - |r| < 1e-9 over the pairwise-complete samples, read only when those hold at least
#' three distinct value pairs (fewer make any two traits perfectly correlated). A
#' constant or all-missing trait is trait_invariant.tsv's and is left to it.
#'
#' @param meta   Processing/tables/metadata.tsv.
#' @param traits Trait columns (metadata columns 5+).
#' @return data.table(predictor, flag = "identical", detail): the trait (named predictor,
#'   as in predictor_flags(), so the chips read both alike) and its twins.
#' @noRd
trait_flags <- function(meta, traits) {
    empty <- data.table::data.table(predictor = character(), flag = character(),
                                    detail = character())
    traits <- intersect(traits, names(meta))
    traits <- traits[vapply(traits, function(t) is.numeric(meta[[t]]), logical(1))]
    if (length(traits) < 2) return(empty)
    twins <- list()
    for (i in seq_along(traits)) for (j in seq_along(traits)) {
        if (j <= i) next
        x <- meta[[traits[i]]]; y <- meta[[traits[j]]]
        ok <- !is.na(x) & !is.na(y)
        if (nrow(unique(data.frame(x = x[ok], y = y[ok]))) < 3) next
        r <- suppressWarnings(stats::cor(x[ok], y[ok]))
        if (is.na(r) || (1 - abs(r)) >= 1e-9) next
        twins[[traits[i]]] <- c(twins[[traits[i]]], traits[j])
        twins[[traits[j]]] <- c(twins[[traits[j]]], traits[i])
    }
    if (length(twins) == 0) return(empty)
    data.table::data.table(predictor = names(twins), flag = "identical",
                           detail = vapply(twins, paste, "", collapse = ", "))
}

# ── Cross-entropy without a minimum ───────────────────────────────────────────

#' The sNMF cross-entropy curve has no minimum in the tested range: the lowest best-run
#' cross-entropy (the point plotted per K) is at the largest K tested, over three or more
#' K values (with two, a minimum inside the range cannot show).
#' @param ce data.table(K, run, cross_entropy, best_run) (cross_entropy_K{s}-{e}.tsv).
#' @return NULL, or list(k_end, k_values, ce) of the best runs.
#' @noRd
cross_entropy_no_minimum <- function(ce) {
    if (is.null(ce) || nrow(ce) == 0 ||
        !all(c("K", "cross_entropy", "best_run") %in% names(ce))) return(NULL)
    best <- ce[as.logical(ce$best_run) %in% TRUE, ]
    best <- best[order(best$K), ]
    if (data.table::uniqueN(best$K) < 3 || anyNA(best$cross_entropy)) return(NULL)
    if (best$K[which.min(best$cross_entropy)] != max(best$K)) return(NULL)
    list(k_end = max(best$K), k_values = best$K, ce = best$cross_entropy)
}

#' Warning for cross_entropy_no_minimum().
#' @noRd
cross_entropy_note <- function(nm) {
    if (is.null(nm)) return(NULL)
    # short: the cross-entropy card's header, which names the curve, is narrow
    warning_note("warn_cross_entropy",
        paste0("No minimum up to K = ", nm$k_end),
        items = paste0("K = ", nm$k_values, ": ", sprintf("%.4f", nm$ce),
                       ifelse(nm$k_values == nm$k_end, " (lowest)", "")))
}

# ── dbMEM: the fallback axis ──────────────────────────────────────────────────

#' How the geography term's dbMEM axes were chosen (dbmem_selected.tsv, `selection`).
#' @param sel data.table(mem, selected[, selection]).
#' @return NULL (forward selection, or nothing to say), "fallback" with the axis, or
#'   "unknown" for a table written before the column (2026-10-30): never read as "no
#'   fallback".
#' @noRd
dbmem_selection_state <- function(sel) {
    if (is.null(sel) || nrow(sel) == 0 || !"mem" %in% names(sel)) return(NULL)
    if (!"selection" %in% names(sel)) return(list(state = "unknown"))
    fb <- sel$mem[sel$selection %in% "fallback_best_marginal_r2"]
    if (length(fb) == 0) return(NULL)
    list(state = "fallback", mems = fb)
}

#' Warning for dbmem_selection_state(): amber for the fallback, neutral for unknown.
#' @noRd
dbmem_selection_note <- function(st) {
    if (is.null(st)) return(NULL)
    if (identical(st$state, "unknown")) {
        return(warning_note("warn_dbmem_fallback", "dbMEM selection unknown",
            items = "dbmem_selected.tsv has no selection column (written before 2026-10-30)",
            class = "bg-secondary"))
    }
    warning_note("warn_dbmem_fallback", "Geography term is a fallback dbMEM axis",
        items = paste0(paste(st$mems, collapse = ", "),
                       ": highest marginal R²; no axis passed forward selection"))
}
