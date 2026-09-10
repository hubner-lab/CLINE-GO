# fct_model_compare.R
# Helpers for on-demand model comparison (compare_offsets.R backend).
#
# A "model key" = "method:::suffix" (triple-colon separator).
# Two-model comparison cache lives under _intermediate/model_compare/{comparison_key}/.
# Novelty cache (ExDet) lives under _intermediate/novelty/{scenario}/ and is shared
# across all comparisons for the same future scenario.

# ─── Discovery ────────────────────────────────────────────────────────────────

#' Return a named vector of all available two-model comparison keys.
#' Each name is "Model A label vs Model B label"; each value is a sorted comparison key.
#' @noRd
find_model_comparisons <- function(project) {
    all_models <- find_mala_models(project)
    if (length(all_models) < 2) return(character(0))
    combos <- utils::combn(seq_along(all_models), 2, simplify = FALSE)
    vals  <- vapply(combos, function(i) {
        model_compare_key(all_models[i[1]], all_models[i[2]])
    }, character(1))
    nms   <- vapply(combos, function(i) {
        paste0(names(all_models)[i[1]], " vs ", names(all_models)[i[2]])
    }, character(1))
    setNames(vals, nms)
}

# ─── Stats loader ─────────────────────────────────────────────────────────────

#' Load stats.json for a comparison (NULL if missing/corrupt)
#' @noRd
load_compare_stats <- function(project, key_a, key_b) {
    p <- model_compare_stats_path(project, key_a, key_b)
    if (!file_ok(p)) return(NULL)
    tryCatch(
        jsonlite::read_json(p, simplifyVector = TRUE),
        error = function(e) NULL
    )
}

# ─── On-demand compute trigger ────────────────────────────────────────────────

#' Resolve method/suffix from a "method:::suffix" key to pipeline file paths.
#' Returns list(site_tsv, map_tsv) or NULL if files missing.
#' @noRd
model_key_to_files <- function(project, model_key) {
    p <- parse_model_key(model_key)
    if (is.null(p)) return(NULL)
    site <- gf_site_table_path(project, p$suffix, method = p$method)
    map  <- gf_map_table_path(project,  p$suffix, method = p$method)
    if (!file_ok(site) || !file_ok(map)) return(NULL)
    list(site_tsv = site, map_tsv = map, method = p$method, suffix = p$suffix)
}

#' Build the compare_offsets.R shell call and run it.
#' Invokes Rscript synchronously (blocks until done).
#' Returns TRUE on success, FALSE otherwise.
#' @noRd
run_compare_offsets <- function(project, key_a, key_b,
                                env_site_pres, env_all_pres, env_all_fut,
                                predictors, raster_tif,
                                top_k = 10L, scenario_label = "future") {
    files_a <- model_key_to_files(project, key_a)
    files_b <- model_key_to_files(project, key_b)
    if (is.null(files_a) || is.null(files_b)) {
        message("ERROR: compare: missing input files for ", key_a, " or ", key_b)
        return(FALSE)
    }

    label_a <- names(find_mala_models(project))[
        match(key_a, find_mala_models(project))]
    label_b <- names(find_mala_models(project))[
        match(key_b, find_mala_models(project))]
    if (is.na(label_a) || is.na(label_b)) {
        label_a <- key_a; label_b <- key_b
    }

    cache_dir   <- model_compare_dir(project, key_a, key_b)
    novelty_dir <- novelty_cache_dir(project, scenario_label)

    pred_str <- paste(predictors, collapse = ",")

    cmd <- paste(
        "Rscript /pipeline/scripts/compare_offsets.R",
        shQuote(cache_dir),
        shQuote(novelty_dir),
        shQuote(files_a$site_tsv),
        shQuote(files_a$map_tsv),
        shQuote(label_a),
        shQuote(files_b$site_tsv),
        shQuote(files_b$map_tsv),
        shQuote(label_b),
        shQuote(env_site_pres),
        shQuote(env_all_pres),
        shQuote(env_all_fut),
        shQuote(pred_str),
        shQuote(raster_tif),
        as.integer(top_k),
        "0"  # N_EXTRA = 0 (no additional models in two-model call)
    )

    message("INFO: Running compare_offsets.R ...")
    ret <- system(cmd, ignore.stdout = FALSE)
    ret == 0L
}

# ─── Data loaders ─────────────────────────────────────────────────────────────

#' Load merged site offsets table for a comparison
#' @noRd
load_compare_site_offsets <- function(project, key_a, key_b) {
    p <- model_compare_site_offsets(project, key_a, key_b)
    if (!file_ok(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, colClasses = c(site = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load rank stability table
#' @noRd
load_compare_rank_stability <- function(project, key_a, key_b) {
    p <- model_compare_rank_stability(project, key_a, key_b)
    if (!file_ok(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, colClasses = c(site = "character")),
        error = function(e) data.table::data.table()
    )
}

# ─── UI helpers ───────────────────────────────────────────────────────────────

#' Render a comparison KPI stat box (value + label + optional p-value).
#' @noRd
compare_stat_box <- function(value, label, p_val = NULL, accent = NULL) {
    p_html <- if (!is.null(p_val) && !is.na(p_val)) {
        p_fmt <- if (p_val < 0.001) "p<0.001" else paste0("p=", signif(p_val, 2))
        htmltools::tags$small(class = "text-muted ms-1", p_fmt)
    } else NULL

    value_html <- htmltools::tags$strong(
        style = if (!is.null(accent)) paste0("color:", accent) else "",
        value
    )
    htmltools::div(
        class = "d-inline-block text-center me-3 mb-2",
        htmltools::div(value_html, p_html),
        htmltools::tags$small(class = "d-block text-muted", label)
    )
}

#' Honesty banner — rendered once at top of comparison view.
#' @noRd
compare_honesty_banner <- function() {
    bslib::card(
        class = "border-warning mb-3",
        bslib::card_body(
            class = "py-2",
            htmltools::div(
                class = "d-flex align-items-center gap-2",
                htmltools::tags$span(
                    style = "color:#d97706; font-size:1.1em;",
                    htmltools::HTML("&#9888;")  # warning triangle
                ),
                htmltools::tags$small(
                    class = "text-muted",
                    "No fitness ground truth exists for these offsets. This view characterises ",
                    htmltools::tags$strong("extrapolation behaviour and cross-model concordance"),
                    " — not which model is correct. Conclusions about which model is ",
                    htmltools::tags$em("better"),
                    " require independent validation data (reciprocal transplant or common-garden experiments)."
                )
            )
        )
    )
}

#' Helper to display an "empty" state (no comparison selected yet).
#' @noRd
compare_empty_state <- function(message_text = "Select Model A and Model B above to begin comparison.") {
    htmltools::div(
        class = "text-center py-5 text-muted",
        htmltools::tags$p(message_text),
        htmltools::tags$small(
            "Novelty (ExDet), rank disagreement, and concordance statistics ",
            "will be computed on first view and cached."
        )
    )
}
