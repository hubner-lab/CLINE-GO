#' RDA Diagnostics & Candidates panel
#'
#' Second accordion panel inside mod_gea.R's per-method accordion, shown only
#' when "RDA" is a configured GEA method. mod_gea.R itself stays fully
#' method-agnostic (it derives methods/traits/colors from discovered files);
#' this module is a separate, RDA-specific detail panel — the same pattern as
#' mod_maladaptation.R's per-method comparison tabs.
#'
#' @param id module namespace id
#' @noRd
mod_rda_details_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("panel"))
}

#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param methods reactive character vector of configured GEA methods (from mod_gea.R)
#' @noRd
mod_rda_details_server <- function(id, project_data, methods) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        has_rda <- shiny::reactive(isTRUE("RDA" %in% methods()))

        diagnostics <- shiny::reactive({
            shiny::req(has_rda())
            load_rda_diagnostics(project_data()$name)
        })
        anova_dt <- shiny::reactive({
            shiny::req(has_rda())
            load_rda_anova(project_data()$name)
        })
        candidates_dt <- shiny::reactive({
            shiny::req(has_rda())
            load_rda_candidates(project_data()$name)
        })

        d <- function(key, default = NA) {
            v <- diagnostics()[[key]]
            if (is.null(v) || identical(v, "")) default else v
        }

        # Numeric formatter for small p-values / thresholds — used across the
        # value boxes, provenance card, and the biplot note below.
        fmt_sci <- function(x) {
            v <- suppressWarnings(as.numeric(x))
            if (is.na(v)) "?" else formatC(v, format = "e", digits = 2)
        }

        # ── Value-box summary row ───────────────────────────────────────────
        output$value_boxes <- shiny::renderUI({
            shiny::req(has_rda())
            fit_mode <- d("fit_mode", "unknown")
            bslib::layout_column_wrap(
                width = 1 / 3, heights_equal = "row",
                bslib::value_box(
                    title = "Fit mode", value = fit_mode,
                    theme = if (identical(fit_mode, "pruned")) "danger" else "success",
                    showcase = bsicons::bs_icon("cpu")
                ),
                bslib::value_box(
                    title = "Markers fitted", value = d("n_markers_fitted", "?"),
                    theme = "info", showcase = bsicons::bs_icon("database-fill")
                ),
                bslib::value_box(
                    title = "RDA axes (K)",
                    value = sprintf("%s / %s", d("rda_axes", "?"), d("rda_axes_max", "?")),
                    theme = "info", showcase = bsicons::bs_icon("bounding-box")
                ),
                bslib::value_box(
                    title = "Adj. R2", value = {
                        v <- suppressWarnings(as.numeric(d("adj_r_squared", NA)))
                        if (is.na(v)) "?" else sprintf("%.4f", v)
                    },
                    theme = "info", showcase = bsicons::bs_icon("graph-up")
                ),
                bslib::value_box(
                    title = "GIF lambda", value = {
                        v <- suppressWarnings(as.numeric(d("gif_lambda", NA)))
                        if (is.na(v)) "?" else sprintf("%.3f", v)
                    },
                    theme = if (!is.na(suppressWarnings(as.numeric(d("gif_lambda", NA)))) &&
                               abs(as.numeric(d("gif_lambda")) - 1) > 0.2) "warning" else "success",
                    showcase = bsicons::bs_icon("rulers")
                ),
                bslib::value_box(
                    title = "Candidates",
                    value = d("n_candidates_total", "?"),
                    htmltools::p(class = "small mb-0", d("candidate_rule", "?"),
                                sprintf(" (p < %s)", fmt_sci(d("candidate_threshold", NA)))),
                    theme = if (identical(d("candidate_threshold_status", "ok"), "ok"))
                                "info" else "warning",
                    showcase = bsicons::bs_icon("crosshair")
                )
            )
        })

        # ── Provenance card: full parameter record + applied significance rule ─
        .kv_table <- function(title, keys) {
            keys_present <- keys[vapply(keys, function(k) !is.na(d(k, NA)), logical(1))]
            if (length(keys_present) == 0) return(NULL)
            htmltools::div(
                htmltools::tags$h6(class = "small text-muted mb-1", title),
                htmltools::tags$table(
                    class = "table table-sm mb-2",
                    htmltools::tags$tbody(lapply(keys_present, function(k) htmltools::tags$tr(
                        htmltools::tags$td(class = "text-muted", style = "width:55%", k),
                        htmltools::tags$td(as.character(d(k, "?")))
                    )))
                )
            )
        }

        .ladder_table <- function() {
            rule <- d("candidate_rule", NA)
            rows <- list(
                list(label = "bonf 0.01/m", key = "n_candidates_bonf_0.01", tag = "bonf_0.01",
                    note = "Capblancq & Forester 2021"),
                list(label = "bonf 0.05/m", key = "n_candidates_bonf_0.05", tag = "bonf_0.05",
                    note = "univariate default"),
                list(label = "q < 0.05",    key = "n_candidates_q_0.05",   tag = "qval_0.05", note = ""),
                list(label = "q < 0.10",    key = "n_candidates_q_0.10",   tag = "qval_0.1",   note = "Capblancq et al. 2018")
            )
            applied_matched <- any(vapply(rows, function(r) identical(r$tag, rule), logical(1)))
            htmltools::div(
                htmltools::tags$h6(class = "small text-muted mb-1",
                                   "Sensitivity ladder — what each published rule would give"),
                htmltools::tags$table(
                    class = "table table-sm mb-1",
                    htmltools::tags$tbody(lapply(rows, function(r) {
                        is_applied <- identical(r$tag, rule)
                        htmltools::tags$tr(
                            htmltools::tags$td(r$label),
                            htmltools::tags$td(d(r$key, "?")),
                            htmltools::tags$td(
                                if (is_applied) htmltools::tags$span(class = "badge bg-primary", "APPLIED"),
                                htmltools::span(class = "text-muted small ms-1", r$note)
                            )
                        )
                    }))
                ),
                if (!applied_matched && !is.na(rule)) htmltools::p(
                    class = "small text-muted mb-1",
                    sprintf("Applied rule (%s) is outside this ladder — top/custom rules are not a published RDA convention.", rule)),
                htmltools::p(class = "small text-muted mb-0", d("sensitivity_ladder_note", ""))
            )
        }

        output$provenance_card <- shiny::renderUI({
            shiny::req(has_rda())
            cfg_rule <- resolve_adjust(project_data()$config, "RDA", "GEA")   # utils_helpers.R
            run_rule <- d("candidate_rule", NA)
            drift <- !is.na(run_rule) && !is.null(cfg_rule) && !identical(run_rule, cfg_rule)
            bslib::card(
                class = if (drift) "border-warning" else NULL,
                bslib::card_header(
                    bsicons::bs_icon("fingerprint"), " Run parameters & significance rule",
                    if (drift) htmltools::span(class = "badge bg-warning ms-2", "stale vs config")
                ),
                bslib::card_body(
                    if (drift) htmltools::div(class = "alert alert-warning small py-1",
                        sprintf(paste0(
                            "These artifacts were produced with rule '%s'. GEA.configs now says ",
                            "'%s'. Re-run GEA to regenerate the candidates table, sig_snps file, ",
                            "QQ and background plots."), run_rule, cfg_rule)),
                    bslib::layout_column_wrap(
                        width = 1 / 2,
                        .kv_table("Model & fit", c(
                            "predictors", "n_predictors", "n_samples", "condition_pcs",
                            "structure_proxy", "fit_mode_requested", "fit_mode",
                            "full_fit_max_snps", "full_fit_max_gb",
                            "n_markers_available", "n_markers_fitted",
                            "n_zero_variance_dropped", "adj_r_squared", "max_vif",
                            "vif_flagged_predictors", "aliased_terms",
                            "axes_requested", "axis_alpha", "rda_axes", "rda_axes_max",
                            "K_selection", "K_floored", "permutations", "seed",
                            "k_best", "cpu", "marker_envelope_status")),
                        htmltools::tagList(
                            .kv_table("Applied significance rule", c(
                                "candidate_rule", "candidate_adjust",
                                "candidate_threshold_value", "candidate_threshold",
                                "candidate_threshold_status", "candidate_n_tested",
                                "candidate_n_na_dropped", "n_candidates_total",
                                "n_candidates_written", "gif_lambda",
                                "qvalue_method", "candidate_qvalue_engine",
                                "pval_hist_flatness_chisq")),
                            .ladder_table()
                        )
                    )
                )
            )
        })

        # ── Warning card: pruned fallback and/or above the validated envelope ─
        output$warning_card <- shiny::renderUI({
            shiny::req(has_rda())
            fit_mode <- d("fit_mode", "full")
            envelope <- d("marker_envelope_status", "")
            msgs <- character(0)
            if (identical(fit_mode, "pruned")) msgs <- c(msgs, trimws(d("fallback_warning", "")))
            if (identical(envelope, "ABOVE_VALIDATED")) msgs <- c(msgs, d("marker_envelope_note", ""))
            msgs <- msgs[nzchar(msgs)]
            if (length(msgs) == 0) return(NULL)
            bslib::card(
                class = "border-danger",
                bslib::card_header(class = "bg-danger-subtle",
                                   bsicons::bs_icon("exclamation-triangle-fill"), " RDA run notice"),
                bslib::card_body(
                    lapply(msgs, function(m) htmltools::tags$pre(class = "small mb-2", m))
                )
            )
        })

        # ── ANOVA table (degenerate rows highlighted) ────────────────────────
        output$anova_table <- DT::renderDT({
            shiny::req(has_rda())
            dt <- anova_dt()
            tbl <- safe_datatable(dt, options = list(pageLength = 15, scrollX = TRUE))
            if ("status" %in% names(dt)) {
                tbl <- DT::formatStyle(tbl, "status",
                    backgroundColor = DT::styleEqual(c("degenerate", "ok"),
                                                     c("rgba(229, 62, 62, 0.12)", "")))
            }
            tbl
        })

        # ── Plot cards ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data(); shiny::req(has_rda())
            mod_image_card_server(
                "screeplot",
                path  = shiny::reactive(rda_plot_path(pd$name, MOD_GEA, "screeplot")),
                title = shiny::reactive("RDA Screeplot"),
                dl_name = shiny::reactive("rda_screeplot"),
                suggestion = shiny::reactive("Run GEA with RDA in GEA.configs"),
                note = shiny::reactive(help_note("rda_screeplot", results = sprintf(
                    "K=%s (max %s), condition_pcs=%s", d("rda_axes", "?"), d("rda_axes_max", "?"),
                    d("condition_pcs", "?"))))
            )
            mod_image_card_server(
                "pval_hist",
                path  = shiny::reactive(rda_plot_path(pd$name, MOD_GEA, "pvalue_histogram")),
                title = shiny::reactive("RDA p-value Histogram"),
                dl_name = shiny::reactive("rda_pvalue_histogram"),
                suggestion = shiny::reactive("Run GEA with RDA in GEA.configs"),
                note = shiny::reactive(help_note("rda_pval_hist", results = sprintf(
                    "GIF lambda=%s, qvalue method=%s", d("gif_lambda", "?"), d("qvalue_method", "?"))))
            )
            mod_image_card_server(
                "biplot",
                path  = shiny::reactive(rda_plot_path(pd$name, MOD_GEA, "loadings_biplot")),
                title = shiny::reactive("RDA Loadings Biplot"),
                dl_name = shiny::reactive("rda_loadings_biplot"),
                suggestion = shiny::reactive("Run GEA with RDA in GEA.configs"),
                note = shiny::reactive(help_note("rda_biplot", results = {
                    st <- d("candidate_threshold_status", "ok")
                    if (!identical(st, "ok"))
                        sprintf("No candidates — threshold unavailable (%s) for rule %s",
                               st, d("candidate_rule", "?"))
                    else
                        sprintf("%s candidate(s) at %s (raw p < %s)",
                               d("n_candidates_total", "?"), d("candidate_rule", "?"),
                               fmt_sci(d("candidate_threshold", NA)))
                }))
            )
        })

        # ── Candidates table (server-side paged — can be large at WGS density) ─
        output$candidates_table <- DT::renderDT({
            shiny::req(has_rda())
            safe_datatable(
                candidates_dt(),
                filter = "top", server = TRUE,
                options = list(pageLength = 25, scrollX = TRUE,
                               order = list(list(which(names(candidates_dt()) == "rank") - 1, "asc")))
            )
        }, server = TRUE)

        # ── Candidates card header — states the rule that produced the table ──
        output$candidates_card_header <- shiny::renderUI({
            shiny::req(has_rda())
            sprintf("Candidates at %s — LD-unclumped (see Region Explorer below for clumped regions)",
                   d("candidate_rule", "?"))
        })

        # ── Panel shell — hidden entirely when RDA isn't configured ──────────
        output$panel <- shiny::renderUI({
            if (!has_rda()) return(NULL)
            htmltools::tagList(
                shiny::uiOutput(ns("value_boxes")),
                shiny::uiOutput(ns("provenance_card")),
                shiny::uiOutput(ns("warning_card")),
                bslib::layout_column_wrap(
                    width = 1 / 3,
                    mod_image_card_ui(ns("screeplot")),
                    mod_image_card_ui(ns("pval_hist")),
                    mod_image_card_ui(ns("biplot"))
                ),
                bslib::card(
                    bslib::card_header("ANOVA (full / by-axis / by-margin)"),
                    DT::DTOutput(ns("anova_table"))
                ),
                bslib::card(
                    bslib::card_header(shiny::uiOutput(ns("candidates_card_header"), inline = TRUE)),
                    DT::DTOutput(ns("candidates_table"))
                )
            )
        })
    })
}
