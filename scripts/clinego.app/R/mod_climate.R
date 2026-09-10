

#' Climate tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_climate_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        render_tsv_table <- function(path) {
            dt <- if (!file_ok(path)) data.table::data.table() else
                tryCatch(data.table::fread(path, sep = "\t", header = TRUE),
                         error = function(e) data.table::data.table())
            safe_datatable(as.data.frame(dt),
                           options = list(dom = "t", scrollX = TRUE, pageLength = -1),
                           rownames = FALSE)
        }

        # ── Predictors ───────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run mode=climate.")

            # Collinearity-drop badge, next to the standard help_note badge —
            # the old rda_predictor_collinearity.png is retired (redundant
            # with this heatmap); which predictors the |r| pre-screen dropped
            # (PreGEA.Advanced.collinearity_r) is still on disk in
            # rda_predictor_collinearity.tsv and shown here instead.
            collin_note <- shiny::reactive({
                path <- pregea_table_path(pd$name, "rda", "rda_predictor_collinearity")
                if (!file_ok(path)) return(help_note("climate_heatmap"))
                dt <- tryCatch(data.table::fread(path, sep = "\t", header = TRUE),
                               error = function(e) data.table::data.table())
                dropped <- if (nrow(dt) > 0 && "action" %in% names(dt)) dt[dt$action == "dropped", ] else dt[0]
                extra_badge <- if (nrow(dropped) > 0) {
                    filter_note(
                        paste0(nrow(dropped), " dropped (PreGEA)"),
                        htmltools::p(htmltools::strong(nrow(dropped), "predictor(s) dropped for collinearity"),
                            " before the PreGEA RDA fit (PreGEA.Advanced.collinearity_r): ",
                            paste(dropped$predictor, collapse = ", "), ". ",
                            "Never auto-applied to this heatmap — informational only."),
                        class = "bg-warning text-dark"
                    )
                } else NULL
                htmltools::tagList(help_note("climate_heatmap"), extra_badge)
            })

            mod_image_card_server("climate_heatmap",
                path       = shiny::reactive(climate_heatmap_path(pd$name)),
                title      = shiny::reactive("Climate Correlation Heatmap"),
                dl_name    = shiny::reactive("climate_heatmap"),
                suggestion = suggestion,
                note       = collin_note
            )
            mod_image_card_server("climate_density",
                path       = shiny::reactive(climate_density_path(pd$name)),
                title      = shiny::reactive("Climate Density (Present)"),
                dl_name    = shiny::reactive("climate_density_present"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_density"))
            )
            output$climate_invariant_warning <- shiny::renderUI({
                p <- climate_invariant_path(pd$name)
                if (!file.exists(p)) return(NULL)
                dt <- tryCatch(
                    data.table::fread(p, sep = "\t", header = TRUE),
                    error = function(e) data.table::data.table()
                )
                if (nrow(dt) == 0 || !"predictor" %in% names(dt)) return(NULL)
                preds <- dt$predictor
                items <- lapply(preds, function(pr) {
                    reason <- if ("reason" %in% names(dt)) {
                        r <- dt[dt$predictor == pr, ]$reason
                        if (length(r) > 0) paste0(" — ", r[1]) else ""
                    } else ""
                    htmltools::tags$li(htmltools::tags$code(pr), reason)
                })
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    filter_note(
                        paste0(length(preds), " excluded"),
                        htmltools::p(
                            htmltools::tags$strong(length(preds),
                                if (length(preds) == 1) "climate predictor" else "climate predictors",
                                "excluded due to zero variance at sample sites"),
                            " — not used in GEA or Gradient Forest analyses. ",
                            "Consider removing from ", htmltools::tags$code("Climate.predictors"), ".",
                            htmltools::tags$ul(class = "mb-0 mt-1", items)
                        ),
                        class = "bg-warning text-dark"
                    )
                )
            })
        })

        # ── Variance partitioning ────────────────────────────────────────────
        # Headline number for the section: the Venn plot below deliberately
        # does NOT draw "Unexplained" (standard Venn convention — it's
        # everything outside the circles), so without this badge a user has
        # no way to see how much of the response was explained at all. Same
        # data-derived-title badge pattern as relatedness_note()/
        # wza_collapse_note() (utils_ui.R) — the number itself is the visible
        # label, full context is the hover body.
        output$variance_explained_badge <- shiny::renderUI({
            pd <- project_data()
            p <- climate_table_path(pd$name, "varpart", "variance_partition")
            if (!file_ok(p)) return(NULL)
            dt <- tryCatch(data.table::fread(p, sep = "\t", header = TRUE),
                           error = function(e) data.table::data.table())
            if (nrow(dt) == 0 || !all(c("component", "variance_pct") %in% names(dt))) return(NULL)
            unexp <- dt[dt$component == "Unexplained", ]
            if (nrow(unexp) == 0) return(NULL)
            unexplained_pct <- unexp$variance_pct[1]
            explained_pct   <- 100 - unexplained_pct
            model <- if ("model" %in% names(dt) && nrow(dt) > 0) dt$model[1] else NA_character_

            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    sprintf("%.1f%% variance explained", explained_pct),
                    htmltools::p(
                        htmltools::strong(sprintf("%.1f%%", explained_pct)),
                        " of the genomic response's variance is jointly explained by",
                        if (!is.na(model)) paste0(" ", model) else " climate/structure/geography",
                        "; the remaining ", sprintf("%.1f%%", unexplained_pct),
                        " is unexplained residual — typically large at this sample size ",
                        "(adjusted R² penalizes many predictors relative to N). ",
                        "The Venn below breaks the explained share down by which factor(s) ",
                        "contribute; Unexplained itself isn't drawn there (it's everything ",
                        "outside the circles, the standard Venn convention) — this badge is ",
                        "where its size is reported."
                    ),
                    class = "bg-secondary"
                )
            )
        })

        # dbMEM has no config knob for the too-few-sites case (always
        # site-level, skip-with-warning below 3 sites — see
        # scripts/pregea_dbmem.R). Surface that skip as a visible badge here;
        # dbmem_diagnostics.tsv's `status` otherwise sits unseen inside a
        # collapsed accordion table (CLAUDE.md Rule 8 — advisory text belongs
        # in the Shiny note, not left undiscoverable).
        output$dbmem_skip_warning <- shiny::renderUI({
            pd <- project_data()
            p <- climate_table_path(pd$name, "spatial", "dbmem_diagnostics")
            if (!file_ok(p)) return(NULL)
            dt <- tryCatch(data.table::fread(p, sep = "\t", header = TRUE),
                           error = function(e) data.table::data.table())
            if (nrow(dt) == 0 || !all(c("key", "value") %in% names(dt))) return(NULL)
            dv <- setNames(dt$value, dt$key)
            status <- dv[["status"]]
            if (is.null(status) || identical(status, "ok")) return(NULL)

            n_sites <- dv[["n_sites"]] %||% "?"
            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    "dbMEM skipped",
                    htmltools::p(
                        htmltools::strong("dbMEM skipped"), " (", status, ") — only ",
                        n_sites, " site(s) with coordinates. ",
                        "The geography fraction was not computed; variance partitioning ",
                        "fell back to climate vs. structure. Expand sampling to ≥ 3 ",
                        "distinct sites, or run Gradient Forest without spatial correction (",
                        htmltools::code("Maladaptation.methods.gradient_forest.spatial_correction: without"),
                        ")."
                    ),
                    class = "bg-warning text-dark"
                )
            )
        })

        # Confounding check is its own dedicated 2-table climate-vs-geography
        # fit, always computed when geography is available, independent of
        # which variance-partition tree (2-way/3-way) the donut below shows
        # (docs/rda_research.md A.2 point 4) — surfaced as a small badge, not
        # mixed into the variance-partition table (CLAUDE.md Rule 8).
        output$confounding_badge <- shiny::renderUI({
            pd <- project_data()
            p <- climate_table_path(pd$name, "varpart", "climate_confounding")
            if (!file_ok(p)) return(NULL)
            dt <- tryCatch(data.table::fread(p, sep = "\t", header = TRUE),
                           error = function(e) data.table::data.table())
            if (nrow(dt) == 0 || !"confounded" %in% names(dt)) return(NULL)
            confounded <- isTRUE(as.logical(dt$confounded[1]))
            shared_pct <- dt$shared_pct[1]; max_unique_pct <- dt$max_unique_pct[1]
            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    if (confounded) "Climate/Geography: confounded" else "Climate/Geography: not confounded",
                    htmltools::p(
                        htmltools::strong("Climate – geography overlap check"), " (independent of the ",
                        "structure covariate, always climate-vs-geography only): the shared fraction (",
                        sprintf("%.1f%%", shared_pct), ") ",
                        if (confounded) "EXCEEDS" else "does not exceed",
                        " the larger of the two unique fractions (", sprintf("%.1f%%", max_unique_pct), "). ",
                        if (confounded) "High confounding means climate and geography are hard to tell apart in this dataset — interpret unique fractions cautiously." else "Climate and geography contribute distinguishable signal here."
                    ),
                    class = if (confounded) "bg-warning text-dark" else "bg-secondary"
                )
            )
        })

        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run mode=climate.")

            mod_image_card_server("dbmem_screeplot",
                path       = shiny::reactive(climate_plot_path(pd$name, "spatial", "dbmem_screeplot")),
                title      = shiny::reactive("dbMEM Screeplot"),
                dl_name    = shiny::reactive("dbmem_screeplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_dbmem_screeplot"))
            )
            mod_image_card_server("dbmem_selection_path",
                path       = shiny::reactive(climate_plot_path(pd$name, "varpart", "dbmem_selection_path")),
                title      = shiny::reactive("dbMEM Forward-selection Path"),
                dl_name    = shiny::reactive("dbmem_selection_path"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_dbmem_selection_path"))
            )
            mod_image_card_server("varpart_venn",
                path       = shiny::reactive(climate_plot_path(pd$name, "varpart", "varpart_venn")),
                title      = shiny::reactive("Variance Partition"),
                dl_name    = shiny::reactive("variance_partition"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_varpart_venn"))
            )
            mod_image_card_server("px_barplot",
                path       = shiny::reactive(climate_plot_path(pd$name, "varpart", "px_barplot")),
                title      = shiny::reactive("Px per Predictor (Lasky et al. 2012)"),
                dl_name    = shiny::reactive("px_barplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_px_barplot"))
            )
        })
        output$varpart_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "varpart", "variance_partition"))
        })
        output$px_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "varpart", "px_per_variable"))
        })
        output$dbmem_diagnostics_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "spatial", "dbmem_diagnostics"))
        })
    })
}

#' climate tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_climate_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "climate",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: both headline figures, in focus ─────────────────────────
            htmltools::div(
                class = "lab-hero-col",

                shiny::uiOutput(ns("climate_invariant_warning")),
                lab_hero(mod_image_card_ui(ns("climate_heatmap"))),

                htmltools::div(
                    class = "lab-varpart-badges",
                    shiny::uiOutput(ns("variance_explained_badge")),
                    shiny::uiOutput(ns("dbmem_skip_warning")),
                    shiny::uiOutput(ns("confounding_badge"))
                ),
                lab_hero(mod_image_card_ui(ns("varpart_venn")))
            ),

            # ── RIGHT: plenty of plots ────────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",

                lab_section_header("Predictor distributions", icon = "graph-up"),
                # density_plot_present is a multi-panel grid (5250x3600, 18.9 MP).
                # A full-width row keeps the per-predictor panels legible; a
                # standard tile height caps the width too and collapses it.
                htmltools::div(
                    class = "lab-density",
                    lab_thumb_grid(
                        cols = 1,
                        lab_thumb(mod_image_card_ui(ns("climate_density")))
                    )
                ),

                lab_section_header("Variance partitioning & spatial structure",
                                   icon = "pie-chart"),
                # Three ~1.40 plots -> three columns, matching their aspect ratio.
                lab_thumb_grid(
                    cols = 3,
                    lab_thumb(mod_image_card_ui(ns("px_barplot"))),
                    lab_thumb(mod_image_card_ui(ns("dbmem_screeplot"))),
                    lab_thumb(mod_image_card_ui(ns("dbmem_selection_path")))
                )
            )
        ),

        # Below the fold: the three varpart/dbMEM tables. Collapsed -> suspended.
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE, multiple = TRUE,
            bslib::accordion_panel("Variance partition", value = "vp",
                DT::DTOutput(ns("varpart_table"))),
            bslib::accordion_panel("Px per variable", value = "px",
                DT::DTOutput(ns("px_table"))),
            bslib::accordion_panel("dbMEM diagnostics", value = "dbmem",
                DT::DTOutput(ns("dbmem_diagnostics_table")))
        )
    )
}

