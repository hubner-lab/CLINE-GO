#' Home tab UI
#'
#' Pipeline overview dashboard: value boxes, module status, summary table.
#'
#' @param id module namespace id
#' @noRd
mod_home_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("home_content"))
}

#' Read-only regime badge
#'
#' The regime is declared at project creation and is not editable afterwards, so
#' this is a statement of fact rather than a control.
#' @noRd
.regime_badge <- function(regime) {
    gwas_only <- identical(regime, "gwas_only")
    htmltools::div(
        class = "mb-2",
        htmltools::span(
            class = paste("badge", if (gwas_only) "bg-secondary" else "bg-primary"),
            if (gwas_only) "GWAS-only project" else "Standard project"
        ),
        htmltools::span(
            class = "text-muted small ms-2",
            if (gwas_only)
                "No coordinates — climate, GEA and maladaptation are not part of this project."
            else
                "Geography and climate — all modules available."
        )
    )
}

.mod_home_dashboard_ui <- function(ns, regime = "standard") {
    htmltools::tagList(
        .regime_badge(regime),

        # Value box row: compact, like every tab's (I7)
        lab_kpi_tiles(
            lab_value_box("Samples", shiny::textOutput(ns("n_samples")), "primary", "people"),
            lab_value_box("SNPs (filtered)", shiny::textOutput(ns("n_snps")), "info", "database"),
            lab_value_box("K best", shiny::textOutput(ns("k_best")), "success", "diagram-3"),
            lab_value_box("Assoc. methods", shiny::textOutput(ns("n_methods")), "warning",
                          "bar-chart"),
            lab_value_box("Regions", shiny::textOutput(ns("n_regions")), "danger", "geo-alt")
        ),

        # Pipeline status card — click a completed module to see its summary
        bslib::card(
            full_screen = TRUE,
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span(
                    bsicons::bs_icon("check-circle"),
                    " Pipeline Status"
                ),
                lab_download_button(ns("dl_summary"), "CSV")
            ),
            bslib::card_body(
                shiny::uiOutput(ns("module_status"))
            )
        )
    )
}

#' Home tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param has_project  reactive: TRUE once a project is selected
#' @noRd
mod_home_server <- function(id, project_data, has_project) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        summary_data <- shiny::reactive({
            pd <- project_data()
            load_pipeline_summary(pd$name)
        })

        # ── Dashboard vs Getting Started ───────────────────────────────────────
        # Declared at project creation; drives which modules this tab talks about.
        regime <- shiny::reactive(config_regime(project_data()$config))

        output$home_content <- shiny::renderUI({
            # with no project every tab would stay blank (I11)
            if (!has_project()) {
                return(plot_placeholder("No project yet",
                    "Choose \u201c+ New Project\u201d in the project list at the top left."))
            }
            pd <- project_data()
            st <- check_module_status(pd$name, summary_data())
            gwas_only <- identical(regime(), "gwas_only")
            if (!any(st)) {
                # Empty project — show getting started guide
                bslib::card(
                    bslib::card_header(
                        bsicons::bs_icon("rocket-takeoff"),
                        " Getting Started"
                    ),
                    bslib::card_body(
                        .regime_badge(regime()),
                        htmltools::tags$p(
                            class = "text-muted mb-3",
                            "Your project is set up. Follow these steps to run the pipeline:"
                        ),
                        htmltools::tags$ol(
                            class = "getting-started-steps",
                            htmltools::tags$li(
                                htmltools::tags$strong("Review your config"),
                                " \u2014 Open the sidebar on the left of any tab to check input paths and parameters; Save, top right, writes them to the project's config."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Run Processing"),
                                " \u2014 Click \u2018Run Processing\u2019 in the ", htmltools::tags$strong("Processing"), " sidebar. It filters your VCF, computes sample QC and prepares intermediate files; its tab shows the QC plots."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Run PreStructure"),
                                " \u2014 In the ", htmltools::tags$strong("PreStructure"), " tab, run sNMF across the K range (sNMF.k_start to sNMF.k_end) and read the cross-entropy plot and the PCA."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Run Structure"),
                                if (gwas_only)
                                    " \u2014 Set sNMF.k_best in the Structure sidebar, then run Structure."
                                else
                                    " \u2014 Set sNMF.k_best in the Structure sidebar, then run Structure: the ancestry piemaps and the climate extraction."
                            ),
                            if (gwas_only) htmltools::tags$li(
                                htmltools::tags$strong("Run Phenotypic, then GWAS"),
                                " \u2014 Characterize the traits in the ", htmltools::tags$strong("Phenotypic"), " tab, then pick traits and methods in the ", htmltools::tags$strong("GWAS"), " sidebar and run."
                            ) else htmltools::tags$li(
                                htmltools::tags$strong("Run Environmental and Phenotypic, then GEA"),
                                " \u2014 Characterize the predictors and traits in the ", htmltools::tags$strong("Environmental"), " and ", htmltools::tags$strong("Phenotypic"), " tabs, then pick predictors and methods in the ", htmltools::tags$strong("GEA"), " sidebar and run."
                            ),
                            if (gwas_only) htmltools::tags$li(
                                htmltools::tags$strong("Explore results"),
                                " \u2014 Drill into significant regions, genes, GO enrichment and haplotypes from the GWAS Manhattan plots."
                            ) else htmltools::tags$li(
                                htmltools::tags$strong("Explore results"),
                                " \u2014 Continue with GWAS, GEAxGWAS and Maladaptation as your research requires; regions, genes, GO enrichment and haplotypes open from any Manhattan plot."
                            )
                        )
                    )
                )
            } else {
                # Normal dashboard
                .mod_home_dashboard_ui(ns, regime())
            }
        })

        # ── Value boxes ────────────────────────────────────────────────────────
        # pipeline_summary.tsv's value for (step, metric); call it inside a render. A
        # data.frame on purpose: inside a data.table's [ ], `step` and `metric` name its
        # columns, so the filter compared each column with itself and every box showed
        # the first row's value (bffd38).
        summary_val <- function(step, metric, default = "—") {
            s <- as.data.frame(summary_data())
            if (nrow(s) == 0) return(default)
            row <- s[s$step == step & s$metric == metric, , drop = FALSE]
            if (nrow(row) == 0) default else as.character(row$value[1])
        }
        # The association module whose results this regime's Home reports
        assoc_step <- shiny::reactive(if (identical(regime(), "gwas_only")) "gwas" else "gea")

        output$n_samples <- shiny::renderText(summary_val("processing", "samples_after_filtering"))
        output$n_snps    <- shiny::renderText(summary_val("processing", "snps_after_filtering"))
        output$k_best    <- shiny::renderText({
            pd <- project_data()
            if (!is.na(pd$k_best)) as.character(pd$k_best) else "—"
        })
        output$n_methods <- shiny::renderText({
            pd <- project_data()
            methods <- find_assoc_methods(pd$name,
                                          if (identical(assoc_step(), "gwas")) MOD_GWAS else MOD_GEA)
            if (length(methods) > 0) as.character(length(methods)) else "—"
        })
        output$n_regions <- shiny::renderText(summary_val(assoc_step(), "regions_combined"))

        # ── Module status — clickable rows ────────────────────────────────────
        # Every module but Home, in the module bar's order and names (N3: the
        # Environmental row was missing); a module's pipeline_summary.tsv step is its id.
        reg <- module_registry()
        reg <- reg[reg$id != "home", , drop = FALSE]
        module_labels <- stats::setNames(reg$label, reg$id)

        selected_module <- shiny::reactiveVal(NULL)

        output$module_status <- shiny::renderUI({
            st  <- check_module_status(project_data()$name, summary_data())
            sel <- selected_module()
            s   <- summary_data()

            # Out-of-regime modules are omitted, not shown as "—": a dash reads as
            # "not run yet", which is wrong when the module does not exist here.
            shown <- intersect(names(module_labels), module_ids_for_regime(regime()))

            rows <- lapply(shown, function(nm) {
                ok     <- isTRUE(st[[nm]])
                active <- identical(sel, nm)

                row <- if (ok) {
                    htmltools::div(
                        class = paste("mode-status mode-status-clickable", if (active) "mode-status-active" else ""),
                        onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("module_click"), nm),
                        htmltools::span(class = "status-icon status-ok", bsicons::bs_icon("check-lg")),
                        htmltools::span(module_labels[[nm]]),
                        htmltools::span(
                            class = paste("ms-auto status-chevron", if (!active) "text-muted" else ""),
                            bsicons::bs_icon(if (active) "chevron-down" else "chevron-right")
                        )
                    )
                } else {
                    htmltools::div(
                        class = "mode-status",
                        htmltools::span(class = "status-icon status-missing", "\u2014"),
                        htmltools::span(module_labels[[nm]])
                    )
                }

                # Inline detail panel immediately after this row when active
                detail <- if (active) {
                    metrics   <- if (identical(nm, "gea_x_gwas")) {
                        # no summary rows: the outputs that make it done
                        files <- geaxgwas_outputs(project_data()$name)
                        data.frame(metric = rep("output", length(files)), value = files)
                    } else if (nrow(s) > 0) {
                        s[s$step == nm, c("metric", "value")]
                    } else data.frame(metric = character(), value = character())
                    tbl_rows  <- lapply(seq_len(nrow(metrics)), function(i) {
                        htmltools::tags$tr(
                            htmltools::tags$td(class = "detail-metric", as.character(metrics$metric[i])),
                            htmltools::tags$td(class = "detail-value",  as.character(metrics$value[i]))
                        )
                    })
                    htmltools::div(
                        class = "module-detail-panel",
                        htmltools::tags$table(
                            class = "module-detail-table",
                            htmltools::tags$tbody(tbl_rows)
                        )
                    )
                } else NULL

                htmltools::tagList(row, detail)
            })
            htmltools::tagList(rows)
        })

        # Toggle selected module on click
        shiny::observeEvent(input$module_click, {
            nm <- input$module_click
            if (identical(selected_module(), nm)) selected_module(NULL) else selected_module(nm)
        })

        output$dl_summary <- shiny::downloadHandler(
            filename = function() {
                paste0("pipeline_summary_", project_data()$name, ".csv")
            },
            content = function(file) {
                s <- summary_data()
                utils::write.csv(s, file, row.names = FALSE)
            }
        )
    })
}
