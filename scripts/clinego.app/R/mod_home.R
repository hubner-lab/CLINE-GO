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

        # Value box row
        bslib::layout_column_wrap(
            width = "200px",
            fill  = FALSE,
            bslib::value_box(
                title = "Samples",
                value = shiny::textOutput(ns("n_samples")),
                theme = "primary",
                showcase = bsicons::bs_icon("people")
            ),
            bslib::value_box(
                title = "SNPs (filtered)",
                value = shiny::textOutput(ns("n_snps")),
                theme = "info",
                showcase = bsicons::bs_icon("database")
            ),
            bslib::value_box(
                title = "K best",
                value = shiny::textOutput(ns("k_best")),
                theme = "success",
                showcase = bsicons::bs_icon("diagram-3")
            ),
            bslib::value_box(
                title = "Assoc. methods",
                value = shiny::textOutput(ns("n_methods")),
                theme = "warning",
                showcase = bsicons::bs_icon("bar-chart")
            ),
            bslib::value_box(
                title = "Regions",
                value = shiny::textOutput(ns("n_regions")),
                theme = "danger",
                showcase = bsicons::bs_icon("geo-alt")
            )
        ),

        shiny::br(),

        # Pipeline status card — click a completed module to see its summary
        bslib::card(
            full_screen = TRUE,
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span(
                    bsicons::bs_icon("check-circle"),
                    " Pipeline Status"
                ),
                shiny::downloadButton(ns("dl_summary"), "CSV",
                                      class = "btn-sm btn-outline-secondary")
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
#' @noRd
mod_home_server <- function(id, project_data) {
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
                                " \u2014 Open the sidebar on the left to check input paths and filtering parameters. Adjust MAF, missingness thresholds, and K range as needed."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Run Processing"),
                                " \u2014 Click \u2018Run Processing\u2019 in the sidebar. This filters your VCF, computes sample QC, and prepares intermediate files."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Check QC"),
                                " \u2014 Switch to the ", htmltools::tags$strong("Processing"), " tab to review QC plots (missingness, MAF, SNP density, filtering attrition)."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Run Structure"),
                                " \u2014 Go to the ", htmltools::tags$strong("Structure"), " tab and run sNMF to estimate population structure across K values."
                            ),
                            htmltools::tags$li(
                                htmltools::tags$strong("Select K"),
                                if (gwas_only)
                                    " \u2014 Review the cross-entropy plot, set snmf.k_best in the sidebar, then run Structure K."
                                else
                                    " \u2014 Review the cross-entropy plot, set snmf.k_best in the sidebar, then run Structure K to generate piemaps."
                            ),
                            if (gwas_only) htmltools::tags$li(
                                htmltools::tags$strong("Run GWAS"),
                                " \u2014 Pick traits and association methods in the ", htmltools::tags$strong("GWAS"), " tab, then run."
                            ) else htmltools::tags$li(
                                htmltools::tags$strong("Run Association"),
                                " \u2014 Configure climate predictors and association methods in the ", htmltools::tags$strong("Association"), " tab, then run."
                            ),
                            if (gwas_only) htmltools::tags$li(
                                htmltools::tags$strong("Explore results"),
                                " \u2014 Drill into significant regions, genes, GO enrichment and haplotypes from the GWAS Manhattan plots."
                            ) else htmltools::tags$li(
                                htmltools::tags$strong("Explore results"),
                                " \u2014 Continue with Phenotype Association, Overlapping Regions, Maladaptation, and Haplotype Analysis as your research requires."
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
        summary_val <- function(step, metric, default = "—") {
            shiny::reactive({
                s <- summary_data()
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            })
        }

        output$n_samples <- shiny::renderText(summary_val("processing", "n_samples")())
        output$n_snps    <- shiny::renderText(summary_val("processing", "n_snps_filtered")())
        output$k_best    <- shiny::renderText({
            pd <- project_data()
            if (!is.na(pd$k_best)) as.character(pd$k_best) else "—"
        })
        output$n_methods <- shiny::renderText({
            pd <- project_data()
            methods <- find_assoc_methods(pd$name)
            if (length(methods) > 0) as.character(length(methods)) else "—"
        })
        output$n_regions <- shiny::renderText(
            summary_val("gea", "n_regions_combined")()
        )

        # ── Module status — clickable rows ────────────────────────────────────
        module_labels <- c(
            processing    = "Processing",
            prestructure  = "PreStructure",
            structure     = "Structure",
            traits        = "Phenotypic",
            pregea        = "PreGEA",
            gea           = "GEA",
            gwas          = "GWAS",
            gea_x_gwas    = "GEAxGWAS",
            maladaptation = "Maladaptation"
        )
        # Step name used in pipeline_summary.tsv per module
        module_steps <- c(
            processing    = "processing",
            prestructure  = "prestructure",
            structure     = "structure",
            traits        = "traits",
            pregea        = "pregea",
            gea           = "gea",
            gwas          = "gwas",
            gea_x_gwas    = "gea_x_gwas",
            maladaptation = "maladaptation"
        )

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
                        htmltools::span(class = "status-icon status-ok", "\u2713"),
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
                    step      <- module_steps[[nm]]
                    metrics   <- if (nrow(s) > 0) s[s$step == step, c("metric", "value")] else data.frame(metric = character(), value = character())
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
