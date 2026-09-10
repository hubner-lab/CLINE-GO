

#' Maladaptation tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_maladaptation_server <- function(id, project_data, snp_sets_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # All maladaptation models across all methods, keyed as "method:::suffix"
        mala_models <- shiny::reactive({
            if (!is.null(snp_sets_trigger)) snp_sets_trigger()
            find_mala_models(project_data()$name)
        })

        output$model_selector <- shiny::renderUI({
            models <- mala_models()
            if (length(models) == 0)
                return(shiny::p(
                    "No maladaptation runs found. Save a SNP set in the GEA tab and run mode=maladaptation.",
                    class = "text-muted small"))
            shiny::selectInput(ns("model_key"), "Model (method / SNP set)",
                               choices = models, selected = models[1])
        })

        selected_model_key <- shiny::reactive({
            k <- input$model_key
            if (is.null(k)) {
                models <- mala_models()
                if (length(models) == 0) return(NULL)
                models[1]
            } else k
        })

        selected_method <- shiny::reactive({
            p <- parse_model_key(selected_model_key())
            if (is.null(p)) "gradient_forest" else p$method
        })

        selected_suffix <- shiny::reactive({
            p <- parse_model_key(selected_model_key())
            if (is.null(p)) return(NULL)
            p$suffix
        })

        # ── Resolve suffix → set name via manifest (NOT string-stripping) ──────
        # Pipeline emits {name}_spatial | {name}_nospatial | {name} (no spatial).
        # Underscores + 'both' mode make suffix-stripping ambiguous.
        selected_set_name <- shiny::reactive({
            suf <- selected_suffix(); if (is.null(suf)) return(NULL)
            if (!is.null(snp_sets_trigger)) snp_sets_trigger()
            set_names <- list_snp_sets(project_data()$name)$name
            hit <- set_names[suf == set_names |
                             suf == paste0(set_names, "_spatial") |
                             suf == paste0(set_names, "_nospatial")]
            if (length(hit) == 0) NULL else hit[which.max(nchar(hit))]
        })

        # ── Manifest entry for the selected set ───────────────────────────────
        selected_set_entry <- shiny::reactive({
            nm <- selected_set_name(); if (is.null(nm)) return(NULL)
            man <- read_snp_sets_manifest(project_data()$name)
            entries <- Filter(function(x) identical(x$name, nm), man)
            if (length(entries) == 0) NULL else entries[[1]]
        })

        # ── View-details popover ───────────────────────────────────────────────
        output$set_details_btn <- shiny::renderUI({
            e <- selected_set_entry(); if (is.null(e)) return(NULL)
            body <- htmltools::tags$ul(
                class = "small mb-0 ps-3",
                htmltools::tags$li(htmltools::strong("Threshold: "),
                    paste0(e$threshold_type %||% "?", " = ", e$threshold_value %||% "?")),
                htmltools::tags$li(htmltools::strong("Strategy: "),  e$strategy %||% "?"),
                htmltools::tags$li(htmltools::strong("Regime: "),    e$regime %||% "?"),
                htmltools::tags$li(htmltools::strong("Source: "),    e$source_module %||% "?"),
                htmltools::tags$li(htmltools::strong("N SNPs: "),    e$n_snps %||% "?"),
                htmltools::tags$li(htmltools::strong("Created: "),   e$created %||% "?")
            )
            bslib::popover(
                shiny::actionButton(ns("set_details"),
                    label = htmltools::tagList(bsicons::bs_icon("info-circle"), " View details"),
                    class = "btn btn-outline-secondary btn-sm ms-2"),
                title = paste("SNP set:", e$name %||% ""),
                body
            )
        })

        # ── Delete button + confirm modal ──────────────────────────────────────
        output$set_delete_btn <- shiny::renderUI({
            if (is.null(selected_set_name())) return(NULL)
            shiny::actionButton(ns("set_delete"),
                label = htmltools::tagList(bsicons::bs_icon("trash3"), " Delete"),
                class = "btn btn-outline-danger btn-sm ms-1")
        })

        shiny::observeEvent(input$set_delete, {
            nm <- selected_set_name(); if (is.null(nm)) return()
            shiny::showModal(shiny::modalDialog(
                title = "Delete SNP set",
                htmltools::p(
                    "Delete set ", htmltools::strong(nm),
                    " and its Gradient Forest results? This cannot be undone."
                ),
                footer = htmltools::tagList(
                    shiny::modalButton("Cancel"),
                    shiny::actionButton(ns("set_delete_confirm"), "Delete",
                                        class = "btn btn-danger")
                ),
                easyClose = TRUE
            ))
        })

        shiny::observeEvent(input$set_delete_confirm, {
            nm <- selected_set_name()
            pd <- project_data()
            if (is.null(nm) || is.null(pd)) { shiny::removeModal(); return() }
            delete_snp_set(pd$name, nm, remove_gf_results = TRUE)
            shiny::removeModal()
            shiny::showNotification(
                sprintf("Deleted SNP set '%s'.", nm),
                type = "message", duration = 4)
            if (!is.null(snp_sets_trigger))
                snp_sets_trigger(snp_sets_trigger() + 1L)
        })

        # ── G1: SNP count badge with denominator ──────────────────────────────
        output$snp_count_badge <- shiny::renderUI({
            suf <- selected_suffix()
            pd  <- project_data()
            if (is.null(suf) || is.null(pd)) return(NULL)
            # Look up via curated set name (manifest-resolved)
            set_nm    <- selected_set_name()
            snps_path <- if (!is.null(set_nm)) snp_set_path(pd$name, set_nm) else NULL
            if (is.null(snps_path) || !file_ok(snps_path)) return(NULL)
            n <- tryCatch({
                nrow(data.table::fread(snps_path, select = 1L))
            }, error = function(e) NULL)
            if (is.null(n)) return(NULL)
            # Get total filtered SNP count for context
            summ <- load_pipeline_summary(pd$name)
            n_total_row <- summ[summ$step == "processing" & summ$metric == "snps_after_filtering", ]
            n_total <- if (nrow(n_total_row) > 0) suppressWarnings(as.integer(n_total_row$value[1])) else NA_integer_
            label <- if (!is.na(n_total) && n_total > 0) {
                pct <- round(100 * n / n_total, 1)
                paste0("\u25c6 ", n, " / ", format(n_total, big.mark = ","), " adaptive SNPs (", pct, "%)")
            } else {
                paste0("\u25c6 ", n, " adaptive SNPs")
            }
            htmltools::tags$span(
                class = "badge bg-secondary ms-2 align-self-center",
                style = "font-size:0.85rem; vertical-align:middle;",
                label
            )
        })

        # ── Imputation-sensitivity badge (Gradient Forest only) ──────────────
        # Gradient Forest is always fitted on an imputed matrix (gradientForest() has no
        # na.action). On a panel with heavy missingness that invites the obvious reviewer
        # question, so mode=maladaptation refits the same SNP set on observed-only site
        # allele frequencies and this badge reports whether the two agree. The other
        # engines get no badge: LEA::genetic.gap requires individual genotypes with
        # mandatory imputation and cannot be refit this way.
        output$sensitivity_badge <- shiny::renderUI({
            suf    <- selected_suffix()
            method <- selected_method()
            pd     <- project_data()
            if (is.null(suf) || is.null(pd) || !identical(method, "gradient_forest")) return(NULL)
            gf_sensitivity_note(load_gf_sensitivity(pd$name, suf, method))
        })

        # ── G2: Method run parameters badge (method-aware) ───────────────────
        output$method_params_badge <- shiny::renderUI({
            suf    <- selected_suffix()
            method <- selected_method()
            pd     <- project_data()
            if (is.null(suf) || is.null(pd)) return(NULL)
            cfg <- pd$config
            if (method == "gradient_forest") {
                ntree   <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                      "ntree", default = 500L)
                cor_thr <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                      "cor_threshold", default = 0.5)
                spatial <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                      "spatial_correction", default = "with")
                rand    <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                      "random_model", default = TRUE)
                badges <- list(
                    config_badge("ntree", ntree),
                    config_badge("cor.thr", cor_thr),
                    config_badge("spatial", spatial),
                    if (isTRUE(rand)) config_badge("random model", "yes") else NULL
                )
            } else if (method == "geometric_offset") {
                k     <- config_get(cfg, "Maladaptation", "methods", "geometric_offset",
                                    "k", default = "(sNMF k_best)")
                scale <- config_get(cfg, "Maladaptation", "methods", "geometric_offset",
                                    "scale", default = TRUE)
                badges <- list(
                    config_badge("K", if (is.null(k) || k == "") "k_best" else as.character(k)),
                    config_badge("scale", if (isTRUE(scale)) "yes" else "no"),
                    config_badge("spatial", "no (nospatial only)")
                )
            } else if (method == "rda_offset") {
                axes    <- config_get(cfg, "Maladaptation", "methods", "rda_offset",
                                      "axes", default = "auto")
                perms   <- config_get(cfg, "Maladaptation", "methods", "rda_offset",
                                      "permutations", default = 999L)
                cond    <- config_get(cfg, "Maladaptation", "methods", "rda_offset",
                                      "condition_pcs", default = 0L)
                badges <- list(
                    config_badge("axes", axes),
                    config_badge("permutations", perms),
                    config_badge("condition PCs", cond),
                    config_badge("spatial", "no (nospatial only)")
                )
            } else {
                badges <- list()
            }
            htmltools::div(
                class = "ms-2 align-self-center",
                do.call(config_badges_bar, badges)
            )
        })

        # ── G4: Future climate scenario badge ─────────────────────────────────
        output$climate_scenario_badge <- shiny::renderUI({
            suf <- selected_suffix()
            pd  <- project_data()
            if (is.null(suf) || is.null(pd)) return(NULL)
            cfg    <- pd$config
            ssp    <- config_get(cfg, "Future", "ssp",    default = NULL)
            year   <- config_get(cfg, "Future", "year",   default = NULL)
            models <- config_get(cfg, "Future", "models", default = NULL)
            if (is.null(ssp) && is.null(year)) return(NULL)
            n_models <- if (is.null(models)) 0L else length(models)
            label_parts <- c(
                if (!is.null(ssp))  paste0("SSP", ssp),
                if (!is.null(year)) as.character(year),
                if (n_models > 0)   paste0(n_models, " GCMs")
            )
            htmltools::div(
                class = "ms-1 align-self-center",
                htmltools::tags$span(
                    class = "badge bg-info text-dark align-self-center",
                    style = "font-size:0.75rem;",
                    paste0("\u2601 Future: ", paste(label_parts, collapse = " "))
                )
            )
        })

        # ── G3: Offset range summary ──────────────────────────────────────────
        output$offset_summary <- shiny::renderUI({
            dt <- site_data()
            if (is.null(dt) || nrow(dt) == 0 || !"genetic_offset" %in% names(dt))
                return(NULL)
            offs <- dt$genetic_offset
            offs <- offs[!is.na(offs)]
            if (length(offs) == 0) return(NULL)
            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    sprintf("%.4f\u2013%.4f", min(offs), max(offs)),
                    htmltools::p(
                        htmltools::tags$strong("Genetic offset range: "),
                        sprintf("%.4f", min(offs)), " \u2013 ", sprintf("%.4f", max(offs)),
                        " (mean: ", sprintf("%.4f", mean(offs)), ", n=", length(offs), " sites)"
                    ),
                    class = "bg-secondary"
                )
            )
        })

        # ── Piemap variant selector (only shows variants that exist) ───────────
        output$piemap_variant_ui <- shiny::renderUI({
            suf    <- shiny::req(selected_suffix())
            method <- selected_method()
            pd  <- project_data()
            choices <- c("Genetic Offset" = "base")
            if (file_ok(gf_offset_piemap_path(pd$name, suf, "tajima_d", method = method)))
                choices <- c(choices, c("Offset \u00d7 Tajima's D" = "tajima_d"))
            if (file_ok(gf_offset_piemap_path(pd$name, suf, "pi_diversity", method = method)))
                choices <- c(choices, c("Offset \u00d7 Pi Diversity" = "pi_diversity"))
            if (length(choices) == 1) return(NULL)  # only base — no selector needed
            shiny::selectInput(ns("piemap_variant"), "Piemap Type",
                               choices = choices, selected = "base")
        })

        # ── Zoom selector (only shows if zoom maps exist) ──────────────────────
        output$zoom_selector <- shiny::renderUI({
            suf    <- shiny::req(selected_suffix())
            method <- selected_method()
            pd    <- project_data()
            zooms <- find_gf_zooms(pd$name, suf, method = method)
            if (length(zooms) == 0) return(NULL)
            shiny::selectInput(ns("zoom"), "Zoom Region",
                choices  = c("Full view" = "none", setNames(zooms, zooms)),
                selected = "none"
            )
        })

        selected_variant <- shiny::reactive(input$piemap_variant %||% "base")
        selected_zoom    <- shiny::reactive(input$zoom %||% "none")
        selected_points  <- shiny::reactive(isTRUE(input$points))

        # geometric_offset writes a conditioning diagnostic (see section 10b of
        # scripts/geometric_offset.R). Surface its headline number here rather
        # than leaving it unseen inside a TSV — same reason mod_climate.R lifts
        # dbmem_diagnostics.tsv's `status` into the tab. When it fires, BOTH the
        # importance bars and the offset ranking are artefacts of the predictor
        # block, so the text goes on both cards.
        geometric_conditioning_note <- shiny::reactive({
            if (!identical(selected_method(), "geometric_offset")) return(NULL)
            suf <- selected_suffix()
            if (is.null(suf)) return(NULL)
            d <- load_geometric_offset_diagnostics(project_data()$name, suf)
            if (length(d) == 0) return(NULL)
            status <- d[["status"]] %||% ""
            if (!identical(status, "ok")) {
                return(paste0("Conditioning diagnostic: ", status,
                              " (geometric_offset_diagnostics.tsv)."))
            }
            share <- suppressWarnings(as.numeric(d[["offset_share_from_null_directions"]]))
            share_site <- suppressWarnings(as.numeric(d[["offset_share_from_null_directions_site"]]))
            if (!is.finite(share)) return(NULL)
            base <- paste0(
                round(100 * share, 1), "% of the realised offset (",
                round(100 * share_site, 1), "% at the sampled sites) comes from eigen-directions ",
                "carrying < ", 100 * as.numeric(d[["env_var_tolerance"]] %||% 1e-3),
                "% of the environmental variance each — ",
                d[["n_predictors"]], " predictors on ", d[["n_distinct_env_rows"]],
                " distinct environment rows, condition number ",
                d[["env_cov_condition_number"]], "."
            )
            if (share > 0.5 || (is.finite(share_site) && share_site > 0.5)) {
                paste0("ILL-CONDITIONED PREDICTOR BLOCK. ", base,
                       " The ranking and the importance bars are then set by the environment, ",
                       "not by the markers. Reduce/decorrelate Climate.predictors before ",
                       "interpreting either. Full breakdown: geometric_offset_diagnostics.tsv.")
            } else {
                paste0("Conditioning OK. ", base)
            }
        })

        # ── Importance images ──────────────────────────────────────────────────
        shiny::observe({
            suf    <- shiny::req(selected_suffix())
            method <- selected_method()
            pd     <- project_data()

            mod_image_card_server("overall_importance",
                path    = shiny::reactive(gf_importance_path(pd$name, suf, "overall", method = method)),
                title   = shiny::reactive("Overall Variable Importance"),
                dl_name = shiny::reactive(paste0("overall_importance_", suf)),
                note    = shiny::reactive(help_note("overall_importance",
                                                    extra = geometric_conditioning_note()))
            )
            mod_image_card_server("cumulative_importance",
                path    = shiny::reactive(gf_importance_path(pd$name, suf, "cumulative", method = method)),
                title   = shiny::reactive("Cumulative Importance"),
                dl_name = shiny::reactive(paste0("cumulative_importance_", suf)),
                note    = shiny::reactive(help_note("cumulative_importance"))
            )
        })

        # ── Single selector-driven offset piemap ───────────────────────────────
        # Points ("clear map") is trait/metric-independent and has no zoom
        # companion — when on, it always shows the project-level points file,
        # ignoring the zoom/variant selectors.
        piemap_path <- shiny::reactive({
            suf     <- shiny::req(selected_suffix())
            method  <- selected_method()
            pd      <- project_data()
            if (selected_points()) {
                return(gf_offset_piemap_path(pd$name, suf, "points", method = method))
            }
            variant <- selected_variant()
            zoom    <- selected_zoom()
            if (zoom != "none") {
                mod_path(pd$name, MOD_MALAD, "plots", method, suf, "zoom",
                         paste0(zoom, ".png"))
            } else {
                gf_offset_piemap_path(pd$name, suf, variant, method = method)
            }
        })

        piemap_title <- shiny::reactive({
            if (selected_points()) return("Genetic Offset (points)")
            variant <- selected_variant()
            zoom    <- selected_zoom()
            label   <- c(
                base         = "Genetic Offset",
                tajima_d     = "Offset \u00d7 Tajima's D",
                pi_diversity = "Offset \u00d7 Pi Diversity"
            )[variant]
            if (zoom != "none") paste0(label, " \u2014 Zoom: ", zoom) else label
        })

        mod_image_card_server("offset_piemap",
            path        = piemap_path,
            title       = piemap_title,
            dl_name     = shiny::reactive(paste0("offset_piemap_",
                                                  if (selected_points()) "points" else selected_variant(),
                                                  "_", selected_suffix() %||% "gf")),
            placeholder = shiny::reactive("Run mode=maladaptation to generate offset results"),
            note        = shiny::reactive(help_note("offset_piemap",
                extra = if (identical(selected_method(), "rda_offset"))
                    paste0("RDA offset: report as a RELATIVE SITE RANKING only, not an ",
                           "absolute magnitude — eigenvalue weighting + candidate-set size ",
                           "make the raw number non-comparable across runs/datasets. ",
                           "Unvalidated absent common-garden data (docs/rda_research.md B20–B21).")
                else geometric_conditioning_note()))
        )

        # ── Site table ─────────────────────────────────────────────────────────
        site_data <- shiny::reactive({
            suf    <- selected_suffix()
            method <- selected_method()
            if (is.null(suf)) return(data.table::data.table())
            p  <- gf_site_table_path(project_data()$name, suf, method = method)
            # mtime fingerprint: re-running for the same suffix otherwise serves stale data
            fp <- if (file.exists(p)) as.character(file.info(p)$mtime) else "missing"
            load_cached(paste0("mala_site_", project_data()$name, "_", method, "_", suf),
                        function() {
                p <- gf_site_table_path(project_data()$name, suf, method = method)
                if (!file_ok(p)) return(data.table::data.table())
                dt <- data.table::fread(p, colClasses = c("site" = "character",
                                                            "sample" = "character"))
                # Join pop stats (TajimaD, PI) if available
                tajima_p <- mod_path(project_data()$name, MOD_STRUCT,
                                     "tables", "pop_stats", "tajima_d_by_pop.tsv")
                pi_p     <- mod_path(project_data()$name, MOD_STRUCT,
                                     "tables", "pop_stats", "pi_diversity_by_pop.tsv")
                if (file_ok(tajima_p)) {
                    taj <- data.table::fread(tajima_p, colClasses = c("site" = "character"))
                    dt  <- taj[dt, on = "site"]
                }
                if (file_ok(pi_p)) {
                    pi_dt <- data.table::fread(pi_p, colClasses = c("site" = "character"))
                    dt    <- pi_dt[dt, on = "site"]
                }
                dt
            }, fingerprint = fp)
        })

        output$site_table <- DT::renderDataTable({
            dt <- site_data()
            safe_datatable(
                as.data.frame(dt),
                extensions = "Buttons",
                options    = list(
                    dom       = "Bfrtip",
                    buttons   = list("csv"),
                    scrollX   = TRUE,
                    pageLength = 20
                ),
                rownames = FALSE
            )
        })

        output$dl_site <- shiny::downloadHandler(
            filename = function() {
                paste0("genetic_offset_site_", project_data()$name,
                       "_", selected_method(), "_", selected_suffix(), ".csv")
            },
            content = function(file) {
                dt <- site_data()
                utils::write.csv(as.data.frame(dt), file, row.names = FALSE)
            }
        )

        # ── Cross-model comparison ──────────────────────────────────────────────
        # Track user's A/B picks so they survive renderUI re-renders.
        compare_sel_a <- shiny::reactiveVal("")
        compare_sel_b <- shiny::reactiveVal("")
        shiny::observeEvent(input$compare_model_a, compare_sel_a(input$compare_model_a),
                            ignoreNULL = TRUE, ignoreInit = TRUE)
        shiny::observeEvent(input$compare_model_b, compare_sel_b(input$compare_model_b),
                            ignoreNULL = TRUE, ignoreInit = TRUE)

        # Render model pickers from mala_models() — reacts to snp_sets_trigger().
        output$compare_model_a_ui <- shiny::renderUI({
            models <- mala_models()
            choices <- c("— select —" = "", models)
            cur <- compare_sel_a()
            sel <- if (nzchar(cur) && cur %in% models) cur else ""
            shiny::selectInput(ns("compare_model_a"), NULL,
                               choices = choices, selected = sel)
        })
        output$compare_model_b_ui <- shiny::renderUI({
            models <- mala_models()
            choices <- c("— select —" = "", models)
            cur <- compare_sel_b()
            sel <- if (nzchar(cur) && cur %in% models) cur else ""
            shiny::selectInput(ns("compare_model_b"), NULL,
                               choices = choices, selected = sel)
        })

        # Spatial-tag mismatch note
        output$compare_spatial_note <- shiny::renderUI({
            key_a <- input$compare_model_a
            key_b <- input$compare_model_b
            if (!nzchar(key_a) || !nzchar(key_b)) return(NULL)
            pa <- parse_model_key(key_a); pb <- parse_model_key(key_b)
            if (is.null(pa) || is.null(pb)) return(NULL)
            tag_a <- gsub("^.*_(spatial|nospatial)$", "\\1", pa$suffix)
            tag_b <- gsub("^.*_(spatial|nospatial)$", "\\1", pb$suffix)
            if (!identical(tag_a, tag_b)) {
                htmltools::div(
                    class = "alert alert-warning py-1 px-2 small mt-1 mb-0",
                    htmltools::HTML("&#9888; Spatial tags differ (",
                                   htmltools::strong(tag_a),
                                   " vs ", htmltools::strong(tag_b),
                                   "). Rank comparison is still valid; absolute magnitudes are not directly comparable.")
                )
            } else NULL
        })

        # Track whether cache exists for current selection
        compare_stats <- shiny::reactiveVal(NULL)

        output$compare_cache_badge <- shiny::renderUI({
            if (!is.null(compare_stats()))
                htmltools::span(class = "badge text-bg-success", "Cached")
            else NULL
        })

        # Run comparison on button click
        shiny::observeEvent(input$run_compare, {
            key_a <- input$compare_model_a
            key_b <- input$compare_model_b
            if (!nzchar(key_a) || !nzchar(key_b) || identical(key_a, key_b)) {
                shiny::showNotification("Select two different models.", type = "warning")
                return()
            }
            # Try loading from cache first
            stats <- load_compare_stats(project_data()$name, key_a, key_b)
            if (!is.null(stats)) {
                compare_stats(stats)
                shiny::showNotification("Loaded from cache.", type = "message", duration = 2)
                return()
            }

            # Find climate inputs from config
            cfg    <- project_data()$config
            preds  <- config_get(cfg, "Climate", "predictors") %||% character(0)
            env_sp <- mod_path(project_data()$name, MOD_CLIMATE, "tables", "present",
                               "climate_present_site.tsv")
            env_ap <- mod_path(project_data()$name, MOD_CLIMATE, "tables", "present",
                               "climate_present_all.tsv")

            # Use first future scenario in config
            ssps <- config_get(cfg, "Future", "ssp")
            yrs  <- config_get(cfg, "Future", "year")
            ssp  <- if (length(ssps) > 0) ssps[[1]] else "585"
            yr   <- if (length(yrs)  > 0) yrs[[1]]  else "2080"
            env_af <- mod_path(project_data()$name, MOD_CLIMATE, "tables", "future",
                               paste0("climate_future_year", yr, "_ssp", ssp, "_all.tsv"))
            scen_label <- paste0("ssp", ssp, "_", yr)

            # Canonical present-climate rasterstack = spatial template for ExDet/disagree rasters.
            # All predictor bands, aligned to climate_present_all.tsv cell IDs.
            # (offset_raster.tif is single-band and cannot serve as template.)
            rast_template <- mod_path(project_data()$name, MOD_CLIMATE,
                                      "rasters", "present", "climate_present_rasterstack.tif")
            if (!file_ok(rast_template)) {
                shiny::showNotification(
                    "Present climate raster not found. Run mode=maladaptation (climate step) first.",
                    type = "error"
                )
                return()
            }
            if (!file_ok(env_af)) {
                shiny::showNotification(
                    paste0("Future climate table not found for ", scen_label,
                           ". Check Future.ssp / Future.year config."),
                    type = "error"
                )
                return()
            }

            shiny::showNotification("Running comparison (may take ~30 s)...",
                                    id = "cmp_run", type = "message", duration = NULL)
            ok <- tryCatch(
                run_compare_offsets(
                    project      = project_data()$name,
                    key_a        = key_a,
                    key_b        = key_b,
                    env_site_pres = env_sp,
                    env_all_pres  = env_ap,
                    env_all_fut   = env_af,
                    predictors    = preds,
                    raster_tif    = rast_template,
                    top_k         = 10L,
                    scenario_label = scen_label
                ),
                error = function(e) { message("compare_offsets ERROR: ", e$message); FALSE }
            )
            shiny::removeNotification("cmp_run")
            if (ok) {
                compare_stats(load_compare_stats(project_data()$name, key_a, key_b))
                shiny::showNotification("Comparison complete.", type = "message")
            } else {
                shiny::showNotification("Comparison failed. Check logs.", type = "error")
            }
        })

        # ── Render comparison results tabs ────────────────────────────────────
        output$compare_results_ui <- shiny::renderUI({
            stats <- compare_stats()
            key_a <- input$compare_model_a
            key_b <- input$compare_model_b

            if (!nzchar(key_a %||% "") || !nzchar(key_b %||% "")) {
                return(compare_empty_state())
            }
            if (is.null(stats)) {
                return(htmltools::div(
                    class = "text-muted small fst-italic py-3",
                    "Click Compare to compute extrapolation diagnostics."
                ))
            }

            # KPI strip
            kpi_strip <- htmltools::div(
                class = "mb-3",
                compare_stat_box(
                    value  = stats$spearman_rho,
                    label  = "Spearman ρ (sites)",
                    p_val  = stats$spearman_p
                ),
                compare_stat_box(
                    value  = stats$kendall_tau,
                    label  = "Kendall τ (sites)",
                    p_val  = stats$kendall_p
                ),
                compare_stat_box(
                    value  = paste0(round(stats$jaccard_top_k * 100), "%"),
                    label  = paste0("Top-", stats$top_k, " Jaccard")
                ),
                if (!is.null(stats$kendall_w) && !is.na(stats$kendall_w))
                    compare_stat_box(
                        value = stats$kendall_w,
                        label = "Kendall W (N-way)"
                    )
            )

            bslib::navset_card_tab(
                id = ns("compare_tabs"),
                bslib::nav_panel(
                    "Novelty surface",
                    bslib::card_body(
                        kpi_strip,
                        shiny::uiOutput(ns("compare_novelty_ui"))
                    )
                ),
                bslib::nav_panel(
                    "Disagreement × Novelty",
                    bslib::card_body(
                        htmltools::p(
                            class = "small text-muted",
                            "Rank-transformed offset rasters; |Δrank| divergence map. High disagreement in non-analog climate cells (ExDet ≥ 0) indicates the two models extrapolate differently — as expected for linear (Geometric) vs plateau-bounded (Gradient Forest) offsets."
                        ),
                        shiny::uiOutput(ns("compare_disagree_ui"))
                    )
                ),
                bslib::nav_panel(
                    "Rank concordance",
                    bslib::card_body(
                        htmltools::div(
                            class = "alert alert-info py-1 px-2 small mb-2",
                            "Per-site concordance only — site agreement is a weaker question than extrapolation behaviour. See tabs 1-2 for the extrapolation-focused view."
                        ),
                        shiny::uiOutput(ns("compare_concordance_ui"))
                    )
                ),
                bslib::nav_panel(
                    "Site stability",
                    bslib::card_body(
                        shiny::uiOutput(ns("compare_stability_ui"))
                    )
                ),
                bslib::nav_panel(
                    "N-way concordance",
                    bslib::card_body(
                        shiny::uiOutput(ns("compare_nway_ui"))
                    )
                )
            )
        })

        # ── Tab sub-renderers ──────────────────────────────────────────────────
        output$compare_novelty_ui <- shiny::renderUI({
            stats <- compare_stats(); if (is.null(stats)) return(NULL)
            pd <- project_data()
            key_a <- input$compare_model_a; key_b <- input$compare_model_b
            scen  <- "future"  # default; could be parameterised
            nov_p <- novelty_raster_path(pd$name, scen)
            if (!file_ok(nov_p)) {
                return(htmltools::div(class = "text-muted small", "Novelty raster not computed yet."))
            }
            # Show cached ExDet table if present
            nov_t <- file.path(dirname(nov_p), "exdet_novelty.tsv")
            nov_info <- if (file_ok(nov_t)) {
                tryCatch({
                    dt <- data.table::fread(nov_t)
                    htmltools::div(
                        compare_stat_box(paste0(dt$pct_cells_novel, "%"), "Novel cells (ExDet ≥ 0)"),
                        compare_stat_box(dt$max_nt2, "Max NT2 (Mahalanobis)")
                    )
                }, error = function(e) NULL)
            } else NULL
            htmltools::tagList(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    htmltools::span(class = "fw-semibold", "Novelty surface"),
                    help_note("compare_novelty")
                ),
                nov_info,
                htmltools::p(
                    class = "small text-muted",
                    "ExDet NT1: univariate novelty (positive = within training range). ",
                    "NT2: multivariate combinatorial novelty (Mahalanobis to training centroid). ",
                    "Reference = present climate at sampled sites."
                )
            )
        })

        output$compare_disagree_ui <- shiny::renderUI({
            stats <- compare_stats(); if (is.null(stats)) return(NULL)
            cor_txt <- if (!is.null(stats$disagree_nt2_spear) && !is.na(stats$disagree_nt2_spear))
                paste0("Spearman(Δrank, NT2) = ", stats$disagree_nt2_spear,
                       " — higher = disagreement tracks climate novelty")
            else "Disagreement-novelty Spearman not available."
            htmltools::div(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    htmltools::span(class = "fw-semibold", "Disagreement × Novelty"),
                    help_note("compare_disagree")
                ),
                htmltools::p(class = "small", cor_txt),
                htmltools::p(class = "small text-muted",
                    "Disagree raster: ", file.path(
                        model_compare_dir(project_data()$name,
                                          input$compare_model_a, input$compare_model_b),
                        "disagree_rank.tif"))
            )
        })

        output$compare_concordance_ui <- shiny::renderUI({
            stats <- compare_stats(); if (is.null(stats)) return(NULL)
            dt <- load_compare_site_offsets(project_data()$name,
                                            input$compare_model_a, input$compare_model_b)
            if (nrow(dt) == 0) return(htmltools::div(class = "text-muted small", "No data."))
            htmltools::tagList(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    htmltools::span(class = "fw-semibold", "Rank concordance"),
                    help_note("compare_concordance")
                ),
                htmltools::div(
                    compare_stat_box(stats$spearman_rho, "Spearman ρ", stats$spearman_p,
                                     accent = "#1B7A6E"),
                    compare_stat_box(stats$kendall_tau,  "Kendall τ",  stats$kendall_p),
                    if (!is.null(stats$dutilleul_p) && !is.na(stats$dutilleul_p))
                        compare_stat_box(
                            signif(stats$dutilleul_p, 3),
                            paste0("Dutilleul p (nₑₓₓ≈",
                                   round(stats$dutilleul_n_eff), ")"),
                            accent = "#5b9bd5"
                        )
                ),
                # Scatter: Model A vs Model B offset per site
                shiny::renderPlot({
                    ggplot2::ggplot(as.data.frame(dt),
                                    ggplot2::aes(x = offset_a, y = offset_b,
                                                  label = site)) +
                        ggplot2::geom_point(color = "#1B7A6E", alpha = 0.7, size = 2.5) +
                        ggplot2::geom_smooth(method = "lm", color = "#888", se = TRUE, linewidth = 0.8) +
                        ggplot2::labs(
                            x = stats$model_a_label,
                            y = stats$model_b_label
                        ) +
                        ggplot2::theme_minimal(base_size = 12)
                }, height = 320)
            )
        })

        output$compare_stability_ui <- shiny::renderUI({
            stats <- compare_stats(); if (is.null(stats)) return(NULL)
            dt    <- load_compare_rank_stability(project_data()$name,
                                                  input$compare_model_a, input$compare_model_b)
            htmltools::tagList(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    htmltools::span(class = "fw-semibold", "Site stability"),
                    help_note("compare_stability")
                ),
                htmltools::div(
                    compare_stat_box(
                        paste0(round(stats$jaccard_top_k * 100), "%"),
                        paste0("Jaccard top-", stats$top_k, " high-risk sites"),
                        accent = "#d97706"
                    )
                ),
                if (nrow(dt) > 0) {
                    shiny::renderPlot({
                        dd <- as.data.frame(dt)
                        dd$site <- factor(dd$site, levels = dd$site[order(dd$rank_diff, decreasing = TRUE)])
                        ggplot2::ggplot(utils::head(dd[order(-dd$rank_diff), ], 20),
                                        ggplot2::aes(x = site, y = rank_diff)) +
                            ggplot2::geom_col(fill = "#d97706", alpha = 0.8) +
                            ggplot2::coord_flip() +
                            ggplot2::labs(x = NULL, y = "|Rank A − Rank B|") +
                            ggplot2::theme_minimal(base_size = 11)
                    }, height = 320)
                }
            )
        })

        output$compare_nway_ui <- shiny::renderUI({
            stats <- compare_stats(); if (is.null(stats)) return(NULL)
            if (is.null(stats$kendall_w) || is.na(stats$kendall_w)) {
                return(htmltools::div(
                    class = "text-muted small",
                    "N-way concordance requires vegan::kendall.global (vegan package) installed in the Docker container."
                ))
            }
            htmltools::tagList(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    htmltools::span(class = "fw-semibold", "N-way concordance"),
                    help_note("compare_nway")
                ),
                compare_stat_box(stats$kendall_w, "Kendall's W", stats$kendall_w_p,
                                  accent = "#1B7A6E"),
                htmltools::p(
                    class = "small text-muted",
                    "Kendall's W = 0 (no concordance) to 1 (perfect concordance) across ",
                    stats$n_models_nway, " models."
                )
            )
        })
    })
}

#' maladaptation tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_maladaptation_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "malad",

        # Model selector + the badge cluster the server renders. offset_summary
        # joins them here (it used to sit above the site table); it is the live
        # result for the selected model, so it reads as the bar's headline.
        htmltools::div(
            class = "control-bar lab-malad-modelbar d-flex align-items-center gap-2 flex-wrap",
            shiny::uiOutput(ns("model_selector")),
            shiny::uiOutput(ns("snp_count_badge")),
            shiny::uiOutput(ns("method_params_badge")),
            shiny::uiOutput(ns("climate_scenario_badge")),
            shiny::uiOutput(ns("sensitivity_badge")),
            htmltools::div(class = "ms-auto d-flex align-items-center gap-2",
                           shiny::uiOutput(ns("set_details_btn")),
                           shiny::uiOutput(ns("set_delete_btn")))
        ),

        bslib::layout_columns(
            col_widths = c(4, 4, 4),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── COLUMN 1: the offset map ──────────────────────────────────────
            htmltools::div(
                class = "lab-malad-col lab-malad-col-map",
                # One inline strip, not a control-bar box. The old version wrapped
                # the two selectors in a layout_columns and put the Points switch
                # on its own row beneath: on a project with neither a piemap
                # variant nor a zoom map (both renderUI's return NULL when there
                # is only one choice) that was a 64px bordered box containing two
                # zero-height grid columns and a 23px switch. Flex + gap means the
                # row is exactly as tall as whatever actually rendered.
                htmltools::div(
                    class = "lab-mapbar lab-malad-mapbar",
                    shiny::uiOutput(ns("piemap_variant_ui")),
                    shiny::uiOutput(ns("zoom_selector")),
                    bslib::input_switch(ns("points"), "Points", value = FALSE),
                    # The offset range/mean/n-sites badge sits with the map
                    # rather than in the model bar: it describes what this map
                    # is showing, and the map strip is where the eye already is
                    # when reading it.
                    htmltools::div(class = "lab-malad-offsetnote ms-auto",
                                   shiny::uiOutput(ns("offset_summary")))
                ),
                lab_hero(htmltools::div(class = "piemap-container",
                                        mod_image_card_ui(ns("offset_piemap"))))
            ),

            # ── COLUMN 2: predictor importance, stacked ───────────────────────
            htmltools::div(
                class = "lab-malad-col lab-malad-col-imp",
                mod_image_card_ui(ns("overall_importance")),
                mod_image_card_ui(ns("cumulative_importance"))
            ),

            # ── COLUMN 3: per-site offsets ────────────────────────────────────
            htmltools::div(
                class = "lab-malad-col lab-malad-col-tbl",
                bslib::card(
                    class = "lab-dock-card lab-dock-offset",
                    bslib::card_header(
                        class = "d-flex justify-content-between align-items-center",
                        htmltools::span(bsicons::bs_icon("geo-alt"), " Site offsets"),
                        shiny::downloadButton(ns("dl_site"), "CSV",
                                              class = "btn-sm btn-outline-secondary")
                    ),
                    bslib::card_body(DT::DTOutput(ns("site_table")))
                )
            )
        ),

        # Cross-model comparison: opt-in and expensive -> collapsed.
        bslib::accordion(
            id = ns("malad_sections"),
            class = "lab-below-fold",
            open = FALSE, multiple = TRUE,
            bslib::accordion_panel(
                value = "compare",
                title = htmltools::tagList(bsicons::bs_icon("shuffle"),
                                           " Cross-model comparison"),
                compare_honesty_banner(),
                bslib::layout_columns(
                    col_widths = c(6, 6),
                    htmltools::div(class = "control-bar lab-cmp-a",
                                   shiny::uiOutput(ns("compare_model_a_ui"))),
                    htmltools::div(class = "control-bar lab-cmp-b",
                                   shiny::uiOutput(ns("compare_model_b_ui")))
                ),
                shiny::uiOutput(ns("compare_spatial_note")),
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    shiny::actionButton(ns("run_compare"), "Compare",
                                        class = "btn-sm btn-primary"),
                    shiny::uiOutput(ns("compare_cache_badge"))
                ),
                shiny::uiOutput(ns("compare_results_ui"))
            )
        )
    )
}

