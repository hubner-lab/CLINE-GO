# mod_dashboard_outputs.R — outputs that exist only for the dashboard layout.
#
# Promoted from scripts/layout_lab/ on 2026-08-29. Each function ADDS outputs the
# original module server never defined (the KPI strip everywhere, plus
# PreStructure's Q-matrix table and GEA's suspend override). They run as a second
# moduleServer on the SAME id, so they share the namespace and can read the real
# server's inputs -- which is how lab_90_dispatch.R wired them, and why they stay
# separate functions rather than being spliced into the module servers.
#
# app_server.R calls each one directly after its module server.


mod_prestructure_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {

        selected_k <- shiny::reactive({
            k <- input$k
            if (!is.null(k)) return(as.integer(k))
            pd <- project_data()
            ks <- find_k_values(pd$name)
            if (length(ks) == 0) return(NULL)
            if (!is.na(pd$k_best) && pd$k_best %in% ks) pd$k_best else ks[1]
        })

        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))

            sv <- function(step, metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }

            # NOT find_k_range(): it parses whichever cross_entropy_K*.png sorts
            # first, and a stale K2-6 file sitting beside the current K2-7.0 one
            # makes it report the WRONG range (logged in
            # docs/pipeline_improvement_requests.md). The K* directories on disk
            # are the K values actually computed, so derive from those and fall
            # back to the summary table, which is authoritative.
            ks <- find_k_values(pd$name)
            k_range_txt <- if (length(ks) > 0) {
                paste0(min(ks), " – ", max(ks))
            } else sv("prestructure", "K_range")

            k_best_txt <- if (!is.null(pd$k_best) && !is.na(pd$k_best))
                as.character(pd$k_best) else "—"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("K best",  k_best_txt,  "primary", "diagram-3"),
                lab_value_box("K tested", k_range_txt, "info",   "list-ol"),
                lab_value_box("Samples",
                              sv("processing", "samples_after_filtering"),
                              "success", "people-fill"),
                lab_value_box("SNPs (LD-pruned)",
                              sv("processing", "snps_after_ld_pruning"),
                              "info", "scissors")
            )
        })

        output$clusters_table <- DT::renderDataTable({
            pd <- project_data()
            k  <- selected_k()
            shiny::req(k)
            df <- as.data.frame(load_clusters(pd$name, k))
            shiny::validate(shiny::need(nrow(df) > 0, "No Q-matrix for this K"))
            # C1..CK are ancestry proportions -- 3 dp is plenty and keeps the
            # table narrow enough for a 5/12 column.
            num <- vapply(df, is.numeric, logical(1))
            df[num] <- lapply(df[num], function(x) round(x, 3))
            safe_datatable(
                df,
                options  = list(dom = "tp", scrollX = TRUE, pageLength = 6),
                rownames = FALSE
            )
        })
    })
}


mod_structure_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))

            sv <- function(step, metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))
            bp <- function(x) {
                v <- num(x)
                if (is.na(v)) "—"
                else if (v >= 1e6) paste0(round(v / 1e6, 2), " Mb")
                else if (v >= 1e3) paste0(round(v / 1e3, 1), " kb")
                else paste0(v, " bp")
            }

            k_best <- sv("structure", "K_best")
            if (is.na(k_best)) k_best <- if (!is.na(pd$k_best)) pd$k_best else "—"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("K best", as.character(k_best), "primary", "diagram-3"),
                lab_value_box("Climate predictors",
                              as.character(sv("structure", "n_climate_variables", "—")),
                              "info", "thermometer-half"),
                lab_value_box("LD half-decay",
                              bp(sv("structure", "ld_decay_half_decay_genome_wide")),
                              "success", "bar-chart-line"),
                lab_value_box("LD r² = 0.2",
                              bp(sv("structure", "ld_decay_r2_02_genome_wide")),
                              "info", "rulers")
            )
        })
    })
}


mod_climate_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {

        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "climate" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            unexp     <- num(sv("varpart_unexplained_R2adj"))
            unexp_txt <- if (is.na(unexp)) "—" else paste0(round(unexp, 1), "%")
            # Nearly all variance unexplained is the interesting/bad case.
            unexp_theme <- if (is.na(unexp)) "secondary"
                           else if (unexp >= 90) "danger"
                           else if (unexp >= 75) "warning" else "success"

            # Climate UNIQUE variance is the headline number of this module:
            # what climate explains once structure and geography are removed.
            # No yes/no boxes here -- the KPI strip carries numbers, and the
            # confounding call is already made by confounding_badge below.
            clim      <- num(sv("varpart_climate_R2adj"))
            clim_txt  <- if (is.na(clim)) "—" else paste0(round(clim, 2), "%")
            clim_theme <- if (is.na(clim)) "secondary"
                          else if (clim >= 5) "success"
                          else if (clim >= 1) "warning" else "danger"

            n_inv       <- num(sv("n_invariant_predictors"))
            n_inv_theme <- if (!is.na(n_inv) && n_inv > 0) "warning" else "success"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Climate unique (R²adj)", clim_txt, clim_theme, "thermometer-half"),
                lab_value_box("Unexplained", unexp_txt, unexp_theme, "pie-chart"),
                lab_value_box("Invariant predictors",
                              as.character(sv("n_invariant_predictors", "—")),
                              n_inv_theme, "exclamation-triangle"),
                lab_value_box("dbMEM vectors",
                              as.character(sv("dbMEM_n_mem_positive", "—")),
                              "info", "diagram-2")
            )
        })
    })
}


mod_traits_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "traits" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            # Worst missingness across traits: not in pipeline_summary.tsv, so
            # derive it from trait_summary.tsv, which the app already loads.
            miss_txt <- "—"; miss_theme <- "secondary"
            p <- trait_table_path(pd$name, "trait_summary")
            if (file_ok(p)) {
                ts <- tryCatch(as.data.frame(data.table::fread(p)),
                               error = function(e) NULL)
                if (!is.null(ts) && "pct_missing" %in% names(ts) && nrow(ts) > 0) {
                    mx <- max(num(ts$pct_missing), na.rm = TRUE)
                    if (is.finite(mx)) {
                        miss_txt   <- paste0(round(mx, 1), "%")
                        miss_theme <- if (mx >= 20) "danger"
                                      else if (mx >= 5) "warning" else "success"
                    }
                }
            }

            n_inv       <- num(sv("n_invariant_traits"))
            n_inv_theme <- if (!is.na(n_inv) && n_inv > 0) "warning" else "success"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Traits", as.character(sv("n_traits", "—")),
                              "primary", "rulers"),
                lab_value_box("With missing data",
                              as.character(sv("n_traits_with_missing", "—")),
                              "info", "question-octagon"),
                lab_value_box("Worst missingness", miss_txt, miss_theme, "percent"),
                lab_value_box("Invariant traits",
                              as.character(sv("n_invariant_traits", "—")),
                              n_inv_theme, "exclamation-triangle")
            )
        })
    })
}


mod_gea_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        for (out_id in c("threshold_bar", "filter_bar", "wza_collapse_note")) {
            ok <- tryCatch({
                shiny::outputOptions(output, out_id, suspendWhenHidden = FALSE)
                TRUE
            }, error = function(e) {
                message("[lab] could not un-suspend '", out_id, "': ", conditionMessage(e))
                FALSE
            })
            if (!ok) next
        }
    })
}


mod_gwas_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "gwas" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            meth <- if (nrow(s) == 0) character(0) else
                s$metric[s$step == "gwas" & grepl("^sig_snps_", s$metric)]
            n_snps     <- num(sv("selected_snps_total"))
            snps_theme <- if (is.na(n_snps)) "secondary"
                          else if (n_snps == 0) "danger" else "primary"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Significant SNPs",
                              as.character(sv("selected_snps_total", "—")),
                              snps_theme, "bullseye"),
                lab_value_box("Regions", as.character(sv("regions_combined", "—")),
                              "info", "geo-alt"),
                lab_value_box("Genes", as.character(sv("genes_found", "—")),
                              "success", "diagram-2"),
                lab_value_box("Methods",
                              if (length(meth) > 0) as.character(length(meth)) else "—",
                              "info", "layers")
            )
        })
    })
}


mod_geaxgwas_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(step, metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }
            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("GEA SNPs",  sv("gea",  "selected_snps_total"),
                              "primary", "arrow-up"),
                lab_value_box("GEA regions", sv("gea", "regions_combined"),
                              "info", "geo-alt"),
                lab_value_box("GWAS SNPs", sv("gwas", "selected_snps_total"),
                              "success", "arrow-down"),
                lab_value_box("GWAS regions", sv("gwas", "regions_combined"),
                              "info", "geo-alt")
            )
        })
    })
}


mod_pregea_dashboard_outputs <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "pregea" & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }
            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("LFMM K", sv("LFMM_K_recommended"),
                              "primary", "diagram-3"),
                lab_value_box("EMMAX #PCs", sv("EMMAX_n_pcs_recommended"),
                              "info", "sliders"),
                lab_value_box("RDA Condition() PCs", sv("RDA_condition_pcs_recommended"),
                              "info", "funnel"),
                lab_value_box("RDA axes", sv("RDA_axes_recommended"),
                              "success", "compass")
            )
        })
    })
}

