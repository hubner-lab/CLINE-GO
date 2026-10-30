

#' Is the GWAS missing-data alert a warning (154411)?
#'
#' Dropping the samples that miss a trait (DROP, complete case) is the standard
#' handling of a missing outcome, while filling it with the trait's MEAN or MEDIAN
#' shrinks the effect estimates. So a trait filled where values were missing, or one
#' with under half its samples, is the warning; DROP alone is information.
#' prepare_phenotypes.R counts n_available before any fill.
#' @param dt phenotype_missing_summary.tsv (trait, n_total, n_available, strategy)
#' @noRd
pheno_missing_alert_warns <- function(dt) {
    strat_col <- if ("missing_strategy" %in% names(dt)) "missing_strategy" else "strategy"
    low    <- any(dt$n_available / dt$n_total < 0.5, na.rm = TRUE)
    filled <- any(toupper(dt[[strat_col]]) %in% c("MEAN", "MEDIAN") &
                  dt$n_available < dt$n_total, na.rm = TRUE)
    low || filled
}

#' Phenotype Association tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param config_state reactiveValues from app_server.R ($working/$saved/$project),
#'   or NULL to omit the "Apply rules to config" button.
#' @noRd
mod_gwas_server <- function(id, project_data, run_trigger = NULL, config_state = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_GWAS

        # ── Data loading ───────────────────────────────────────────────────────
        methods <- shiny::reactive(find_assoc_methods(project_data()$name, module))

        # method -> list(adjust=,threshold=,family=) — pins non-univariate
        # methods to their registry rule in the matrix unless a cell override
        # exists (see R/fct_threshold_rules.R). GWAS excludes RDA entirely
        # (gwas_method_registry()), so this is a no-op fallback in practice —
        # kept for consistency with mod_gea.R and in case that ever changes.
        registry_defaults <- shiny::reactive({
            gea_method_significance_defaults(gwas_method_registry())
        })

        # config_rules: METHOD -> resolve_adjust() string ("bonf_0.05") — the
        # rule the PIPELINE FILES on disk were actually built with.
        config_rules <- shiny::reactive({
            pd <- project_data(); ms <- methods()
            if (length(ms) == 0) return(list())
            stats::setNames(
                lapply(ms, function(m) resolve_adjust(pd$config, m, module) %||% ""), ms)
        })

        all_method_pvalues <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_pvalues_cached(pd$name, module, pd$k_best)
        })

        all_method_wza_pvalues <- shiny::reactive({
            if (!regime_wza()) return(list())
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_wza_pvalues_cached(pd$name, module, pd$k_best)
        })

        effective_method_pvalues <- shiny::reactive({
            if (regime_wza()) all_method_wza_pvalues() else all_method_pvalues()
        })

        all_trait_names_rv <- shiny::reactive({
            pv <- all_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        traits <- shiny::reactive({
            sig <- effective_method_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        # ── Filter bar ─────────────────────────────────────────────────────────

        trait_colors <- shiny::reactive({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(character(0))
            trait_color_map(tr)
        })

        method_shapes <- shiny::reactive({
            ms <- methods()
            if (length(ms) == 0) return(character(0))
            method_shape_map(ms)
        })

        combo_counts <- shiny::reactive({
            all_snps <- effective_method_sigsnps()
            counts <- list()
            for (m in names(all_snps)) {
                dt <- all_snps[[m]]
                if (nrow(dt) > 0) {
                    tc <- dt[, .(n = .N), by = "trait"]
                    for (i in seq_len(nrow(tc)))
                        counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
                }
            }
            counts
        })

        # File fingerprint reactive — changes when pipeline regenerates pvalue TSVs.
        pvalue_fingerprint <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate on pipeline completion
            pd <- project_data()
            pvalues_file_fingerprint(pd$name, module, pd$k_best,
                                     regime = if (regime_wza()) "wza" else "snp")
        })

        # Per-cell significance threshold: "trait::method" -> cutoff p-value (all cells, NA if unavailable)
        combo_thresholds <- shiny::reactive({
            compute_method_thresholds(
                pvalues_list      = effective_method_pvalues(),
                type              = threshold_type(),
                value             = threshold_value(),
                overrides         = threshold_overrides(),
                registry_defaults = registry_defaults()
            )
        }) |> shiny::bindCache(threshold_type(), threshold_value(), regime_wza(),
                               pvalue_fingerprint(),
                               threshold_overrides_key(threshold_overrides(), registry_defaults()))

        # ── Filter-bar state: threshold rule, strategy, regime ─────────────────
        # Saved values else the config (strategy Union), saved only when moved
        # (filter_bar_state()); the GEAxGWAS "Fill from GWAS tab" button reads them.
        bar <- filter_bar_state(input, project_data, module, source = module,
                                registry_defaults = registry_defaults)

        active_strategy <- shiny::reactive(input$combine_strategy %||% bar$strategy())

        # ── SNP clumping distance (single param for both merging and overlap) ───
        # Starts at the saved value, else at the project's own distance (its LD-decay
        # distance unless the config sets bp); saved only when the box is moved. The
        # filter bar renders from this reactive, not from the input (see
        # clumping_distance_state()).
        snp_clumping_distance <- clumping_distance_state(
            input, project_data, "snp_clumping_distance", MOD_GWAS, "snp_clumping_distance", MOD_GWAS)

        # The region merge: the config's region_clumping_distance when it sets one,
        # else the clumping distance above, as the pipeline merges (acff31).
        region_merge <- shiny::reactive({
            pd <- project_data()
            config_region_clumping_distance(pd$config, MOD_GWAS, ld_decay_table_path(pd$name))
        })
        region_clumping_distance <- shiny::reactive(region_merge()$bp %||% snp_clumping_distance())

        # ── Regime (per-SNP vs WZA) ────────────────────────────────────────────
        regime_wza <- shiny::reactive(isTRUE(input$regime %||% bar$regime()))

        # ── WZA collapse note (appears when regime switch is ON) ───────────────
        output$wza_collapse_note <- shiny::renderUI({
            if (!regime_wza()) return(NULL)
            wza_collapse_note(wza_collapse_stats(all_method_wza_pvalues()))
        })

        # ── Threshold type + value ────────────────────────────────────────────
        threshold_type_raw <- shiny::reactive(input$threshold_type %||% bar$threshold_type())
        # Debounce type so switching to "qval" doesn't fire the full 2M-row scan instantly
        threshold_type  <- shiny::debounce(threshold_type_raw, 500)
        threshold_value_raw <- shiny::reactive({
            v <- input$threshold_value
            if (is.null(v) || is.na(v) || v <= 0) bar$threshold_value() else as.numeric(v)
        })
        threshold_value <- shiny::debounce(threshold_value_raw, 500)

        # ── Per-cell threshold overrides, matrix selection, rule popup ─────────
        # Shared with mod_gea.R / mod_gea_x_gwas.R — see
        # R/fct_combine.R::setup_matrix_rules_server().
        matrix_rules <- setup_matrix_rules_server(
            input, output, session, ns, input_prefix = "",
            project_data = project_data, module = module,
            methods = methods, all_traits = all_trait_names_rv,
            threshold_type = threshold_type, threshold_value = threshold_value,
            registry_defaults = registry_defaults, config_rules = config_rules,
            config_state = config_state, config_module_key = "GWAS",
            combo_counts = combo_counts
        )
        threshold_overrides <- matrix_rules$overrides

        output$threshold_hint <- shiny::renderUI({
            t <- threshold_type()
            hint <- switch(t,
                bonf = "α for Bonferroni correction", qval = "FDR q-value target (0–1)",
                top  = "Number of top SNPs per trait",
                toppct = "Percent of the SNPs per trait (rounded up to a count)",
                custom = "Raw p-value cutoff (e.g. 1e-5)", "")
            htmltools::span(class = "text-muted small mt-1", hint)
        })

        # Coerce threshold value when type changes
        shiny::observeEvent(input$threshold_type, {
            v <- input$threshold_value
            if (!threshold_value_valid_for_type(input$threshold_type, v)) {
                shiny::updateNumericInput(
                    session, "threshold_value",
                    value = threshold_value_default_for_type(input$threshold_type)
                )
            }
        }, ignoreInit = TRUE)

        # ── Always-present threshold bar ──────────────────────────────────────
        output$threshold_bar <- shiny::renderUI({
            project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    regime_value          = bar$regime(),
                    threshold_type_value  = bar$threshold_type(),
                    threshold_value_value = bar$threshold_value(),
                    regime_context        = "gwas",
                    show_apply_to_config  = !is.null(config_state)
                )
            })
        })

        # ── Interactive sig SNPs: computed from full pvalue tables + threshold ─
        effective_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            compute_method_sigsnps_cached(
                pvalues_list = effective_method_pvalues(),
                type         = threshold_type(),
                value        = threshold_value(),
                k            = pd$k_best,
                regime       = if (regime_wza()) "wza" else "snp",
                project      = pd$name,
                module       = module,
                cutoffs           = combo_thresholds(),
                overrides         = threshold_overrides(),
                registry_defaults = registry_defaults()
            )
        })

        # Derived from bindCached combo_thresholds — no redundant full-vector scan.
        combined_threshold_y <- shiny::reactive({
            ct <- combo_thresholds()
            if (length(ct) == 0) return(NULL)
            vals <- unlist(ct)
            vals <- vals[!is.na(vals) & vals > 0]
            if (length(vals) == 0) return(NULL)
            -log10(min(vals))  # min raw p = most stringent = highest y line
        })

        per_method_threshold_y <- shiny::reactive({
            m <- active_method(); t <- per_method_trait()
            if (is.null(m) || is.null(t)) return(NULL)
            ct  <- combo_thresholds()
            key <- paste0(t, "::", m)
            thr <- ct[[key]]
            if (is.null(thr) || is.na(thr) || thr <= 0) return(NULL)
            -log10(thr)
        })

        # ── I4: Config parameter badges ────────────────────────────────────────
        output$config_badges <- shiny::renderUI({
            pd <- project_data()
            k  <- pd$k_best
            if (is.na(k) && length(methods()) == 0) return(NULL)
            cfg     <- pd$config
            # Show the strategy this panel ACTUALLY computed with (the interactive
            # filter-bar selector), not GWAS.combine_method from the YAML — that key
            # governs the static pipeline tables only and the two can disagree.
            cmb     <- .normalize_strategy(active_strategy())
            # DROP when unset: the pipeline's default (common.smk), so the badge names
            # the strategy the tables were made with (1cbc36)
            missing_strat <- config_get(cfg, "GWAS", "missing_strategy", default = "DROP")
            config_badges_bar(
                if (!is.na(k)) config_badge("K", k, "bg-primary"),
                config_badge("combine", cmb),
                config_badge("missing", missing_strat),
                if (regime_wza()) config_badge("regime", "WZA", "bg-info"),
                if (!is.null(region_merge()$bp))
                    config_badge("region merge", .format_bp(region_merge()$bp)),
                clumping_issue_note(pd, MOD_GWAS),
                region_clumping_issue_note(region_merge())
            )
        })

        # ── E1/E2: Phenotype missing data alert ────────────────────────────────
        output$pheno_missing_alert <- shiny::renderUI({
            pd <- project_data()
            dt <- load_pheno_missing_summary(pd$name)
            if (nrow(dt) == 0) return(NULL)
            req_cols <- c("trait", "n_total", "n_available", "strategy")
            if (!all(req_cols %in% names(dt))) return(NULL)

            warn        <- pheno_missing_alert_warns(dt)
            alert_class <- if (warn) "alert-warning" else "alert-info"
            icon_name   <- if (warn) "exclamation-triangle-fill" else "info-circle-fill"
            strat_col   <- if ("missing_strategy" %in% names(dt)) "missing_strategy" else "strategy"

            rows <- lapply(seq_len(nrow(dt)), function(i) {
                row      <- dt[i, ]
                pct      <- round(100 * row$n_available / row$n_total)
                strat    <- toupper(as.character(row[[strat_col]]))
                row_class <- if (pct < 50) "text-warning fw-semibold" else ""
                htmltools::tags$li(
                    class = row_class,
                    htmltools::tags$code(row$trait),
                    paste0(" \u2014 ", row$n_available, "/", row$n_total,
                           " samples (", pct, "%), strategy: ", strat)
                )
            })

            htmltools::div(
                class = paste("alert d-flex gap-2 align-items-start mb-2", alert_class),
                bsicons::bs_icon(icon_name, class = "flex-shrink-0 mt-1"),
                htmltools::div(
                    htmltools::tags$strong("Phenotype data availability"),
                    " \u2014 samples with non-missing values per trait:",
                    htmltools::tags$ul(class = "mb-0 mt-1", rows)
                )
            )
        })

        # ── Filter bar UI (matrix + strategy + clumping) ──────────────────────
        output$filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = all_trait_names_rv(),   # full list — always show grid
                methods                     = methods(),
                trait_colors                = trait_colors(),
                combo_counts                = combo_counts(),
                combo_thresholds            = combo_thresholds(),
                # isolate()d: seeded from the user's move, never re-rendered by it
                default_strategy_value      = shiny::isolate(bar$strategy()),
                snp_clumping_distance_value = snp_clumping_distance(),
                region_clumping_distance_value = region_merge()$bp,
                # isolate()d — see setup_matrix_rules_server(): the matrix's
                # on/off selection must survive a re-render triggered by
                # combo_counts/combo_thresholds changing (master threshold edit).
                selected_pairs              = shiny::isolate(matrix_rules$selected_pairs()),
                overrides                   = threshold_overrides(),
                registry_defaults           = registry_defaults(),
                master_type                 = threshold_type(),
                master_value                = threshold_value()
            )
        })

        # ── Interactive sig SNPs ───────────────────────────────────────────────
        interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = effective_method_sigsnps(),
                tm_selection_json   = input$tm_selection,
                combo_counts        = combo_counts(),
                known_traits        = all_trait_names_rv(),
                strategy            = active_strategy(),
                clumping_distance   = snp_clumping_distance(),
                project_name        = project_data()$name,
                module              = module
            )
        })

        # ── Phenotype trait selector (for phenomap) — uses all trait names ──────
        output$trait_selector <- shiny::renderUI({
            tr <- all_trait_names_rv()
            if (length(tr) == 0)
                return(shiny::p("No traits found.", class = "text-muted small"))
            shiny::selectInput(ns("pheno_trait"), "Trait (phenomap)",
                               choices = tr, selected = tr[1])
        })

        selected_pheno_trait <- shiny::reactive(input$pheno_trait %||% all_trait_names_rv()[1])
        selected_points      <- shiny::reactive(isTRUE(input$points))

        # ── Phenomap ───────────────────────────────────────────────────────────
        # Points ("clear map") is trait-independent — ignores the trait selector.
        output$phenomap_content <- shiny::renderUI({
            tr <- selected_pheno_trait()
            if (is.null(tr)) return(plot_placeholder("Select a trait"))
            pd   <- project_data()
            path <- pheno_piemap_path(pd$name, tr, points = selected_points())
            if (file_ok(path)) {
                shiny::imageOutput(ns("phenomap_img"), height = "auto", width = "100%")
            } else {
                plot_placeholder("Phenotype map not available",
                    "Run mode=gwas to generate phenotype piemaps")
            }
        })

        output$phenomap_img <- shiny::renderImage({
            tr <- shiny::req(selected_pheno_trait())
            pd <- project_data()
            p  <- pheno_piemap_path(pd$name, tr, points = selected_points())
            shiny::validate(shiny::need(file_ok(p), "Phenomap not found"))
            list(src = p, contentType = "image/png", width = "100%", alt = paste("Phenomap", tr))
        }, deleteFile = FALSE)

        # ── WZA path overrides for Manhattan ──────────────────────────────────
        combined_wza_bg <- shiny::reactive({
            if (!regime_wza()) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            combined_manhattan_wza_bg_path(pd$name, module, k)
        })
        combined_wza_coords <- shiny::reactive({
            if (!regime_wza()) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            combined_manhattan_wza_coords_path(pd$name, module, k)
        })

        # ── Interactive region explorer ────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = interactive_sigsnps,
            region_distance     = region_clumping_distance,
            regime              = regime_wza
        )

        # SNPs re-stamped with region_id from the live-computed regions (not the
        # static pipeline file) so a Manhattan click always selects the region
        # whose bounds match what the table/rectangles show. copy() is required —
        # assign_region_ids_from_regions() mutates by reference and
        # interactive_sigsnps() is a shared cached reactive.
        plotted_sigsnps <- shiny::reactive({
            snps <- interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(snps)
            assign_region_ids_from_regions(data.table::copy(snps), explorer$computed_regions())
        })

        # ── Combined Manhattan ─────────────────────────────────────────────────
        manhattan_click <- mod_manhattan_overlay_server("combined_manhattan",
            project_data         = project_data,
            module               = module,
            combined             = TRUE,
            title_label          = shiny::reactive({
                if (regime_wza()) "Combined GWAS Manhattan (WZA)" else "Combined GWAS Manhattan"
            }),
            note                 = shiny::reactive(htmltools::tagList(
                help_note("manhattan_gwas_combined"),
                no_sig_traits_note(setdiff(all_trait_names_rv(), traits()))
            )),
            regions              = explorer$computed_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = plotted_sigsnps,
            trait_colors         = trait_colors,
            method_shapes        = method_shapes,
            bg_path_override     = combined_wza_bg,
            coords_path_override = combined_wza_coords,
            threshold_y          = combined_threshold_y
        )

        # Manhattan SNP click → select the enclosing region in the explorer
        shiny::observeEvent(manhattan_click(), {
            rid <- manhattan_click()
            if (!is.null(rid) && nzchar(rid)) {
                explorer$selected_region_id(rid)
            }
        }, ignoreNULL = TRUE)

        # ── Per-method tabs ────────────────────────────────────────────────────
        output$method_tabs_ui <- shiny::renderUI({
            ms <- methods()
            if (length(ms) == 0) return(NULL)
            panels <- lapply(ms, function(m) bslib::nav_panel(m, value = m))
            do.call(bslib::navset_card_underline, c(list(id = ns("method_tab")), panels))
        })

        active_method <- shiny::reactive(input$method_tab %||% methods()[1])

        output$per_method_trait_ui <- shiny::renderUI({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(NULL)
            shiny::selectInput(ns("per_method_trait"), "Trait",
                               choices = tr, selected = tr[1], width = "200px")
        })

        per_method_trait <- shiny::reactive({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(NULL)
            input$per_method_trait %||% tr[1]
        })

        # Per-method Manhattan overlay — use interactive sig SNPs (avoids pipeline file read)
        per_method_sigsnps_override <- shiny::reactive({
            m <- active_method(); t <- per_method_trait()
            all_ms <- effective_method_sigsnps()
            if (is.null(m) || is.null(t) || length(all_ms) == 0) return(data.table::data.table())
            dt <- all_ms[[m]]
            if (is.null(dt) || nrow(dt) == 0) return(data.table::data.table())
            dt[trait == t]
        })
        mod_manhattan_overlay_server("method_manhattan",
            project_data         = project_data,
            module               = module,
            method               = active_method,
            trait                = per_method_trait,
            combined             = FALSE,
            sig_snps_override    = per_method_sigsnps_override,
            title_label          = shiny::reactive({
                m <- active_method()
                t <- per_method_trait()
                if (!is.null(m) && !is.null(t)) paste0(m, " — ", t) else "Method Manhattan"
            }),
            note                 = shiny::reactive(help_note("manhattan_gwas_method")),
            show_regions_control = FALSE,
            threshold_y          = per_method_threshold_y
        )

        qq_path <- shiny::reactive({
            m   <- active_method()
            t   <- per_method_trait()
            pd  <- project_data()
            k   <- pd$k_best
            if (is.null(m) || is.null(t) || is.na(k)) return(NULL)
            adj <- resolve_adjust(pd$config, m, MOD_GWAS)
            if (is.null(adj)) return(NULL)
            qq_plot_path(pd$name, module, m, t, k, adj)
        })

        mod_image_card_server("qq_plot",
            path    = qq_path,
            title   = shiny::reactive("QQ Plot"),
            dl_name = shiny::reactive(paste0("qq_", active_method() %||% "method",
                                              "_", per_method_trait() %||% "trait")),
            note    = shiny::reactive(help_note("qq_plot_gwas"))
        )

        # Live counts for the KPI strip (mod_gwas_dashboard_outputs(), 1c6896)
        list(kpis = assoc_live_kpis(
            sigsnps = interactive_sigsnps, regions = explorer$computed_regions,
            traits_hit = traits, traits_all = all_trait_names_rv,
            methods = methods, wza = regime_wza))
    })
}

#' gwas tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_gwas_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "gwas",

        lab_kpi_row(ns, alert_id = NULL),

        shiny::uiOutput(ns("config_badges")),
        shiny::uiOutput(ns("pheno_missing_alert")),

        htmltools::div(
            class = "control-bar lab-gea-controls",
            shiny::uiOutput(ns("threshold_bar")),
            shiny::uiOutput(ns("wza_collapse_note")),
            shiny::uiOutput(ns("filter_bar"))
        ),

        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("combined_manhattan"), height = "100%")
        ),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # LEFT: phenotype map + the two controls that drive only it
            htmltools::div(
                class = "lab-hero-col",
                lab_section_header("Phenotype map", icon = "geo-alt"),
                htmltools::div(
                    class = "control-bar lab-phenomap-bar d-flex align-items-center gap-3",
                    shiny::uiOutput(ns("trait_selector")),
                    bslib::input_switch(ns("points"), "Points", value = FALSE)
                ),
                bslib::card(
                    class = "lab-phenomap-card",
                    bslib::card_header(bsicons::bs_icon("geo-alt"), " Phenotype map"),
                    bslib::card_body(
                        class = "p-2 text-center",
                        htmltools::div(class = "piemap-container",
                                       shiny::uiOutput(ns("phenomap_content")))
                    )
                )
            ),

            # RIGHT: per-method detail
            htmltools::div(
                class = "lab-multiples-col",
                lab_section_header("Per-method detail", icon = "layers"),
                htmltools::div(
                    class = "lab-gea-methodbar d-flex align-items-center gap-3 flex-wrap",
                    shiny::uiOutput(ns("method_tabs_ui")),
                    shiny::uiOutput(ns("per_method_trait_ui"))
                ),
                bslib::layout_columns(
                    col_widths = c(4, 8),
                    fill = FALSE, fillable = FALSE,
                    gap = "0.5rem",
                    htmltools::div(class = "lab-gea-qq",
                                   mod_image_card_ui(ns("qq_plot"))),
                    htmltools::div(class = "lab-gea-methodman",
                                   mod_manhattan_overlay_ui(ns("method_manhattan"),
                                                            height = "100%"))
                )
            )
        ),

        mod_region_explorer_ui(ns("region_explorer"))
    )
}

