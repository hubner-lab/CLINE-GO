

#' Association tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param module character: MOD_GEA or MOD_GWAS
#' @param config_state reactiveValues from app_server.R ($working/$saved/$project),
#'   or NULL to omit the "Apply rules to config" button (no caller currently
#'   passes NULL, but keeps this module testable/embeddable without it).
#' @noRd
mod_gea_server <- function(id, project_data, run_trigger = NULL, module = MOD_GEA,
                           snp_sets_trigger = NULL, config_state = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Data loading ───────────────────────────────────────────────────────
        methods <- shiny::reactive(find_assoc_methods(project_data()$name, module))

        # Registry (GEA vs GWAS methods differ — GWAS excludes multivariate RDA)
        # — method -> list(adjust=,threshold=,family=). Pins non-univariate
        # methods (RDA) to their registry rule in the trait×method matrix
        # unless a cell override exists (see R/fct_threshold_rules.R).
        method_registry <- shiny::reactive({
            if (identical(module, MOD_GWAS)) gwas_method_registry() else gea_method_registry()
        })
        registry_defaults <- shiny::reactive({
            gea_method_significance_defaults(method_registry())
        })

        # config_rules: METHOD -> resolve_adjust() string ("bonf_0.05") — the rule
        # the PIPELINE FILES on disk (QQ plot, Manhattan background/coords) were
        # actually built with. Independent of the interactive threshold below —
        # shown read-only in the matrix popup so the divergence stays visible.
        config_rules <- shiny::reactive({
            pd <- project_data(); ms <- methods()
            if (length(ms) == 0) return(list())
            stats::setNames(
                lapply(ms, function(m) resolve_adjust(pd$config, m, module) %||% ""), ms)
        })

        # RDA-specific diagnostics/candidates panel — RDA is GEA-only (excluded
        # from GWAS_METHODS via supports_phenotypes=False), so its content
        # naturally renders nothing when this module instance is the GWAS tab.
        mod_rda_details_server("rda_details", project_data = project_data, methods = methods)

        # Full p-value tables (all SNPs, all traits, wide format). Loaded once per
        # project switch and on every successful pipeline run (run_trigger dependency).
        # Kept in reactive memo (not cachem — can exceed 200MB at WGS scale).
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

        # All trait names (threshold-independent; for color stability + no_sig warning)
        all_trait_names_rv <- shiny::reactive({
            pv <- all_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        # Significant traits at current threshold (derived from effective_method_sigsnps)
        traits <- shiny::reactive({
            sig <- effective_method_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        # ── Filter bar ─────────────────────────────────────────────────────────

        # Trait chip colours (Okabe-Ito, keyed on ALL traits for color stability)
        trait_colors <- shiny::reactive({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(character(0))
            trait_color_map(tr)
        })

        # Method shapes (stable across filtering — keyed on full method list)
        method_shapes <- shiny::reactive({
            ms <- methods()
            if (length(ms) == 0) return(character(0))
            method_shape_map(ms)
        })

        # ── Trait x Method combo counts (how many sig SNPs per combo) ─────────
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

        # File fingerprint reactive — changes whenever the pipeline regenerates the
        # pvalue TSVs. Used as part of bindCache keys so stale cached thresholds are
        # automatically invalidated after a pipeline re-run.
        pvalue_fingerprint <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate on pipeline completion
            pd <- project_data()
            pvalues_file_fingerprint(pd$name, module, pd$k_best,
                                     regime = if (regime_wza()) "wza" else "snp")
        })

        # Per-cell significance threshold: "trait::method" -> cutoff p-value (all cells, NA if unavailable)
        # bindCache: the expensive compute_pval_threshold scan runs once per unique
        # (type, value, regime, pipeline-fingerprint) combination, then returns instantly.
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

        # Strategy: defaults to "All"; persisted to region_params.json so the GEAxGWAS
        # "Fill from GEA tab" button can read the user's last-chosen value.
        default_strategy <- shiny::reactive({
            pd    <- project_data()
            rp    <- read_region_params(pd$name)
            saved <- get_global_param(rp, module, "combine_strategy")
            if (!is.null(saved)) .normalize_strategy(saved) else "All"
        })

        shiny::observeEvent(input$combine_strategy, {
            pd <- project_data(); if (is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, module, "combine_strategy", input$combine_strategy)
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        active_strategy <- shiny::reactive(input$combine_strategy %||% default_strategy())

        # ── SNP clumping distance (single param for both merging and overlap) ───
        shiny::observe({
            pd      <- project_data()
            rp      <- read_region_params(pd$name)
            saved_d <- get_global_param(rp, MOD_GEA, "snp_clumping_distance")
            d <- if (!is.null(saved_d)) as.integer(saved_d)
                 else DEFAULT_CLUMPING_DISTANCE
            shiny::updateNumericInput(session, "snp_clumping_distance", value = d)
        })

        shiny::observeEvent(input$snp_clumping_distance, {
            v  <- input$snp_clumping_distance
            pd <- project_data()
            if (is.null(v) || is.na(v) || v < 1000L || is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, MOD_GEA, "snp_clumping_distance", as.integer(v))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        snp_clumping_distance <- shiny::reactive({
            v <- input$snp_clumping_distance
            if (is.null(v) || is.na(v) || v < 1000L) {
                DEFAULT_CLUMPING_DISTANCE
            } else {
                as.integer(v)
            }
        })

        # ── Regime (per-SNP vs WZA) ────────────────────────────────────────────
        shiny::observe({
            pd    <- project_data()
            rp    <- read_region_params(pd$name)
            saved <- get_global_param(rp, MOD_GEA, "regime")
            if (!is.null(saved))
                bslib::update_switch("regime", value = isTRUE(saved), session = session)
        })

        shiny::observeEvent(input$regime, {
            pd <- project_data()
            if (is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, MOD_GEA, "regime", isTRUE(input$regime))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        regime_wza <- shiny::reactive(isTRUE(input$regime))

        # ── WZA collapse note (appears when regime switch is ON) ───────────────
        output$wza_collapse_note <- shiny::renderUI({
            if (!regime_wza()) return(NULL)
            wza_collapse_note(wza_collapse_stats(all_method_wza_pvalues()))
        })

        # ── Threshold type + value (interactive significance level) ───────────
        shiny::observe({
            pd      <- project_data()
            cfg     <- pd$config
            rp      <- read_region_params(pd$name)
            saved_t <- get_global_param(rp, module, "threshold_type")
            saved_v <- get_global_param(rp, module, "threshold_value")
            if (!is.null(saved_t)) {
                shiny::updateSelectInput(session,  "threshold_type",  selected = saved_t)
            } else {
                def <- default_threshold(cfg, module)
                shiny::updateSelectInput(session,  "threshold_type",  selected = def$type)
            }
            if (!is.null(saved_v)) {
                shiny::updateNumericInput(session, "threshold_value", value = as.numeric(saved_v))
            } else if (is.null(saved_t)) {
                def <- default_threshold(cfg, module)
                shiny::updateNumericInput(session, "threshold_value", value = def$value)
            }
        })

        shiny::observeEvent(input$threshold_type, {
            pd <- project_data(); if (is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, module, "threshold_type", input$threshold_type)
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$threshold_value, {
            v  <- input$threshold_value
            pd <- project_data(); if (is.null(pd) || is.null(v) || is.na(v)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, module, "threshold_value", as.numeric(v))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        threshold_type_raw <- shiny::reactive(input$threshold_type %||% "bonf")
        # Debounce type so switching to "qval" doesn't fire the full 2M-row scan instantly
        threshold_type  <- shiny::debounce(threshold_type_raw, 500)
        threshold_value_raw <- shiny::reactive({
            v <- input$threshold_value
            if (is.null(v) || is.na(v) || v <= 0) {
                default_threshold(project_data()$config, module)$value
            } else {
                as.numeric(v)
            }
        })
        # Debounce threshold value to avoid firing qval computation on every keystroke
        threshold_value <- shiny::debounce(threshold_value_raw, 500)

        # ── Per-cell threshold overrides, matrix selection, rule popup ─────────
        # Set via the trait×method matrix's cell/row/column popups (this is
        # the "Apply rules to config" home too) — shared with mod_gwas.R and
        # mod_gea_x_gwas.R so the popup/precedence logic has one implementation.
        # See R/fct_combine.R::setup_matrix_rules_server().
        matrix_rules <- setup_matrix_rules_server(
            input, output, session, ns, input_prefix = "",
            project_data = project_data, module = module,
            methods = methods, all_traits = all_trait_names_rv,
            threshold_type = threshold_type, threshold_value = threshold_value,
            registry_defaults = registry_defaults, config_rules = config_rules,
            config_state = config_state,
            config_module_key = if (identical(module, MOD_GWAS)) "GWAS" else "GEA"
        )
        threshold_overrides <- matrix_rules$overrides

        # Threshold hint text (context-aware help under the value input)
        output$threshold_hint <- shiny::renderUI({
            t <- threshold_type()
            hint <- switch(t,
                bonf   = "α for Bonferroni correction",
                qval   = "FDR q-value target (0–1)",
                top    = "Number of top SNPs per trait",
                custom = "Raw p-value cutoff (e.g. 1e-5)",
                ""
            )
            htmltools::span(class = "text-muted small mt-1", hint)
        })

        # Coerce threshold value when type changes to avoid invalid combos
        shiny::observeEvent(input$threshold_type, {
            v <- input$threshold_value
            if (!threshold_value_valid_for_type(input$threshold_type, v)) {
                shiny::updateNumericInput(
                    session, "threshold_value",
                    value = threshold_value_default_for_type(input$threshold_type)
                )
            }
        }, ignoreInit = TRUE)

        # ── Always-present threshold bar (regime + threshold type/value) ─────
        output$threshold_bar <- shiny::renderUI({
            pd <- project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    regime_value          = isTRUE(input$regime),
                    threshold_type_value  = input$threshold_type  %||% "bonf",
                    threshold_value_value = input$threshold_value %||%
                        default_threshold(pd$config, module)$value,
                    regime_context        = "gea",
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
                # T2.3: pass pre-computed cutoffs from bindCached combo_thresholds so the
                # cold-path threshold scan runs once (inside combo_thresholds) instead of twice.
                cutoffs           = combo_thresholds(),
                overrides         = threshold_overrides(),
                registry_defaults = registry_defaults()
            )
        })

        # ── Threshold y-value for live threshold line on Manhattan ────────────
        # Derived from the bindCached combo_thresholds — no redundant full-vector scan.
        # combined: min raw p threshold across all methods+traits → highest -log10 line
        combined_threshold_y <- shiny::reactive({
            ct <- combo_thresholds()
            if (length(ct) == 0) return(NULL)
            vals <- unlist(ct)
            vals <- vals[!is.na(vals) & vals > 0]
            if (length(vals) == 0) return(NULL)
            -log10(min(vals))  # min raw p = most stringent = highest y line
        })

        # per-method panel: threshold for selected (method, trait)
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
            # Show the strategy this panel ACTUALLY computed with (the interactive
            # filter-bar selector), not GEA.combine_method from the YAML — that key
            # governs the static pipeline tables only and the two can disagree.
            cmb  <- .normalize_strategy(active_strategy())
            # Everything after K mirrors the CURRENT state of the (collapsed)
            # control accordion, so the bar reads as its summary: what the panel
            # computed with, without opening it. All five are live reactives, so
            # the badges move the moment a control does.
            config_badges_bar(
                if (!is.na(k)) config_badge("K", k, "bg-primary"),
                config_badge("strategy", cmb),
                config_badge("significance",
                             format_threshold_rule(threshold_type(), threshold_value())),
                config_badge("clumping", .format_bp(snp_clumping_distance())),
                # Always shown, both ways: the regime is one of two states, and a
                # badge that only appears when WZA is on leaves per-SNP implicit.
                if (regime_wza()) config_badge("regime", "WZA windows", "bg-info")
                else              config_badge("regime", "per-SNP")
            )
        })

        # ── D5: Traits with zero significant SNPs warning ──────────────────────
        # Traits with no significant SNPs/windows at the current threshold.
        # Surfaced as an amber hover badge in the combined Manhattan's card header
        # (see the `note` argument below), not as a full-width alert -- the banner
        # cost ~98px of permanent height for guidance the user acts on rarely.
        missing_traits <- shiny::reactive({
            setdiff(all_trait_names_rv(), traits())
        })

        # ── Filter bar: Trait x Method matrix + Strategy + Clumping ──────────
        output$filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = all_trait_names_rv(),   # full list — always show grid
                methods                     = methods(),
                trait_colors                = trait_colors(),
                combo_counts                = combo_counts(),
                combo_thresholds            = combo_thresholds(),
                default_strategy_value      = default_strategy(),
                snp_clumping_distance_value = snp_clumping_distance(),
                # isolate()d: the matrix's own on/off selection must survive a
                # re-render triggered by combo_counts/combo_thresholds changing
                # (e.g. a master threshold edit) — see setup_matrix_rules_server().
                selected_pairs              = shiny::isolate(matrix_rules$selected_pairs()),
                overrides                   = threshold_overrides(),
                registry_defaults           = registry_defaults(),
                master_type                 = threshold_type(),
                master_value                = threshold_value()
            )
        })

        # ── Interactive sig SNPs (combined + filtered by matrix selection) ─────
        # effective_method_sigsnps() is already threshold-filtered by compute_method_sigsnps_cached.
        # compute_interactive_sigsnps further applies matrix filter, strategy, and clumping.
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

        # ── Save SNP set for maladaptation ────────────────────────────────────
        # Unique-SNP count for the live modal label (long dt → unique by SNPID)
        save_set_n_unique <- shiny::reactive({
            dt <- interactive_sigsnps()
            if (is.null(dt) || nrow(dt) == 0L) return(0L)
            data.table::uniqueN(dt$SNPID)
        })

        shiny::observeEvent(input$save_snp_set, {
            if (save_set_n_unique() == 0L) {
                shiny::showNotification(
                    "No significant SNPs at current threshold/strategy. Adjust filters first.",
                    type = "warning", duration = 5)
                return()
            }
            shiny::showModal(shiny::modalDialog(
                title = "Save SNP set for maladaptation",
                shiny::textInput(ns("snp_set_name"), "Set name",
                                 placeholder = "e.g. bonf05_union"),
                shiny::uiOutput(ns("snp_set_name_feedback")),
                footer = htmltools::tagList(
                    shiny::modalButton("Cancel"),
                    shiny::actionButton(ns("snp_set_save_confirm"), "Save",
                                        class = "btn btn-success")
                ),
                easyClose = TRUE
            ))
        })

        # Live feedback in modal: N SNPs + params summary + format/duplicate validation
        output$snp_set_name_feedback <- shiny::renderUI({
            nm  <- input$snp_set_name %||% ""
            n   <- save_set_n_unique()
            pd  <- shiny::req(project_data())
            valid_fmt <- grepl("^[A-Za-z0-9_.]+$", nm)
            dup       <- nzchar(nm) && valid_fmt && set_exists(pd$name, nm)
            msg <- if (!nzchar(nm)) "Enter a name."
                   else if (!valid_fmt) "Only letters, digits, underscore and dot allowed."
                   else if (dup)       "A set with this name already exists — choose another."
                   else NULL
            n_ov <- length(threshold_overrides())
            htmltools::tagList(
                htmltools::p(class = "small text-muted mt-1",
                    sprintf("%d unique SNPs — threshold %s = %s%s — strategy %s — regime %s — clump %s bp",
                            n,
                            threshold_type(),
                            threshold_value(),
                            if (n_ov > 0) sprintf(" (%d cell%s overridden)", n_ov,
                                                  if (n_ov == 1) "" else "s") else "",
                            active_strategy(),
                            if (regime_wza()) "WZA" else "per-SNP",
                            format(snp_clumping_distance(), big.mark = ","))),
                if (!is.null(msg)) htmltools::div(class = "text-danger small", msg)
            )
        })

        shiny::observeEvent(input$snp_set_save_confirm, {
            nm <- input$snp_set_name %||% ""
            pd <- project_data()
            if (!grepl("^[A-Za-z0-9_.]+$", nm)) {
                shiny::showNotification("Invalid name.", type = "error"); return()
            }
            if (set_exists(pd$name, nm)) {
                shiny::showNotification("A set with this name already exists.", type = "error"); return()
            }
            dt <- interactive_sigsnps()
            if (is.null(dt) || nrow(dt) == 0L) { shiny::removeModal(); return() }

            # Derive traits/methods from WHAT WENT INTO the set (after filter), not the full grid
            params <- list(
                source_module       = module,
                threshold_type      = threshold_type(),
                threshold_value     = threshold_value(),
                threshold_overrides = threshold_overrides(),
                regime              = if (regime_wza()) "wza" else "snp",
                strategy            = active_strategy(),
                clumping_distance   = snp_clumping_distance(),
                traits              = sort(unique(dt$trait)),
                methods             = sort(unique(dt$method))
            )
            n <- save_snp_set(pd$name, nm, dt, params)
            shiny::removeModal()
            shiny::showNotification(
                sprintf("Saved SNP set ‘%s’ (%d SNPs).", nm, n),
                type = "message", duration = 4)
            if (!is.null(snp_sets_trigger))
                snp_sets_trigger(snp_sets_trigger() + 1L)
        })

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
            region_distance     = snp_clumping_distance,
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
                if (regime_wza()) "Combined Manhattan (WZA)" else "Combined Manhattan"
            }),
            note                 = shiny::reactive(htmltools::tagList(
                help_note("manhattan_gea_combined"),
                no_sig_traits_note(missing_traits())
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

        # ── Per-method tabs UI ─────────────────────────────────────────────────
        output$method_tabs_ui <- shiny::renderUI({
            ms <- methods()
            if (length(ms) == 0) return(NULL)
            panels <- lapply(ms, function(m) {
                bslib::nav_panel(m, value = m)
            })
            do.call(bslib::navset_card_underline, c(list(id = ns("method_tab")), panels))
        })

        # Per-method trait selector UI — uses all trait names so user can view any
        # trait regardless of current threshold (useful for QC / comparison)
        output$per_method_trait_ui <- shiny::renderUI({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(NULL)
            shiny::selectInput(ns("per_method_trait"), "Trait",
                               choices = tr, selected = tr[1], width = "200px")
        })

        per_method_trait <- dedupe_reactive(shiny::reactive({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(NULL)
            input$per_method_trait %||% tr[1]
        }))

        # Deduped: both resolve twice on a cold tab -- first the `%||%` fallback,
        # then the real renderUI-backed input carrying the SAME value -- and each
        # repeat forced a full per-method Plotly.newPlot. Measured on SIMDATA at
        # 1680x1050: method_manhattan load renders 3 -> 2. (Deduping the shared
        # upstreams -- regime_wza, the debounced threshold pair,
        # effective_method_sigsnps -- was tried and moved nothing, so it was
        # reverted rather than left in as unproven complexity.)
        active_method <- dedupe_reactive(
            shiny::reactive(input$method_tab %||% methods()[1]))

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
                if (!is.null(m) && !is.null(t)) paste0(m, " — ", t)
                else "Method Manhattan"
            }),
            note                 = shiny::reactive(help_note("manhattan_gea_method")),
            show_regions_control = FALSE,
            threshold_y          = per_method_threshold_y
        )

        # QQ plot for current method + trait
        qq_path <- shiny::reactive({
            m   <- active_method()
            t   <- per_method_trait()
            pd  <- project_data()
            k   <- pd$k_best
            if (is.null(m) || is.null(t) || is.na(k)) return(NULL)
            adj <- resolve_adjust(pd$config, m,
                if (module == MOD_GWAS) MOD_GWAS else MOD_GEA)
            if (is.null(adj)) return(NULL)
            qq_plot_path(pd$name, module, m, t, k, adj)
        })

        mod_image_card_server("qq_plot",
            path    = qq_path,
            title   = shiny::reactive("QQ Plot"),
            dl_name = shiny::reactive({
                paste0("qq_", active_method() %||% "method",
                       "_", per_method_trait() %||% "trait")
            }),
            note    = shiny::reactive(help_note("qq_plot_gea"))
        )

    })
}

#' gea tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_gea_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "gea",

        shiny::uiOutput(ns("config_badges")),

        # ── Controls: collapsed by default, above the plot they drive ─────────
        # Safe to collapse: mod_gea_server() has exactly one req() in 729 lines
        # (on project_data(), inside the save observer), and every control input
        # it reads is guarded by `%||%` fallback -- input$threshold_type %||%
        # "bonf", input$combine_strategy %||% default_strategy(), and so on. A
        # suspended control bar therefore yields pipeline defaults rather than a
        # blocked plot. mod_gea_dashboard_outputs() additionally un-suspends them so
        # the values saved in region_params.json still reach them.
        bslib::accordion(
            id    = ns("params_accordion"),
            class = "lab-gea-params",
            open  = FALSE,
            bslib::accordion_panel(
                value = "params",
                title = htmltools::tagList(
                    bsicons::bs_icon("sliders"),
                    htmltools::span(" Significance, methods & strategy")
                ),
                htmltools::div(
                    class = "lab-gea-controls",
                    shiny::uiOutput(ns("threshold_bar")),
                    shiny::uiOutput(ns("wza_collapse_note")),
                    shiny::uiOutput(ns("filter_bar")),
                    htmltools::div(
                        class = "d-flex justify-content-end mt-1",
                        shiny::actionButton(
                            ns("save_snp_set"),
                            label = htmltools::tagList(
                                bsicons::bs_icon("save"),
                                " Save SNP set for maladaptation"
                            ),
                            class = "btn-sm btn-success"
                        )
                    )
                )
            )
        ),

        # ── HERO: combined Manhattan, full width, no injected filter_ui ───────
        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("combined_manhattan"), height = "100%")
        ),

        # ── Per-method detail ─────────────────────────────────────────────────
        lab_section_header("Per-method detail", icon = "layers"),
        htmltools::div(
            class = "lab-gea-methodbar d-flex align-items-center gap-3 flex-wrap",
            shiny::uiOutput(ns("method_tabs_ui")),
            shiny::uiOutput(ns("per_method_trait_ui"))
        ),
        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",
            # QQ is square-ish -> the narrow column; the per-method Manhattan is
            # wide -> the broad one. Inverted vs the display modules on purpose.
            htmltools::div(class = "lab-gea-qq", mod_image_card_ui(ns("qq_plot"))),
            htmltools::div(class = "lab-gea-methodman",
                           mod_manhattan_overlay_ui(ns("method_manhattan"),
                                                    height = "100%"))
        ),

        # ── RDA diagnostics: collapsed (suspended until opened) ───────────────
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "rda",
                title = htmltools::tagList(bsicons::bs_icon("compass"),
                                           " RDA diagnostics & candidates"),
                mod_rda_details_ui(ns("rda_details"))
            )
        ),

        # ── Region Explorer: unchanged, full width, bottom ────────────────────
        mod_region_explorer_ui(ns("region_explorer"))
    )
}

