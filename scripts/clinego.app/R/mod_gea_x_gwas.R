

#' Overlapping Regions tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_gea_x_gwas_server <- function(id, project_data, run_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_GEAXGWAS

        # ── GEA side data ──────────────────────────────────────────────────────
        gea_methods <- shiny::reactive(find_assoc_methods(project_data()$name, MOD_GEA))

        gea_method_pvalues <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_pvalues_cached(pd$name, MOD_GEA, pd$k_best)
        })
        gea_method_wza_pvalues <- shiny::reactive({
            if (!gea_regime_wza()) return(list())
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_wza_pvalues_cached(pd$name, MOD_GEA, pd$k_best)
        })
        gea_effective_pvalues <- shiny::reactive({
            if (gea_regime_wza()) gea_method_wza_pvalues() else gea_method_pvalues()
        })

        gea_all_trait_names <- shiny::reactive({
            pv <- gea_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID","chr","pos","n_snps","mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        # method -> list(adjust=,threshold=,family=) — pins non-univariate
        # methods (RDA) to their registry rule in the GEA-side matrix unless a
        # cell override exists (see R/fct_threshold_rules.R).
        gea_registry_defaults <- shiny::reactive({
            gea_method_significance_defaults(gea_method_registry())
        })
        gea_config_rules <- shiny::reactive({
            pd <- project_data(); ms <- gea_methods()
            if (length(ms) == 0) return(list())
            stats::setNames(
                lapply(ms, function(m) resolve_adjust(pd$config, m, MOD_GEA) %||% ""), ms)
        })

        gea_effective_sigsnps <- shiny::reactive({
            pd <- project_data()
            compute_method_sigsnps_cached(
                pvalues_list = gea_effective_pvalues(),
                type         = gea_threshold_type(),
                value        = gea_threshold_value(),
                k            = pd$k_best,
                regime       = if (gea_regime_wza()) "wza" else "snp",
                project      = pd$name,
                module       = MOD_GEA,
                cutoffs           = gea_combo_thresholds(),
                overrides         = gea_threshold_overrides(),
                registry_defaults = gea_registry_defaults()
            )
        })

        gea_traits <- shiny::reactive({
            sig <- gea_effective_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        gea_trait_colors  <- shiny::reactive(trait_color_map(gea_all_trait_names()))
        gea_method_shapes <- shiny::reactive(method_shape_map(gea_methods()))

        gea_combo_counts <- shiny::reactive({
            snps <- gea_effective_sigsnps()
            counts <- list()
            for (m in names(snps)) {
                dt <- snps[[m]]
                if (nrow(dt) > 0) {
                    tc <- dt[, .(n = .N), by = "trait"]
                    for (i in seq_len(nrow(tc)))
                        counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
                }
            }
            counts
        })

        # ── GWAS side data ─────────────────────────────────────────────────────
        gwas_methods <- shiny::reactive(find_assoc_methods(project_data()$name, MOD_GWAS))

        gwas_method_pvalues <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_pvalues_cached(pd$name, MOD_GWAS, pd$k_best)
        })
        gwas_method_wza_pvalues <- shiny::reactive({
            if (!gwas_regime_wza()) return(list())
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_wza_pvalues_cached(pd$name, MOD_GWAS, pd$k_best)
        })
        gwas_effective_pvalues <- shiny::reactive({
            if (gwas_regime_wza()) gwas_method_wza_pvalues() else gwas_method_pvalues()
        })

        gwas_all_trait_names <- shiny::reactive({
            pv <- gwas_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID","chr","pos","n_snps","mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        # method -> list(adjust=,threshold=,family=) — GWAS excludes RDA
        # entirely, so this is a no-op fallback in practice.
        gwas_registry_defaults <- shiny::reactive({
            gea_method_significance_defaults(gwas_method_registry())
        })
        gwas_config_rules <- shiny::reactive({
            pd <- project_data(); ms <- gwas_methods()
            if (length(ms) == 0) return(list())
            stats::setNames(
                lapply(ms, function(m) resolve_adjust(pd$config, m, MOD_GWAS) %||% ""), ms)
        })

        gwas_effective_sigsnps <- shiny::reactive({
            pd <- project_data()
            compute_method_sigsnps_cached(
                pvalues_list = gwas_effective_pvalues(),
                type         = gwas_threshold_type(),
                value        = gwas_threshold_value(),
                k            = pd$k_best,
                regime       = if (gwas_regime_wza()) "wza" else "snp",
                project      = pd$name,
                module       = MOD_GWAS,
                cutoffs           = gwas_combo_thresholds(),
                overrides         = gwas_threshold_overrides(),
                registry_defaults = gwas_registry_defaults()
            )
        })

        gwas_traits <- shiny::reactive({
            sig <- gwas_effective_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        gwas_trait_colors  <- shiny::reactive(trait_color_map(gwas_all_trait_names()))
        gwas_method_shapes <- shiny::reactive(method_shape_map(gwas_methods()))

        gwas_combo_counts <- shiny::reactive({
            snps <- gwas_effective_sigsnps()
            counts <- list()
            for (m in names(snps)) {
                dt <- snps[[m]]
                if (nrow(dt) > 0) {
                    tc <- dt[, .(n = .N), by = "trait"]
                    for (i in seq_len(nrow(tc)))
                        counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
                }
            }
            counts
        })

        # File fingerprint reactives — change when the pipeline regenerates pvalue TSVs.
        gea_pvalue_fingerprint <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()
            pd <- project_data()
            pvalues_file_fingerprint(pd$name, MOD_GEA, pd$k_best,
                                     regime = if (gea_regime_wza()) "wza" else "snp")
        })
        gwas_pvalue_fingerprint <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()
            pd <- project_data()
            pvalues_file_fingerprint(pd$name, MOD_GWAS, pd$k_best,
                                     regime = if (gwas_regime_wza()) "wza" else "snp")
        })

        # Per-cell significance thresholds for filter bar cells (bindCached per side)
        gea_combo_thresholds <- shiny::reactive({
            compute_method_thresholds(
                pvalues_list      = gea_effective_pvalues(),
                type              = gea_threshold_type(),
                value             = gea_threshold_value(),
                overrides         = gea_threshold_overrides(),
                registry_defaults = gea_registry_defaults()
            )
        }) |> shiny::bindCache(gea_threshold_type(), gea_threshold_value(),
                               gea_regime_wza(), gea_pvalue_fingerprint(),
                               threshold_overrides_key(gea_threshold_overrides(), gea_registry_defaults()))

        gwas_combo_thresholds <- shiny::reactive({
            compute_method_thresholds(
                pvalues_list      = gwas_effective_pvalues(),
                type              = gwas_threshold_type(),
                value             = gwas_threshold_value(),
                overrides         = gwas_threshold_overrides(),
                registry_defaults = gwas_registry_defaults()
            )
        }) |> shiny::bindCache(gwas_threshold_type(), gwas_threshold_value(),
                               gwas_regime_wza(), gwas_pvalue_fingerprint(),
                               threshold_overrides_key(gwas_threshold_overrides(), gwas_registry_defaults()))

        # ── Per-side clumping distance ──────────────────────────────────────────
        # Each side starts at its saved value, else at its source's own distance (the
        # LD-decay distance unless the config sets bp), as the GEA and GWAS tabs do, and
        # is saved only when moved (clumping_distance_state()).
        gea_snp_clumping_distance <- clumping_distance_state(
            input, project_data, "gea_snp_clumping_distance", module, "gea_snp_clumping_distance", MOD_GEA)
        gwas_snp_clumping_distance <- clumping_distance_state(
            input, project_data, "gwas_snp_clumping_distance", module, "gwas_snp_clumping_distance", MOD_GWAS)

        # Each side's region merge: its config's region_clumping_distance when set, else
        # its clumping distance, as the pipeline merges (acff31).
        .region_merge <- function(src) shiny::reactive({
            pd <- project_data()
            config_region_clumping_distance(pd$config, src, ld_decay_table_path(pd$name))
        })
        gea_region_merge  <- .region_merge(MOD_GEA)
        gwas_region_merge <- .region_merge(MOD_GWAS)

        # ── Threshold rule, strategy and regime per side ──────────────────────
        # Each side's controls start at their saved value, else at their source's
        # config (the strategy at Union), and are saved only when the user moves them
        # (filter_bar_state()). The bars render from these reactives.
        side_keys <- function(prefix) c(
            threshold_type   = paste0(prefix, "threshold_type"),
            threshold_value  = paste0(prefix, "threshold_value"),
            combine_strategy = paste0(prefix, "strategy"),
            regime           = paste0(prefix, "regime"))
        gea_bar  <- filter_bar_state(input, project_data, module, source = MOD_GEA,
                                     registry_defaults = gea_registry_defaults,
                                     prefix = "gea_", rp_keys = side_keys("gea_"))
        gwas_bar <- filter_bar_state(input, project_data, module, source = MOD_GWAS,
                                     registry_defaults = gwas_registry_defaults,
                                     prefix = "gwas_", rp_keys = side_keys("gwas_"))

        gea_regime_wza  <- shiny::reactive(isTRUE(input$gea_regime  %||% gea_bar$regime()))
        gwas_regime_wza <- shiny::reactive(isTRUE(input$gwas_regime %||% gwas_bar$regime()))

        # ── Fill from GEA / GWAS tab buttons ─────────────────────────────────
        # Reads the currently persisted params from the source module's region_params.json
        # keys and pushes them into the prefixed GEAxGWAS filter-bar inputs; each
        # control reports the copied value as a move, so its filter_param_state()
        # saves it.
        # threshold_overrides and the matrix selection have no Shiny input to
        # update*() through (they're reactiveVals inside setup_matrix_rules_server,
        # driven by JS bridges/modals, not plain widgets) — copied directly at
        # the region_params.json level instead, then dest_matrix_rules$reload()
        # forces THIS instance's reactiveVals to pick up the change immediately.
        .fill_from <- function(src_module, prefix, dest_matrix_rules) {
            pd  <- project_data()
            rp  <- read_region_params(pd$name)
            cfg <- pd$config

            src_rule <- default_threshold(cfg, src_module,
                if (identical(src_module, MOD_GWAS)) gwas_registry_defaults() else gea_registry_defaults())
            tt <- get_global_param(rp, src_module, "threshold_type") %||% src_rule$type
            tv <- as.numeric(get_global_param(rp, src_module, "threshold_value") %||% src_rule$value)
            cd <- as.integer(get_global_param(rp, src_module, "snp_clumping_distance") %||%
                      default_clumping_distance(cfg, src_module, ld_decay_table_path(pd$name)))
            rg <- isTRUE(get_global_param(rp, src_module, "regime"))
            st <- get_global_param(rp, src_module, "combine_strategy")
            ov <- get_global_param(rp, src_module, "threshold_overrides")
            sel <- get_global_param(rp, src_module, "snp_matrix_selection")

            shiny::updateSelectInput(session,  paste0(prefix, "threshold_type"),  selected = tt)
            shiny::updateNumericInput(session, paste0(prefix, "threshold_value"), value = tv)
            shiny::updateNumericInput(session, paste0(prefix, "snp_clumping_distance"), value = cd)
            bslib::update_switch(paste0(prefix, "regime"), value = rg, session = session)
            if (!is.null(st))
                shiny::updateRadioButtons(session, paste0(prefix, "combine_strategy"),
                                          selected = .normalize_strategy(st))

            # GEA.configs/GWAS.configs overrides use bare method names (from the
            # source tab, which has no traits of its own to key by here) or
            # "trait::method" pairs — normalize_threshold_overrides() (called by
            # reload()) handles both, migrating bare keys against THIS module's
            # own trait list.
            rp <- set_global_param(rp, module, paste0(prefix, "threshold_overrides"), ov)
            rp <- set_global_param(rp, module, paste0(prefix, "snp_matrix_selection"), sel)
            save_region_params(pd$name, rp)
            dest_matrix_rules$reload()

            shiny::showNotification(
                sprintf("Copied %s tab settings into the GEAxGWAS %s filter bar.", src_module, toupper(prefix)),
                type = "message", duration = 3
            )
        }

        shiny::observeEvent(input$fill_gea,  .fill_from(MOD_GEA,  "gea_",  gea_matrix_rules))
        shiny::observeEvent(input$fill_gwas, .fill_from(MOD_GWAS, "gwas_", gwas_matrix_rules))

        gea_threshold_type_raw  <- shiny::reactive(input$gea_threshold_type  %||% gea_bar$threshold_type())
        gwas_threshold_type_raw <- shiny::reactive(input$gwas_threshold_type %||% gwas_bar$threshold_type())
        # Debounce type so switching to "qval" doesn't fire the full 2M-row scan instantly
        gea_threshold_type  <- shiny::debounce(gea_threshold_type_raw,  500)
        gwas_threshold_type <- shiny::debounce(gwas_threshold_type_raw, 500)

        .make_threshold_value <- function(input_id, side_bar) {
            raw <- shiny::reactive({
                v <- input[[input_id]]
                if (is.null(v) || is.na(v) || v <= 0) side_bar$threshold_value()
                else as.numeric(v)
            })
            shiny::debounce(raw, 500)
        }
        gea_threshold_value  <- .make_threshold_value("gea_threshold_value",  gea_bar)
        gwas_threshold_value <- .make_threshold_value("gwas_threshold_value", gwas_bar)

        # ── Per-cell threshold overrides, matrix selection, rule popup ──────────
        # One call per side — module = MOD_GEAXGWAS (both sides), input_prefix
        # distinguishes "gea_"/"gwas_" so the two popups/reset observers never
        # collide within this one moduleServer. config_state = NULL: neither
        # bar drives a pipeline run, so there's no "Apply rules to config"
        # here — see R/fct_combine.R::setup_matrix_rules_server().
        gea_matrix_rules <- setup_matrix_rules_server(
            input, output, session, ns, input_prefix = "gea_",
            project_data = project_data, module = module,
            methods = gea_methods, all_traits = gea_all_trait_names,
            threshold_type = gea_threshold_type, threshold_value = gea_threshold_value,
            registry_defaults = gea_registry_defaults, config_rules = gea_config_rules,
            config_state = NULL, combo_counts = gea_combo_counts
        )
        gwas_matrix_rules <- setup_matrix_rules_server(
            input, output, session, ns, input_prefix = "gwas_",
            project_data = project_data, module = module,
            methods = gwas_methods, all_traits = gwas_all_trait_names,
            threshold_type = gwas_threshold_type, threshold_value = gwas_threshold_value,
            registry_defaults = gwas_registry_defaults, config_rules = gwas_config_rules,
            config_state = NULL, combo_counts = gwas_combo_counts
        )
        gea_threshold_overrides  <- gea_matrix_rules$overrides
        gwas_threshold_overrides <- gwas_matrix_rules$overrides

        output$gea_threshold_hint <- shiny::renderUI({
            t <- gea_threshold_type()
            hint <- switch(t, bonf="α for Bonferroni", qval="FDR q-value (0–1)",
                           top="Top N SNPs", toppct="Top % of SNPs", custom="Raw p-value (e.g. 1e-5)", "")
            htmltools::span(class="text-muted small mt-1", hint)
        })
        output$gwas_threshold_hint <- shiny::renderUI({
            t <- gwas_threshold_type()
            hint <- switch(t, bonf="α for Bonferroni", qval="FDR q-value (0–1)",
                           top="Top N SNPs", toppct="Top % of SNPs", custom="Raw p-value (e.g. 1e-5)", "")
            htmltools::span(class="text-muted small mt-1", hint)
        })

        # Coerce threshold values when type changes
        shiny::observeEvent(input$gea_threshold_type, {
            v <- input$gea_threshold_value
            if (!threshold_value_valid_for_type(input$gea_threshold_type, v))
                shiny::updateNumericInput(session, "gea_threshold_value",
                    value = threshold_value_default_for_type(input$gea_threshold_type))
        }, ignoreInit = TRUE)
        shiny::observeEvent(input$gwas_threshold_type, {
            v <- input$gwas_threshold_value
            if (!threshold_value_valid_for_type(input$gwas_threshold_type, v))
                shiny::updateNumericInput(session, "gwas_threshold_value",
                    value = threshold_value_default_for_type(input$gwas_threshold_type))
        }, ignoreInit = TRUE)

        # ── Always-present threshold bars ─────────────────────────────────────
        output$gea_threshold_bar <- shiny::renderUI({
            project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    input_prefix          = "gea_",
                    regime_value          = gea_bar$regime(),
                    threshold_type_value  = gea_bar$threshold_type(),
                    threshold_value_value = gea_bar$threshold_value(),
                    regime_context        = "overlap"
                )
            })
        })
        output$gwas_threshold_bar <- shiny::renderUI({
            project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    input_prefix          = "gwas_",
                    regime_value          = gwas_bar$regime(),
                    threshold_type_value  = gwas_bar$threshold_type(),
                    threshold_value_value = gwas_bar$threshold_value(),
                    regime_context        = "overlap"
                )
            })
        })

        # ── Overlap bounds (saved only when moved, filter_param_state()) ──────
        # A file written before the choice values were fixed holds "gea"/"gwas",
        # which is not one of the radio's values any more: .normalize_bounds() maps
        # them onto the canonical names, so restoring never leaves the group unselected.
        bounds_state <- filter_param_state(
            input, project_data, "overlap_bounds", module, "overlap_bounds",
            start_fn = function(pd, saved) {
                tryCatch(.normalize_bounds(saved), error = function(e) NULL) %||% "union"
            },
            coerce = function(v) tryCatch(.normalize_bounds(v), error = function(e) NULL))
        shiny::observe({
            shiny::updateRadioButtons(session, "overlap_bounds", selected = bounds_state())
        })

        overlap_bounds <- shiny::reactive({
            .normalize_bounds(input$overlap_bounds %||% bounds_state())
        })

        # ── GEA filter bar UI (matrix + strategy + clumping) ──────────────────
        output$gea_filter_bar <- shiny::renderUI({
            project_data()
            # isolate()d: seeded from the user's move, never re-rendered by it
            strat_val <- shiny::isolate(gea_bar$strategy())
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = gea_all_trait_names(),
                methods                     = gea_methods(),
                trait_colors                = gea_trait_colors(),
                combo_counts                = gea_combo_counts(),
                combo_thresholds            = gea_combo_thresholds(),
                default_strategy_value      = strat_val,
                snp_clumping_distance_value = gea_snp_clumping_distance(),
                region_clumping_distance_value = gea_region_merge()$bp,
                input_prefix                = "gea_",
                # isolate()d — see setup_matrix_rules_server()
                selected_pairs              = shiny::isolate(gea_matrix_rules$selected_pairs()),
                overrides                   = gea_threshold_overrides(),
                registry_defaults           = gea_registry_defaults(),
                master_type                 = gea_threshold_type(),
                master_value                = gea_threshold_value()
            )
        })

        # ── GWAS filter bar UI (matrix + strategy + clumping) ─────────────────
        output$gwas_filter_bar <- shiny::renderUI({
            project_data()
            # isolate()d: seeded from the user's move, never re-rendered by it
            strat_val <- shiny::isolate(gwas_bar$strategy())
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = gwas_all_trait_names(),
                methods                     = gwas_methods(),
                trait_colors                = gwas_trait_colors(),
                combo_counts                = gwas_combo_counts(),
                combo_thresholds            = gwas_combo_thresholds(),
                default_strategy_value      = strat_val,
                snp_clumping_distance_value = gwas_snp_clumping_distance(),
                region_clumping_distance_value = gwas_region_merge()$bp,
                input_prefix                = "gwas_",
                # isolate()d — see setup_matrix_rules_server()
                selected_pairs              = shiny::isolate(gwas_matrix_rules$selected_pairs()),
                overrides                   = gwas_threshold_overrides(),
                registry_defaults           = gwas_registry_defaults(),
                master_type                 = gwas_threshold_type(),
                master_value                = gwas_threshold_value()
            )
        })

        # ── Interactive sig SNPs (independent per side) ────────────────────────
        gea_interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = gea_effective_sigsnps(),
                tm_selection_json   = input$gea_tm_selection,
                combo_counts        = gea_combo_counts(),
                known_traits        = gea_all_trait_names(),
                strategy            = input$gea_combine_strategy %||% gea_bar$strategy(),
                clumping_distance   = gea_snp_clumping_distance(),
                project_name        = project_data()$name,
                module              = MOD_GEA
            )
        })

        gwas_interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = gwas_effective_sigsnps(),
                tm_selection_json   = input$gwas_tm_selection,
                combo_counts        = gwas_combo_counts(),
                known_traits        = gwas_all_trait_names(),
                strategy            = input$gwas_combine_strategy %||% gwas_bar$strategy(),
                clumping_distance   = gwas_snp_clumping_distance(),
                project_name        = project_data()$name,
                module              = MOD_GWAS
            )
        })

        # ── Independent regions (one set per source) ───────────────────────────
        gea_regions <- shiny::reactive({
            snps <- gea_interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, gea_region_merge()$bp %||% gea_snp_clumping_distance())
        })

        gwas_regions <- shiny::reactive({
            snps <- gwas_interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, gwas_region_merge()$bp %||% gwas_snp_clumping_distance())
        })

        # ── Overlap detection ──────────────────────────────────────────────────
        overlap_pairs <- shiny::reactive({
            compute_region_overlaps(gea_regions(), gwas_regions())
        })

        # ── Overlap regions for the explorer (with bounds strategy applied) ────
        overlap_regions <- shiny::reactive({
            pairs <- overlap_pairs()
            if (is.null(pairs) || nrow(pairs) == 0) return(.empty_overlap_regions())
            compute_all_overlap_regions(
                overlap_pairs  = pairs,
                gea_regions    = gea_regions(),
                gwas_regions   = gwas_regions(),
                strategy       = overlap_bounds(),
                gea_sig_snps   = gea_interactive_sigsnps(),
                gwas_sig_snps  = gwas_interactive_sigsnps()
            )
        })

        # ── Unified sig SNPs for the detail panel (region_id = overlap region_id) ──
        unified_sigsnps <- shiny::reactive({
            gea_snps  <- gea_interactive_sigsnps()
            gwas_snps <- gwas_interactive_sigsnps()

            parts <- list()
            if (!is.null(gea_snps) && nrow(gea_snps) > 0) {
                gea_copy <- data.table::copy(gea_snps)
                gea_copy[, source := "gea"]
                parts[["gea"]] <- gea_copy
            }
            if (!is.null(gwas_snps) && nrow(gwas_snps) > 0) {
                gwas_copy <- data.table::copy(gwas_snps)
                gwas_copy[, source := "gwas"]
                parts[["gwas"]] <- gwas_copy
            }
            if (length(parts) == 0) return(NULL)

            unified <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

            # Assign overlap region_ids so Miami plot click → correct region
            ov_regs <- overlap_regions()
            if (!is.null(ov_regs) && nrow(ov_regs) > 0) {
                unified <- assign_region_ids_from_regions(unified, ov_regs)
            }
            unified
        })

        # ── Region explorer ────────────────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = unified_sigsnps,
            region_distance     = shiny::reactive(1000000L),  # not used when override_regions supplied
            override_regions    = overlap_regions,
            regime              = shiny::reactive(gea_regime_wza() && gwas_regime_wza())
        )

        # ── Miami plot ─────────────────────────────────────────────────────────
        # Combine trait/method maps for both sources (unified legend, stable colors)
        all_traits  <- shiny::reactive(sort(union(gea_all_trait_names(), gwas_all_trait_names())))
        all_methods <- shiny::reactive(union(gea_methods(), gwas_methods()))
        all_trait_colors  <- shiny::reactive(trait_color_map(all_traits()))
        all_method_shapes <- shiny::reactive(method_shape_map(all_methods()))

        # Miami threshold y: min threshold across both sides (most stringent for display)
        # Derived from bindCached combo_thresholds — no redundant full-vector scan.
        miami_threshold_y <- shiny::reactive({
            all_ct <- c(unlist(gea_combo_thresholds()), unlist(gwas_combo_thresholds()))
            all_ct <- all_ct[!is.na(all_ct) & all_ct > 0]
            if (length(all_ct) == 0) return(NULL)
            -log10(min(all_ct))
        })

        miami_wza_bg <- shiny::reactive({
            if (!(gea_regime_wza() && gwas_regime_wza())) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            miami_wza_bg_path(pd$name, k)
        })
        miami_wza_co <- shiny::reactive({
            if (!(gea_regime_wza() && gwas_regime_wza())) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            miami_wza_coords_path(pd$name, k)
        })

        miami_click <- mod_manhattan_overlay_server("miami",
            project_data         = project_data,
            module               = module,
            is_miami             = TRUE,
            combined             = TRUE,
            title_label          = shiny::reactive({
                if (gea_regime_wza() && gwas_regime_wza())
                    "Miami Plot WZA (GEA \u2191 | GWAS \u2193)"
                else
                    "Miami Plot (GEA \u2191 | GWAS \u2193)"
            }),
            note                 = shiny::reactive(help_note("miami_plot")),
            regions              = gea_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = unified_sigsnps,
            trait_colors         = all_trait_colors,
            method_shapes        = all_method_shapes,
            gwas_regions         = gwas_regions,
            overlap_regions      = overlap_regions,
            bg_path_override     = miami_wza_bg,
            coords_path_override = miami_wza_co,
            threshold_y          = miami_threshold_y
        )

        # Miami SNP click → select the enclosing overlap region
        shiny::observeEvent(miami_click(), {
            rid <- miami_click()
            if (!is.null(rid) && nzchar(rid)) {
                # Only select if the clicked region_id is an overlap region
                ov_regs <- overlap_regions()
                if (!is.null(ov_regs) && rid %in% ov_regs$region_id) {
                    explorer$selected_region_id(rid)
                }
                # Clicks on GEA/GWAS-only region IDs are intentionally ignored
                # (they're contextual highlights, not selectable overlap regions)
            }
        }, ignoreNULL = TRUE)

        # ── Pairwise Trait Overlap (on-demand async, driven by filter-bar params) ─
        miami_coords <- shiny::reactive({
            pd   <- project_data()
            k    <- pd$k_best
            if (is.na(k)) return(NULL)
            cp   <- miami_coords_path(pd$name, k)
            fp   <- if (file.exists(cp)) as.character(file.info(cp)$mtime) else "missing"
            load_cached(paste0("coords_miami_", pd$name, "_", k), function() {
                load_coords(cp)
            }, fingerprint = fp)
        })

        # Parameter bundle for pairwise hash — includes a file fingerprint so
        # a pipeline re-run (new pvalue files) invalidates the cache automatically.
        pairwise_params <- shiny::reactive({
            pd <- project_data()
            k  <- pd$k_best
            fp <- pairwise_file_fingerprint(pd$name)
            list(
                k              = k,
                gea_thr_type   = gea_threshold_type(),
                gea_thr_value  = gea_threshold_value(),
                gea_regime     = gea_regime_wza(),
                gwas_thr_type  = gwas_threshold_type(),
                gwas_thr_value = gwas_threshold_value(),
                gwas_regime    = gwas_regime_wza(),
                gea_window     = gea_snp_clumping_distance(),
                gwas_window    = gwas_snp_clumping_distance(),
                min_snps       = as.integer(
                    tryCatch(pd$config$GEAxGWAS$pairwise$min_snps, error = function(e) 2L) %||% 2L
                ),
                fp             = fp
            )
        })

        mod_pairwise_overlap_server("pairwise",
            project_data    = project_data,
            coords          = miami_coords,
            gea_sigsnps     = gea_effective_sigsnps,
            gwas_sigsnps    = gwas_effective_sigsnps,
            pairwise_params = pairwise_params
        )

        # Live counts for the KPI strip (mod_geaxgwas_dashboard_outputs(), 1c6896)
        n_hits <- function(dt) if (is.null(dt) || nrow(dt) == 0) 0L else data.table::uniqueN(dt$SNPID)
        list(kpis = shiny::reactive(list(
            gea_hits     = n_hits(gea_interactive_sigsnps()),
            gea_regions  = nrow(gea_regions()),
            gwas_hits    = n_hits(gwas_interactive_sigsnps()),
            gwas_regions = nrow(gwas_regions()),
            overlaps     = nrow(overlap_regions() %||% data.frame())
        )))
    })
}

#' gea_x_gwas tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_gea_x_gwas_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "geaxgwas",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(6, 6),
            fill = FALSE, fillable = FALSE,
            gap = "0.6rem",

            bslib::card(
                class = "lab-xg-filter lab-xg-gea",
                bslib::card_header(
                    class = "d-flex justify-content-between align-items-center",
                    htmltools::span(bsicons::bs_icon("arrow-up"),
                                    " GEA filter (above the zero line)"),
                    shiny::actionButton(ns("fill_gea"), "Fill from GEA",
                                        class = "btn-sm btn-outline-secondary")
                ),
                bslib::card_body(
                    shiny::uiOutput(ns("gea_threshold_bar")),
                    shiny::uiOutput(ns("gea_filter_bar"))
                )
            ),

            bslib::card(
                class = "lab-xg-filter lab-xg-gwas",
                bslib::card_header(
                    class = "d-flex justify-content-between align-items-center",
                    htmltools::span(bsicons::bs_icon("arrow-down"),
                                    " GWAS filter (below the zero line)"),
                    shiny::actionButton(ns("fill_gwas"), "Fill from GWAS",
                                        class = "btn-sm btn-outline-secondary")
                ),
                bslib::card_body(
                    shiny::uiOutput(ns("gwas_threshold_bar")),
                    shiny::uiOutput(ns("gwas_filter_bar"))
                )
            )
        ),

        # ── HERO: the Miami plot ──────────────────────────────────────────────
        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("miami"), height = "100%")
        ),

        htmltools::div(
            class = "control-bar lab-xg-bounds d-flex align-items-center gap-3",
            htmltools::span(class = "lab-xg-bounds-label", "Overlap bounds"),
            shiny::radioButtons(
                ns("overlap_bounds"), label = NULL,
                # Values are the names .apply_bounds() switches on. They read
                # "gea"/"gwas" until 2026-09-15, which the unnamed switch default
                # silently turned into union bounds.
                choices = c("Union" = "union", "Intersection" = "intersection",
                            "GEA only" = "gea_only", "GWAS only" = "gwas_only"),
                selected = "union", inline = TRUE
            )
        ),

        mod_region_explorer_ui(ns("region_explorer")),

        bslib::accordion(
            id = ns("pairwise_accordion"),
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "pairwise",
                title = htmltools::tagList(bsicons::bs_icon("grid-3x3"),
                                           " Pairwise trait overlap"),
                mod_pairwise_overlap_ui(ns("pairwise"))
            )
        )
    )
}

