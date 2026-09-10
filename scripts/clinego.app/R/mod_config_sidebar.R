#' Config sidebar module
#'
#' A hidden left sidebar that displays and edits pipeline configuration
#' parameters for the active tab. Parameters are split into mandatory
#' (always visible) and optional (in a collapsed accordion).
#'
#' Reads from config_state$working, writes back on input change.
#' Shows a "modified" badge when any field differs from config_state$saved.
#' Provides a Reset button to restore saved values.
#'
#' Phase 2+3: bidirectional sync + run button (runner module is separate).
#'
#' @param id module namespace id
#' @param tab_title display title shown in the sidebar header
#' @noRd
mod_config_sidebar_ui <- function(id, tab_title = "Config", runner_ui = NULL) {
    ns <- shiny::NS(id)

    bslib::sidebar(
        id       = ns("sidebar"),
        title    = htmltools::span(
            class = "d-flex align-items-center gap-2 config-sidebar-title",
            bsicons::bs_icon("sliders2", size = "1em"),
            tab_title
        ),
        open     = "closed",
        width    = "70vw",
        position = "left",
        # Dynamic content rendered by server (fires once on project switch)
        shiny::uiOutput(ns("sidebar_content")),
        # Reset button — shown/hidden via shinyjs
        shinyjs::hidden(
            htmltools::div(
                id    = ns("reset_footer"),
                class = "mt-3 pt-2 border-top",
                shiny::actionButton(
                    ns("reset_btn"),
                    label = htmltools::tagList(
                        bsicons::bs_icon("arrow-counterclockwise", size = "0.85em"),
                        " Reset to saved"
                    ),
                    class = "btn btn-outline-secondary btn-sm w-100 config-reset-btn"
                )
            )
        ),
        # Optional runner UI (Run button + progress) injected from app_ui
        if (!is.null(runner_ui))
            htmltools::div(class = "mt-3 pt-2 border-top", runner_ui)
    )
}

#' @param id module namespace id
#' @param config_state reactiveValues with $saved, $working, $project
#' @param tab_name tab identifier matching config_schema()$tab values
#' @return reactive integer: number of modified fields on this tab
#' @noRd
mod_config_sidebar_server <- function(id, config_state, tab_name, snp_sets_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        entries     <- schema_for_tab(tab_name)
        input_ids   <- setNames(
            gsub("\\.", "_", sapply(entries, `[[`, "key")),
            sapply(entries, `[[`, "key")
        )

        # ── Render UI once per project switch ──────────────────────────────────
        # Uses config_state$working at render time. Does NOT re-render when
        # working changes (that would cause circular update loops).
        # Updates are applied via updateXxxInput() instead.

        last_project <- shiny::reactiveVal(NULL)

        output$sidebar_content <- shiny::renderUI({
            # Depend on project identity, not working config values
            proj <- config_state$project
            last_project(proj)

            cfg <- if (!is.null(proj)) shiny::isolate(config_state$working) else list()

            # Regime filtering happens HERE, on a local copy -- not on `entries`
            # above, which is built once at server init and is what the per-entry
            # observers and input_ids key off. A section dropped here simply never
            # renders; its inputs read NULL, nothing collects them, and the values
            # on disk are left untouched.
            hidden  <- hidden_sidebar_sections(config_regime(cfg))
            entries <- if (length(hidden))
                Filter(function(e) !isTRUE(e$section %in% hidden), entries) else entries

            if (length(entries) == 0)
                return(htmltools::p(class = "text-muted small p-2",
                                    "No configuration parameters for this tab."))

            # Tabs where any entry sets `group` render as flat named top-level
            # sections instead of a single mandatory-core/Advanced split — the
            # mandatory/Advanced split still applies WITHIN each group (nested
            # Advanced accordion) — see config_schema() docs and render_group().
            has_groups <- any(vapply(entries, function(e) !is.null(e$group), logical(1)))

            if (has_groups) {
                groups <- unique(vapply(entries, function(e) e$group %||% "", character(1)))
                return(htmltools::tagList(lapply(groups, function(g) {
                    render_group(ns, g, Filter(function(e) identical(e$group, g), entries), cfg, proj)
                })))
            }

            sp    <- schema_split(entries)
            mand  <- sp$mandatory
            opt   <- sp$optional

            htmltools::tagList(
                if (length(mand) > 0)
                    htmltools::div(
                        class = "config-section",
                        render_section_groups(ns, mand, cfg, proj)
                    ),
                if (length(opt) > 0)
                    bslib::accordion(
                        id   = ns("adv"),
                        open = FALSE,
                        bslib::accordion_panel(
                            title = htmltools::tagList(
                                bsicons::bs_icon("gear", size = "0.9em"),
                                " Advanced"
                            ),
                            value = "advanced",
                            render_section_groups(ns, opt, cfg, proj)
                        )
                    )
            )
        })

        # ── SNP set picker: nested renderUI driven by snp_sets_trigger ──────────
        # Must NOT live inside sidebar_content's renderUI (render-once-per-project).
        # Instead we render a uiOutput placeholder in build_config_input and drive it here.
        picker_entries <- Filter(function(e) e$type == "snp_set_picker", entries)
        if (length(picker_entries) > 0 && !is.null(snp_sets_trigger)) {
            for (pe in picker_entries) {
                local({
                    pe_local <- pe
                    iid      <- gsub("\\.", "_", pe_local$key)
                    dyn_id   <- paste0(iid, "_dynamic")

                    output[[dyn_id]] <- shiny::renderUI({
                        snp_sets_trigger()                        # live refresh on save/delete
                        proj <- config_state$project
                        sets <- list_snp_sets(proj)

                        if (nrow(sets) == 0) {
                            return(htmltools::div(
                                class = "alert alert-warning small mb-0 p-2",
                                bsicons::bs_icon("info-circle"),
                                " No saved SNP sets. Go to the ",
                                htmltools::strong("GEA tab"),
                                " → ",
                                htmltools::strong("Save SNP set for maladaptation"),
                                "."
                            ))
                        }

                        cur <- config_get_by_path(shiny::isolate(config_state$working),
                                                  pe_local$key)
                        if (is.null(cur) || identical(cur, "all")) {
                            selected <- sets$name
                        } else {
                            selected <- intersect(unlist(cur), sets$name)
                            if (length(selected) == 0) selected <- sets$name
                        }

                        shiny::checkboxGroupInput(
                            session$ns(iid), label = NULL,
                            choiceNames  = sprintf("%s (%d SNPs)", sets$name, sets$n_snps),
                            choiceValues = sets$name,
                            selected     = selected
                        )
                    })
                })
            }
        }

        # ── Input observers: write changed values into config_state$working ───
        # Guard: skip if the new value equals the current working value
        # (prevents circular updates from Reset/updateXxxInput calls).

        for (entry in entries) {
            local({
                e        <- entry
                iid      <- gsub("\\.", "_", e$key)
                key_path <- e$key

                if (e$type == "snp_set_picker") {
                    # Dedicated observer — must write "all" or an explicit list, NOT a CSV scalar
                    shiny::observeEvent(input[[iid]], {
                        proj     <- config_state$project
                        all_nms  <- list_snp_sets(proj)$name
                        sel      <- input[[iid]] %||% character(0)
                        # "all" sentinel when every available set is checked; else explicit list
                        new_val <- if (length(all_nms) > 0 && setequal(sel, all_nms)) {
                            "all"
                        } else {
                            as.list(sel)
                        }
                        cur_val <- config_get_by_path(config_state$working, key_path)
                        if (identical(new_val, cur_val)) return()
                        config_state$working <- config_set_by_path(
                            config_state$working, key_path, new_val)
                    }, ignoreInit = TRUE, ignoreNULL = FALSE)
                } else if (e$type == "method_table") {
                    # method_table: observe a hidden JSON bridge input. The JS
                    # row editor collects only {method, params} — adjust/
                    # threshold moved to the interactive matrix's popups
                    # (fct_combine.R). Reconcile them here: carry the existing
                    # working entry's value forward by method name, or seed a
                    # genuinely new method from the registry default. This
                    # server-side reconcile is the ONLY writer of these two
                    # fields other than "Apply rules to config" (mod_gea.R) —
                    # if the JS emitted them itself from stale DOM state, the
                    # next unrelated sidebar edit (add/remove any method)
                    # would silently clobber whatever Apply-to-config had set,
                    # since this editor's DOM is never re-rendered after
                    # project load.
                    json_id         <- paste0(iid, "_json")
                    is_gwas_configs <- identical(key_path, "GWAS.configs")
                    shiny::observeEvent(input[[json_id]], {
                        raw <- input[[json_id]]
                        if (is.null(raw) || nchar(raw) == 0) return()
                        tryCatch({
                            incoming <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
                            # Dedup by method (pipeline enforces one entry per method)
                            if (length(incoming) > 0) {
                                methods_in <- vapply(incoming, function(x) x$method %||% "", character(1))
                                dup_idx <- duplicated(methods_in)
                                if (any(dup_idx)) {
                                    shiny::showNotification(
                                        sprintf("Duplicate method(s) removed: %s",
                                                paste(unique(methods_in[dup_idx]), collapse = ", ")),
                                        type = "warning", duration = 4
                                    )
                                    incoming <- incoming[!dup_idx]
                                }
                            }

                            cur_val      <- config_get_by_path(config_state$working, key_path)
                            registry     <- if (is_gwas_configs) gwas_method_registry() else gea_method_registry()
                            sig_defaults <- gea_method_significance_defaults(registry)
                            cur_by_method <- stats::setNames(
                                cur_val %||% list(),
                                vapply(cur_val %||% list(), function(x) x$method %||% "", character(1))
                            )
                            new_val <- lapply(incoming, function(m) {
                                existing <- cur_by_method[[m$method %||% ""]]
                                sd <- sig_defaults[[m$method]] %||% list(adjust = "bonf", threshold = "0.05")
                                m$adjust    <- existing$adjust    %||% sd$adjust
                                m$threshold <- existing$threshold %||% sd$threshold
                                m
                            })

                            if (identical(new_val, cur_val)) return()
                            config_state$working <- config_set_by_path(
                                config_state$working, key_path, new_val
                            )
                        }, error = function(e) NULL)
                    }, ignoreInit = TRUE, ignoreNULL = TRUE)
                } else {
                    shiny::observeEvent(input[[iid]], {
                        raw     <- input[[iid]]
                        new_val <- input_to_config_value(raw, e$type)
                        cur_val <- config_get_by_path(config_state$working, key_path)

                        # No-op if equal (prevents circular re-triggers)
                        if (config_values_equal(new_val, cur_val)) return()

                        config_state$working <- config_set_by_path(
                            config_state$working, key_path, new_val
                        )
                    }, ignoreInit = TRUE)
                }
            })
        }

        # ── Reset button: restore working from saved for this tab's fields ────
        shiny::observeEvent(input$reset_btn, {
            new_working <- config_state$working
            for (e in entries) {
                saved_val <- config_get_by_path(config_state$saved, e$key)
                new_working <- config_set_by_path(new_working, e$key, saved_val)
            }
            config_state$working <- new_working

            # Push saved values back to inputs (including method tables via JS)
            update_sidebar_inputs(session, entries, config_state$saved)
        })

        # ── Dirty count: compare working vs saved for this tab's fields ───────
        dirty_count <- shiny::reactive({
            working <- config_state$working
            saved   <- config_state$saved
            sum(vapply(entries, function(e) {
                w <- config_get_by_path(working, e$key)
                s <- config_get_by_path(saved,   e$key)
                if (e$type == "method_table") {
                    !identical(w, s)
                } else {
                    !config_values_equal(w, s)
                }
            }, logical(1)))
        })

        # ── Badge + Reset visibility via JS/shinyjs ────────────────────────────
        shiny::observe({
            n <- dirty_count()
            # Send custom message to config-dirty.js
            session$sendCustomMessage("config_dirty_update", list(
                sidebar_id  = ns("sidebar"),
                dirty_count = n
            ))
            # Show/hide reset footer
            if (n > 0) shinyjs::show("reset_footer")
            else       shinyjs::hide("reset_footer")
        })

        # ── Return dirty count for use by pipeline runner ─────────────────────
        return(dirty_count)
    })
}

# ── Internal helpers ───────────────────────────────────────────────────────────

#' Render one top-level config group (used by tabs that opt into `group`,
#' e.g. Processing's "SNP Filtering" / "Climate's "Variance partitioning").
#' Entries with `section == NULL` render flat; entries with a distinct
#' `section` (e.g. "LD Pruning", "Response matrix (Y)") render as a labeled
#' subsection within the group, via the same helper used elsewhere.
#'
#' `mandatory` IS honored here (unlike the old flat-group behavior): any
#' entry with `mandatory = FALSE` is pulled out of the normal section flow
#' and rendered inside a nested, collapsed "Advanced" accordion at the
#' bottom of the group — the same accordion markup used by non-grouped tabs.
#' @noRd
render_group <- function(ns, group_title, entries, config, project = NULL) {
    sp   <- schema_split(entries)
    mand <- sp$mandatory
    opt  <- sp$optional

    htmltools::div(
        class = "config-section",
        htmltools::p(class = "config-section-header", group_title),
        if (length(mand) > 0) render_section_groups(ns, mand, config, project),
        if (length(opt) > 0)
            bslib::accordion(
                id   = ns(paste0("adv_", gsub("[^a-zA-Z0-9]+", "_", group_title))),
                open = FALSE,
                bslib::accordion_panel(
                    title = htmltools::tagList(
                        bsicons::bs_icon("gear", size = "0.9em"),
                        " Advanced"
                    ),
                    value = "advanced",
                    render_section_groups(ns, opt, config, project)
                )
            )
    )
}

#' Render inputs grouped by their section label
#'
#' `section` may be NULL (flat field, no subsection header — used by tabs
#' with `group`, e.g. Processing). NULL is normalized to "" for comparison
#' so it never reaches `==`/`Filter` as a bare NULL (which would error on a
#' zero-length comparison); "" is only ever suppressed from rendering a
#' header, never treated as a distinct blank-titled subsection.
#' @noRd
render_section_groups <- function(ns, entries, config, project = NULL) {
    sec_key <- function(e) e$section %||% ""
    sections <- unique(vapply(entries, sec_key, character(1)))
    multi_section <- length(sections) > 1

    htmltools::tagList(lapply(sections, function(sec) {
        sec_entries <- Filter(function(e) identical(sec_key(e), sec), entries)
        htmltools::tagList(
            if (multi_section && nzchar(sec))
                htmltools::p(class = "config-subsection-label", sec),
            lapply(sec_entries, function(e) render_config_field(ns, e, config, project))
        )
    }))
}

#' Render a single config field (label + input + help)
#' @noRd
render_config_field <- function(ns, entry, config, project = NULL) {
    value    <- config_get_by_path(config, entry$key)
    input_id <- ns(gsub("\\.", "_", entry$key))
    help_el  <- if (!is.null(entry$help))
        htmltools::tags$small(class = "form-text text-muted config-help", entry$help)
    else NULL
    input_el <- build_config_input(input_id, entry, value, project, config)
    label_el <- if (!entry$type %in% c("checkbox", "checkbox_invert"))
        htmltools::tags$label(
            class = "form-label config-field-label",
            `for` = input_id,
            entry$label
        )
    else NULL

    field <- htmltools::div(class = "config-field", label_el, input_el, help_el)

    # Optional reveal-on-check: only render while the referenced checkbox is
    # checked (client-side conditionalPanel, no server round-trip needed).
    if (!is.null(entry$show_if)) {
        cond_id <- gsub("\\.", "_", entry$show_if)
        field <- shiny::conditionalPanel(
            condition = sprintf("input.%s == true", cond_id),
            ns = ns,
            field
        )
    }

    field
}

#' Build the correct Shiny input widget for a schema entry
#' @noRd
build_config_input <- function(input_id, entry, value, project = NULL, config = NULL) {
    type        <- entry$type
    display_val <- normalize_display_value(value, type)

    switch(type,
        "numeric" = shiny::numericInput(
            input_id, label = NULL, value = display_val,
            min = entry$min %||% NA, max = entry$max %||% NA,
            step = entry$step %||% 1, width = "100%"
        ),
        "text" = shiny::textInput(
            input_id, label = NULL, value = display_val,
            placeholder = entry$placeholder %||% "", width = "100%"
        ),
        "select" = shiny::selectInput(
            input_id, label = NULL,
            choices = entry$choices %||% character(0),
            selected = display_val, width = "100%"
        ),
        "checkbox" = shiny::checkboxInput(
            input_id, label = entry$label, value = isTRUE(display_val)
        ),
        "checkbox_invert" = shiny::checkboxInput(
            input_id, label = entry$label, value = !isTRUE(display_val)
        ),
        "textarea" = shiny::textAreaInput(
            input_id, label = NULL, value = display_val,
            rows = 2, width = "100%",
            placeholder = entry$placeholder %||% ""
        ),
        "method_table" = {
            # input_id is already ns()-wrapped; extract the local id for textInput
            # by stripping the namespace prefix (everything up to and incl. last "-")
            local_id <- sub("^.*-", "", input_id)
            # GEA.configs -> GEA_METHODS (includes RDA); GWAS.configs -> GWAS_METHODS
            # (GEA_METHODS filtered to supports_phenotypes=True — excludes RDA).
            is_gea_configs <- identical(entry$key, "GEA.configs")
            registry <- if (identical(entry$key, "GWAS.configs")) gwas_method_registry()
                        else gea_method_registry()
            render_method_editor(input_id, local_id, display_val,
                                registry = registry,
                                k_best = config_k_best(config),
                                # PreGEA only explores climate/GEA traits — its
                                # recommendations don't apply to GWAS's phenotype
                                # method editor, so only wire the project (and
                                # hence the recs lookup) in for GEA.configs.
                                project = if (is_gea_configs) project else NULL)
        },
        "snp_set_picker" = {
            # Render a placeholder uiOutput; actual checkboxGroupInput is rendered
            # server-side and driven by snp_sets_trigger (not sidebar_content).
            # input_id is already ns()-wrapped; append "_dynamic" for the nested output.
            shiny::uiOutput(paste0(input_id, "_dynamic"))
        },
        "bio_chips" = {
            inv <- character(0)
            if (!is.null(project)) {
                inv_path <- climate_invariant_path(project)
                if (file.exists(inv_path)) {
                    inv <- tryCatch(
                        data.table::fread(inv_path, sep = "\t", header = TRUE)$predictor,
                        error = function(e) character(0)
                    )
                }
            }
            render_bio_chips(input_id, display_val, invariant = inv)
        },
        "pheno_chips" = {
            render_pheno_chips(input_id, display_val, project)
        },
        shiny::textInput(input_id, label = NULL,
                          value = as.character(display_val %||% ""), width = "100%")
    )
}

#' Push saved values back to already-rendered inputs (called on Reset)
#' @noRd
update_sidebar_inputs <- function(session, entries, saved_config) {
    for (e in entries) {
        iid <- gsub("\\.", "_", e$key)
        val <- config_get_by_path(saved_config, e$key)

        if (e$type == "method_table") {
            # Reset method editor via JS: overwrite with saved JSON
            # input_id must match the full_json_id used in the editor's JS handler
            json_val <- tryCatch(as.character(jsonlite::toJSON(val %||% list(), auto_unbox = FALSE)),
                                 error = function(e) "[]")
            session$sendCustomMessage("method_editor_reset", list(
                input_id = paste0(session$ns(iid), "_json"),
                configs  = json_val
            ))
            next
        }

        if (e$type == "snp_set_picker") {
            # The nested renderUI re-reads config_state$working on snp_sets_trigger,
            # so just update the working config to the saved value; the renderUI will
            # pick it up on the next tick (or on next snp_sets_trigger bump).
            # No direct updateCheckboxGroupInput needed — the renderUI handles it.
            next
        }

        if (e$type == "bio_chips") {
            session$sendCustomMessage("bio_chips_reset", list(
                container_id = paste0(session$ns(iid), "_container"),
                input_id     = session$ns(iid),
                value        = as.character(val %||% "")
            ))
            next
        }

        if (e$type == "pheno_chips") {
            session$sendCustomMessage("pheno_chips_reset", list(
                container_id = paste0(session$ns(iid), "_container"),
                input_id     = session$ns(iid),
                value        = as.character(val %||% "")
            ))
            next
        }

        dv  <- normalize_display_value(val, e$type)
        switch(e$type,
            "numeric"         = shiny::updateNumericInput(session, iid, value = dv),
            "text"            = shiny::updateTextInput(session, iid, value = dv %||% ""),
            "select"          = shiny::updateSelectInput(session, iid, selected = dv %||% ""),
            "checkbox"        = shiny::updateCheckboxInput(session, iid, value = isTRUE(dv)),
            "checkbox_invert" = shiny::updateCheckboxInput(session, iid, value = !isTRUE(dv)),
            "textarea"        = shiny::updateTextAreaInput(session, iid, value = dv %||% "")
        )
    }
}

#' Convert raw config value to a display-ready scalar
#' @noRd
normalize_display_value <- function(value, type) {
    if (is.null(value)) return(switch(type,
        "numeric"         = NA_real_,
        "checkbox"        = FALSE,
        "checkbox_invert" = FALSE,
        "method_table"    = list(),
        "bio_chips"       = "",
        "pheno_chips"     = "",
        ""
    ))
    switch(type,
        "checkbox" = ,
        "checkbox_invert" = {
            if (is.logical(value)) value
            else value %in% c("T", "TRUE", "true", "1", "yes")
        },
        "textarea" = {
            if (is.list(value) || length(value) > 1)
                paste(value, collapse = ", ")
            else as.character(value)
        },
        "method_table" = value,
        "numeric" = {
            n <- suppressWarnings(as.numeric(value))
            if (is.na(n)) NA_real_ else n
        },
        {
            if (is.logical(value)) as.character(value)
            else if (length(value) > 1) paste(value, collapse = ", ")
            else as.character(value %||% "")
        }
    )
}

#' Compare two config values for equality, tolerating integer/double coercion
#'
#' Sidebar inputs always return doubles for numeric fields, but YAML may
#' store them as integers (e.g. 5L). Using identical() would always report
#' these as dirty. This helper normalizes numeric types before comparing.
#' @noRd
config_values_equal <- function(a, b) {
    if (identical(a, b)) return(TRUE)
    # Numeric type coercion: integer vs double
    if (is.numeric(a) && is.numeric(b)) {
        return(isTRUE(all.equal(as.double(a), as.double(b),
                                tolerance = .Machine$double.eps ^ 0.5)))
    }
    # NULL vs empty string (both mean "not set")
    if ((is.null(a) || identical(a, "")) && (is.null(b) || identical(b, ""))) return(TRUE)
    FALSE
}

#' Recommendation badge/apply-button for a param widget, shared by the R and
#' JS renderers (kept as one function so the two stay in sync by construction).
#' NULL when there is nothing to show (no PreGEA recommendation for this param).
#' @noRd
render_recommendation_badge <- function(spec, value) {
    if (is.null(spec$recommended)) return(NULL)
    matches <- identical(as.character(spec$recommended), as.character(value))
    title <- if (!is.null(spec$evidence)) paste0("PreGEA recommends ", spec$recommended, " (", spec$evidence, ")")
             else paste0("PreGEA recommends ", spec$recommended)
    if (matches) {
        htmltools::tags$span(class = "pregea-match-badge", title = title, "✓ preGEA")
    } else {
        htmltools::tags$button(
            type = "button", class = "pregea-apply-btn", `data-role` = "apply-recommendation",
            `data-apply-value` = as.character(spec$recommended), title = title,
            paste0("→ ", spec$recommended)
        )
    }
}

#' Render one hyperparameter widget (label + input) for a method_table row.
#' @param name param name (e.g. "condition_pcs")
#' @param spec registry param spec: list(type=, default=, min=, max=, choices=, help=,
#'   recommended=, evidence= — the last two set only when PreGEA has a matching
#'   recommendation, see render_method_editor())
#' @param value current resolved value to pre-fill
#' @noRd
render_param_widget <- function(name, spec, value) {
    help_attr <- spec$help %||% ""
    badge <- render_recommendation_badge(spec, value)
    if (identical(spec$type, "enum")) {
        opts <- paste(sapply(spec$choices %||% list(), function(ch) {
            sel <- if (identical(as.character(ch), as.character(value))) " selected" else ""
            sprintf('<option value="%s"%s>%s</option>', ch, sel, ch)
        }), collapse = "")
        htmltools::tags$label(
            class = "method-param",
            name,
            htmltools::tags$select(
                class = "form-select form-select-sm", `data-param` = name,
                `data-param-type` = "enum", title = help_attr,
                htmltools::HTML(opts)
            ),
            badge
        )
    } else if (identical(spec$type, "bool")) {
        input_attrs <- list(type = "checkbox", `data-param` = name,
                            `data-param-type` = "bool", title = help_attr)
        if (isTRUE(value)) input_attrs$checked <- NA
        htmltools::tags$label(
            class = "method-param",
            do.call(htmltools::tags$input, input_attrs),
            name,
            badge
        )
    } else {
        input_type <- if (spec$type %in% c("int", "float")) "number" else "text"
        step <- if (identical(spec$type, "int")) "1" else "any"
        input_attrs <- list(
            type = input_type, class = "form-control form-control-sm",
            `data-param` = name, `data-param-type` = spec$type,
            value = as.character(value %||% ""), step = step, title = help_attr
        )
        if (!is.null(spec$min)) input_attrs$min <- spec$min
        if (!is.null(spec$max)) input_attrs$max <- spec$max
        htmltools::tags$label(
            class = "method-param",
            name,
            do.call(htmltools::tags$input, input_attrs),
            badge
        )
    }
}

#' Inline editable rows for GEA.configs/GWAS.configs method list
#'
#' Renders a compact editor: one row per method (method select, read-only
#' pipeline-rule badge, remove button, collapsed per-method hyperparameter
#' widgets), an "Add method" button, and a hidden textInput JSON bridge that
#' Shiny observers watch.
#'
#' adjust/threshold are NOT editable here — that control moved to the
#' interactive trait×method matrix (fct_combine.R's cell/row/column popups),
#' the only place with per-trait granularity. This editor shows the current
#' PIPELINE rule (what genes/enrichment/selected_snps.tsv will be built with)
#' as a read-only badge; changing it happens via the matrix's "Apply rules to
#' config" button (mod_gea.R), which writes straight into config_state$working.
#' The JS row here never collects adjust/threshold — see the server-side
#' reconcile in mod_config_sidebar_server()'s method_table observer, which
#' carries the existing value forward by method name (or seeds a genuinely
#' new method from the registry default). This is deliberate: if this row
#' collected and re-emitted adjust/threshold itself, they'd hold whatever was
#' true when the page loaded, and the next unrelated sidebar edit (add/remove
#' any method) would silently clobber a value Apply-to-config had since set —
#' this editor's DOM is never re-rendered after project load, so a stale
#' client-side value has no way to self-correct.
#'
#' Registry-driven (registry = gea_method_registry(), read from
#' workflow/methods/gea.py) — the method list AND each method's hyperparameter
#' widgets come entirely from the registry, so adding a new GEA method needs
#' no change here.
#' @noRd
render_method_editor <- function(input_id, local_id, configs_list,
                                 registry = gea_method_registry(), k_best = NULL,
                                 project = NULL) {
    if (is.null(configs_list)) configs_list <- list()

    METHOD_CHOICES <- names(registry)

    # PreGEA's recommended value per (method, param) — same registry param
    # names (LFMM.K, EMMAX.n_pcs, RDA.condition_pcs/axes/predictor_set), since
    # pregea_recommend.R writes against this vocabulary. Merged into the specs
    # below so BOTH the R-rendered rows and the JS-templated (dynamically
    # added/switched) rows show the same pre-fill/Apply badge for free.
    pregea_recs <- if (!is.null(project)) pregea_recommendations_lookup(project) else list()

    # Resolved defaults per method (registry default, sentinel-resolved against
    # k_best) — used both for R-rendered rows and injected into JS for
    # dynamically added/switched rows.
    param_specs_resolved <- Map(function(method, cfg) {
        Map(function(pname, spec) {
            spec$default <- gea_param_default(spec, k_best)
            rec <- pregea_recs[[method]][[pname]]
            if (!is.null(rec)) {
                spec$recommended <- rec$value
                spec$evidence    <- rec$evidence
            }
            spec
        }, names(cfg$params %||% list()), cfg$params %||% list())
    }, names(registry), registry)

    # method -> list(adjust=, threshold=, family=) from the registry's
    # adjust_default/threshold_default/significance_family fields — RDA seeds
    # bonf/0.01, everything else bonf/0.05 (workflow/methods/gea.py). Used
    # only as the DISPLAY fallback for a method with no saved rule yet (the
    # actual reconciliation happens server-side, mirroring this exactly).
    sig_defaults <- gea_method_significance_defaults(registry)

    # Normalise each entry — params merges registry defaults with any saved override.
    configs_norm <- lapply(configs_list, function(m) {
        method <- m$method %||% m$METHOD %||% "EMMAX"
        sd     <- sig_defaults[[method]] %||% list(adjust = "bonf", threshold = "0.05")
        list(
            method    = method,
            adjust    = m$adjust    %||% m$ADJUST    %||% sd$adjust,
            threshold = as.character(m$threshold %||% m$THRESHOLD %||% sd$threshold),
            params    = resolve_row_params(method, m$params %||% m$PARAMS, registry, k_best)
        )
    })

    # Initial JSON for the bridge (used by JS on page load)
    init_json <- tryCatch(jsonlite::toJSON(configs_norm, auto_unbox = TRUE),
                          error = function(e) "[]")

    # Build one HTML row per config entry (pure HTML, no Shiny widgets)
    make_row <- function(idx, m) {
        row_id <- paste0(input_id, "_row_", idx)

        method_opts <- paste(sapply(METHOD_CHOICES, function(ch) {
            sel <- if (ch == m$method) " selected" else ""
            sprintf('<option value="%s"%s>%s</option>', ch, sel, ch)
        }), collapse = "")

        spec_dict <- registry[[m$method]]$params %||% list()
        params_body <- if (length(spec_dict) == 0) {
            htmltools::tags$span(class = "text-muted small", "No parameters")
        } else {
            htmltools::tagList(lapply(names(spec_dict), function(pname) {
                render_param_widget(pname, param_specs_resolved[[m$method]][[pname]],
                                    m$params[[pname]])
            }))
        }

        htmltools::tags$div(
            class = "method-row",
            id    = row_id,
            `data-idx` = idx,
            htmltools::tags$select(
                class    = "form-select form-select-sm method-select",
                `data-role` = "method",
                htmltools::HTML(method_opts)
            ),
            htmltools::tags$span(
                class = "text-muted small ms-1 method-sig-badge",
                title = paste(
                    "Rule the PIPELINE uses for this method (genes, GO enrichment,",
                    "selected_snps.tsv). Set from the GEA/GWAS tab's trait×method",
                    "matrix — \"Apply rules to config\" — or by editing",
                    "config_{project}.yaml directly."
                ),
                sprintf("%s %s", m$adjust, m$threshold)
            ),
            htmltools::tags$button(
                class    = "method-remove-btn",
                type     = "button",
                `data-role` = "remove",
                bsicons::bs_icon("trash3", size = "0.85em")
            ),
            htmltools::tags$details(
                class = "method-params",
                htmltools::tags$summary("Params"),
                htmltools::tags$div(class = "method-params-body", params_body)
            )
        )
    }

    rows <- if (length(configs_norm) > 0)
        mapply(make_row, seq_along(configs_norm), configs_norm, SIMPLIFY = FALSE)
    else
        list()

    full_json_id <- paste0(input_id, "_json")   # ns already applied to input_id

    htmltools::tagList(
        # Hidden JSON bridge — use full namespaced id so Shiny module routing works
        shiny::textInput(full_json_id, label = NULL, value = init_json) |>
            shinyjs::hidden(),

        # Visible editor container
        htmltools::div(
            class           = "method-editor",
            id              = paste0(input_id, "_editor"),
            `data-json-id`  = full_json_id,
            rows,
            htmltools::tags$button(
                class    = "btn btn-outline-secondary btn-sm method-add-btn",
                type     = "button",
                `data-role` = "add",
                bsicons::bs_icon("plus-circle", size = "0.85em"),
                " Add method"
            )
        ),

        # JS: wire up delegated event handlers for this editor.
        # Split into two sprintf() calls concatenated below — sprintf()'s fmt
        # argument has an 8192-char hard limit (R docs), and this one IIFE
        # body (still ONE function scope — the split does not close it, just
        # breaks the R-side string literal) crossed that after the PreGEA
        # recommendation badge/apply-button JS was added.
        htmltools::tags$script(htmltools::HTML(paste0(sprintf(
'(function() {
  var editorId  = "%s";
  var jsonInputId = "%s";
  var PARAM_SPECS = %s;   // method -> {param_name -> {type, default, min, max, choices, help}}

  function escapeAttr(s) {
    return String(s).replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;");
  }

  // Mirrors render_recommendation_badge() (R side) so both renderers agree.
  function recommendationBadgeHTML(spec, value) {
    if (spec.recommended === undefined || spec.recommended === null) return "";
    var matches = (String(spec.recommended) === String(value));
    var title = spec.evidence
      ? ("PreGEA recommends " + spec.recommended + " (" + spec.evidence + ")")
      : ("PreGEA recommends " + spec.recommended);
    if (matches) {
      return "<span class=\\"pregea-match-badge\\" title=\\"" + escapeAttr(title) + "\\">\\u2713 preGEA</span>";
    }
    return "<button type=\\"button\\" class=\\"pregea-apply-btn\\" data-role=\\"apply-recommendation\\" data-apply-value=\\"" +
           escapeAttr(spec.recommended) + "\\" title=\\"" + escapeAttr(title) + "\\">\\u2192 " + spec.recommended + "</button>";
  }

  function widgetHTML(name, spec, value) {
    var help = spec.help ? (" title=\\"" + escapeAttr(spec.help) + "\\"") : "";
    var badge = recommendationBadgeHTML(spec, value);
    if (spec.type === "enum") {
      var opts = (spec.choices || []).map(function(c) {
        var sel = (String(c) === String(value)) ? " selected" : "";
        return "<option value=\\"" + c + "\\"" + sel + ">" + c + "</option>";
      }).join("");
      return "<label class=\\"method-param\\">" + name +
             "<select class=\\"form-select form-select-sm\\" data-param=\\"" + name +
             "\\" data-param-type=\\"enum\\"" + help + ">" + opts + "</select>" + badge + "</label>";
    } else if (spec.type === "bool") {
      var checked = value ? " checked" : "";
      return "<label class=\\"method-param\\"><input type=\\"checkbox\\" data-param=\\"" + name +
             "\\" data-param-type=\\"bool\\"" + checked + help + "> " + name + badge + "</label>";
    } else {
      var inputType = (spec.type === "int" || spec.type === "float") ? "number" : "text";
      var step = (spec.type === "int") ? "1" : "any";
      var minmax = "";
      if (spec.min !== undefined && spec.min !== null) minmax += " min=\\"" + spec.min + "\\"";
      if (spec.max !== undefined && spec.max !== null) minmax += " max=\\"" + spec.max + "\\"";
      var v = (value === undefined || value === null) ? "" : value;
      return "<label class=\\"method-param\\">" + name +
             "<input type=\\"" + inputType + "\\" class=\\"form-control form-control-sm\\" data-param=\\"" +
             name + "\\" data-param-type=\\"" + spec.type + "\\" value=\\"" + v +
             "\\" step=\\"" + step + "\\"" + minmax + help + ">" + badge + "</label>";
    }
  }

  function buildParamsHTML(method, values) {
    var specs = PARAM_SPECS[method] || {};
    var names = Object.keys(specs);
    if (names.length === 0) return "<span class=\\"text-muted small\\">No parameters</span>";
    return names.map(function(n) {
      var v = (values && values[n] !== undefined) ? values[n] : specs[n].default;
      return widgetHTML(n, specs[n], v);
    }).join("");
  }
',
            paste0(input_id, "_editor"),  # %s 1: editorId
            full_json_id,                 # %s 2: jsonInputId (DOM id, namespaced)
            jsonlite::toJSON(param_specs_resolved, auto_unbox = TRUE)  # %s 3: PARAM_SPECS
        ), sprintf(
'
  function collectRows(editor) {
    var rows = [];
    editor.querySelectorAll(".method-row").forEach(function(row) {
      var params = {};
      row.querySelectorAll("[data-param]").forEach(function(w) {
        var name = w.dataset.param;
        var type = w.dataset.paramType || "str";
        var raw;
        if (type === "bool") raw = w.checked;
        else if (type === "int") raw = parseInt(w.value, 10);
        else if (type === "float") raw = parseFloat(w.value);
        else raw = w.value;
        params[name] = raw;
      });
      // adjust/threshold are intentionally NOT collected here — the server
      // reconciles them by method name against the existing working config
      // (see mod_config_sidebar_server()). They are set via the interactive
      // matrix\'s "Apply rules to config" button, never from this row.
      rows.push({
        method: row.querySelector("[data-role=method]").value,
        params: params
      });
    });
    return rows;
  }

  function pushToShiny(editor) {
    var json = JSON.stringify(collectRows(editor));
    var el = document.getElementById(jsonInputId);
    if (el) {
      el.value = json;
      el.dispatchEvent(new Event("change", {bubbles: true}));
      Shiny.setInputValue(jsonInputId, json, {priority: "event"});
    }
  }

  var METHOD_CHOICES = %s;
  var SIG_DEFAULTS   = %s;   // method -> {adjust, threshold, family} — display-only preview

  function sigDefault(method) {
    return SIG_DEFAULTS[method] || {adjust: "bonf", threshold: "0.05"};
  }

  function sigBadgeHTML(method) {
    var sd = sigDefault(method);
    return sd.adjust + " " + sd.threshold;
  }
',
            jsonlite::toJSON(METHOD_CHOICES, auto_unbox = FALSE),  # %s 1: METHOD_CHOICES
            jsonlite::toJSON(sig_defaults, auto_unbox = TRUE)  # %s 2: SIG_DEFAULTS
        ), sprintf(
'
  function makeRow(m) {
    var div = document.createElement("div");
    div.className = "method-row";
    div.setAttribute("data-role-container", "row");
    var def = m ? m.method : METHOD_CHOICES[0];
    var sel = document.createElement("select");
    sel.className = "form-select form-select-sm method-select";
    sel.setAttribute("data-role", "method");
    METHOD_CHOICES.forEach(function(ch) {
      var opt = document.createElement("option");
      opt.value = ch; opt.text = ch;
      if (ch === def) opt.selected = true;
      sel.appendChild(opt);
    });
    div.appendChild(sel);
    var badge = document.createElement("span");
    badge.className = "text-muted small ms-1 method-sig-badge";
    badge.title = "Rule the pipeline uses for this method. Set via the matrix\'s \\"Apply rules to config\\" button, or edit config_{project}.yaml directly.";
    // A row created here is either brand new (m == null) or just had its
    // method switched — either way there is no saved rule for it yet, so
    // the registry default IS the value the server will reconcile it to.
    badge.textContent = sigBadgeHTML(def);
    div.appendChild(badge);
    var btn = document.createElement("button");
    btn.type = "button"; btn.className = "method-remove-btn";
    btn.setAttribute("data-role", "remove");
    btn.innerHTML = "<svg xmlns=\\"http://www.w3.org/2000/svg\\" width=\\"12\\" height=\\"12\\" fill=\\"currentColor\\" viewBox=\\"0 0 16 16\\"><path d=\\"M5.5 5.5A.5.5 0 0 1 6 6v6a.5.5 0 0 1-1 0V6a.5.5 0 0 1 .5-.5zm2.5 0a.5.5 0 0 1 .5.5v6a.5.5 0 0 1-1 0V6a.5.5 0 0 1 .5-.5zm3 .5a.5.5 0 0 0-1 0v6a.5.5 0 0 0 1 0V6z\\"/><path fill-rule=\\"evenodd\\" d=\\"M14.5 3a1 1 0 0 1-1 1H13v9a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V4h-.5a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1H6a1 1 0 0 1 1-1h2a1 1 0 0 1 1 1h3.5a1 1 0 0 1 1 1v1zM4.118 4 4 4.059V13a1 1 0 0 0 1 1h6a1 1 0 0 0 1-1V4.059L11.882 4H4.118zM2.5 3V2h11v1h-11z\\"/></svg>";
    div.appendChild(btn);
    var details = document.createElement("details");
    details.className = "method-params";
    // Auto-open Params for a row with no saved config yet (new or just
    // method-switched) — a saved row from page load goes through the R
    // renderer (make_row() above), which never sets `open`, so those stay
    // collapsed. A project with many methods all pre-expanded is a wall.
    if (!m) details.open = true;
    var summary = document.createElement("summary");
    summary.textContent = "Params";
    details.appendChild(summary);
    var body = document.createElement("div");
    body.className = "method-params-body";
    body.innerHTML = buildParamsHTML(def, m ? m.params : null);
    details.appendChild(body);
    div.appendChild(details);
    return div;
  }

  var editor = document.getElementById(editorId);
  if (!editor) return;

  editor.addEventListener("click", function(e) {
    var role = e.target.closest("[data-role]");
    if (!role) return;
    if (role.getAttribute("data-role") === "remove") {
      role.closest(".method-row").remove();
      pushToShiny(editor);
    } else if (role.getAttribute("data-role") === "add") {
      var addBtn = editor.querySelector("[data-role=add]");
      editor.insertBefore(makeRow(null), addBtn);
      pushToShiny(editor);
    } else if (role.getAttribute("data-role") === "apply-recommendation") {
      // Copy the PreGEA-recommended value into the own input for this param,
      // then drop the button — nothing left to apply once it matches. A
      // stale checkmark badge is not worth re-deriving here; the input value
      // is now the source of truth (visible + persisted via pushToShiny below).
      var label = role.closest(".method-param");
      var input = label ? label.querySelector("[data-param]") : null;
      if (input) {
        var v = role.getAttribute("data-apply-value");
        if (input.tagName === "SELECT") input.value = v;
        else if (input.type === "checkbox") input.checked = (v === "true" || v === "1");
        else input.value = v;
      }
      role.remove();
      pushToShiny(editor);
    }
  });

  editor.addEventListener("change", function(e) {
    // Method switch: rebuild the params block AND the read-only rule badge
    // for THIS row using the new method — a row switched RDA -> LFMM must
    // not carry RDAs axes/condition_pcs params into an LFMM entry (unknown
    // params, rejected at Snakemake parse time, confusingly far from this
    // UI), nor show RDAs bonf/0.01 badge for what the server will actually
    // reconcile as a fresh LFMM entry (registry default: bonf/0.05).
    if (e.target.matches("[data-role=method]")) {
      var row = e.target.closest(".method-row");
      var body = row.querySelector(".method-params-body");
      if (body) body.innerHTML = buildParamsHTML(e.target.value, null);
      var badge = row.querySelector(".method-sig-badge");
      if (badge) badge.textContent = sigBadgeHTML(e.target.value);
      var details = row.querySelector(".method-params");
      if (details) details.open = true;
      pushToShiny(editor);
      return;
    }
    if (e.target.getAttribute("data-role") || e.target.hasAttribute("data-param")) {
      pushToShiny(editor);
    }
  });

  // Reset handler: overwrite editor contents from new JSON (params included —
  // makeRow() already builds the params block from each saved rows m.params).
  Shiny.addCustomMessageHandler("method_editor_reset", function(msg) {
    if (msg.input_id !== "%s") return;
    var editor = document.getElementById(editorId);
    if (!editor) return;
    editor.querySelectorAll(".method-row").forEach(function(r) { r.remove(); });
    var configs = JSON.parse(msg.configs);
    var addBtn = editor.querySelector("[data-role=add]");
    configs.forEach(function(m) { editor.insertBefore(makeRow(m), addBtn); });
    // Update hidden bridge silently (no Shiny event — avoid round-trip)
    var el = document.getElementById(jsonInputId);
    if (el) el.value = msg.configs;
  });
})();
',
            full_json_id                  # %s: method_editor_reset id check
        ))))
    )
}

#' Render a grid of bio_1..bio_19 toggle chips
#'
#' Replaces the plain text input for climate.predictors with clickable pills.
#' A hidden textInput (input_id) carries the comma-separated value for Shiny.
#' @noRd
render_bio_chips <- function(input_id, display_val, invariant = character(0)) {
    all_bios <- paste0("bio_", 1:19)

    # Parse current value; default to all selected when empty/NULL
    selected <- if (!is.null(display_val) && nzchar(trimws(display_val))) {
        trimws(strsplit(as.character(display_val), ",")[[1]])
    } else {
        all_bios
    }
    # Force-remove invariant predictors from selected (they cannot be active)
    selected <- setdiff(selected, invariant)

    chips <- lapply(all_bios, function(b) {
        is_inv <- b %in% invariant
        cls <- if (is_inv) {
            "bc-chip bc-invariant"
        } else if (b %in% selected) {
            "bc-chip bc-active"
        } else {
            "bc-chip"
        }
        htmltools::tags$button(
            type               = "button",
            class              = cls,
            `data-bio`         = b,
            `data-invariant`   = if (is_inv) "true" else NULL,
            title              = if (is_inv) paste0(b, ": no variance in current map extent") else NULL,
            b
        )
    })

    container_id <- paste0(input_id, "_container")

    js_code <- sprintf('
(function() {
    var container = document.getElementById("%s");
    var inputEl   = document.getElementById("%s");
    if (!container || !inputEl) return;

    function syncValue() {
        var active = container.querySelectorAll(".bc-chip.bc-active");
        var vals = Array.from(active).map(function(el) { return el.dataset.bio; });
        var csv = vals.join(",");
        inputEl.value = csv;
        Shiny.setInputValue("%s", csv, {priority: "event"});
    }

    container.addEventListener("click", function(e) {
        var chip = e.target.closest(".bc-chip");
        if (!chip) return;
        if (chip.dataset.invariant === "true") return;
        chip.classList.toggle("bc-active");
        syncValue();
    });
})();
', container_id, input_id, input_id)

    # Wrap chips + hidden bridge
    htmltools::tagList(
        htmltools::div(
            id    = container_id,
            class = "bio-chips"
        ) |> htmltools::tagAppendChildren(.list = chips),
        # Hidden bridge: display:none via CSS class
        htmltools::div(
            class = "bio-chips-bridge",
            shiny::textInput(input_id, label = NULL,
                             value = paste(selected, collapse = ","), width = "100%")
        ),
        htmltools::tags$script(htmltools::HTML(js_code))
    )
}

#' JS message handler for bio_chips_reset (registered in global JS)
#'
#' Called by update_sidebar_inputs() when Reset is clicked.
#' Sent via session$sendCustomMessage("bio_chips_reset", list(...)).
#' The matching handler lives in www/config-dirty.js.
#' @noRd
NULL  # handler registered in www/config-dirty.js

#' Render a grid of phenotype trait toggle chips
#'
#' Analogous to render_bio_chips but for phenotype traits discovered from
#' processing/tables/metadata.tsv columns 5+.
#' Shows a placeholder if the metadata file doesn't exist yet (processing not run).
#' @noRd
render_pheno_chips <- function(input_id, display_val, project) {
    # Discover trait names from processed metadata
    all_traits <- character(0)
    if (!is.null(project)) {
        meta_path <- hap_metadata_path(project)
        if (file.exists(meta_path)) {
            hdr <- tryCatch({
                con <- file(meta_path, "r")
                on.exit(close(con))
                line <- readLines(con, n = 1)
                trimws(strsplit(line, "\t")[[1]])
            }, error = function(e) character(0))
            if (length(hdr) > 4) all_traits <- hdr[5:length(hdr)]
        }
    }

    # No metadata yet — show placeholder
    if (length(all_traits) == 0) {
        return(htmltools::tagList(
            htmltools::div(
                class = "text-muted small fst-italic",
                "Run processing mode first to discover phenotype traits."
            ),
            # Hidden bridge with empty value so field is skipped on collect
            htmltools::div(
                class = "bio-chips-bridge",
                shiny::textInput(input_id, label = NULL, value = "", width = "100%")
            )
        ))
    }

    # Parse current selection; default to all selected when empty/NULL
    selected <- if (!is.null(display_val) && nzchar(trimws(as.character(display_val)))) {
        trimws(strsplit(as.character(display_val), ",")[[1]])
    } else {
        all_traits
    }
    # Remove any saved traits no longer in metadata
    selected <- intersect(selected, all_traits)
    if (length(selected) == 0) selected <- all_traits

    chips <- lapply(all_traits, function(tr) {
        cls <- if (tr %in% selected) "bc-chip bc-active" else "bc-chip"
        htmltools::tags$button(
            type       = "button",
            class      = cls,
            `data-trait` = tr,
            tr
        )
    })

    container_id <- paste0(input_id, "_container")

    js_code <- sprintf('
(function() {
    var container = document.getElementById("%s");
    var inputEl   = document.getElementById("%s");
    if (!container || !inputEl) return;

    function syncValue() {
        var active = container.querySelectorAll(".bc-chip.bc-active");
        var vals = Array.from(active).map(function(el) { return el.dataset.trait; });
        var csv = vals.join(",");
        inputEl.value = csv;
        Shiny.setInputValue("%s", csv, {priority: "event"});
    }

    container.addEventListener("click", function(e) {
        var chip = e.target.closest(".bc-chip");
        if (!chip) return;
        chip.classList.toggle("bc-active");
        syncValue();
    });
})();
', container_id, input_id, input_id)

    htmltools::tagList(
        htmltools::div(
            id    = container_id,
            class = "bio-chips"
        ) |> htmltools::tagAppendChildren(.list = chips),
        htmltools::div(
            class = "bio-chips-bridge",
            shiny::textInput(input_id, label = NULL,
                             value = paste(selected, collapse = ","), width = "100%")
        ),
        htmltools::tags$script(htmltools::HTML(js_code))
    )
}

#' JS message handler for pheno_chips_reset (registered in www/config-dirty.js)
#' @noRd
NULL
