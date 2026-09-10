#' Pipeline runner module
#'
#' Provides a Run button, progress bar, collapsible rule list, and log viewer
#' for each tab. When the user clicks Run:
#'   1. Working config is written to YAML (via fct_config_writer.R)
#'   2. Snakemake is launched as a background processx process
#'   3. Progress is polled from .progress.jsonl every second
#'   4. Rule list is shown as rules start/finish (expandable toggle)
#'   5. On success: bar turns solid green at 100%; config_state$saved updated
#'   6. On failure: bar turns solid red; failed rules highlighted in list
#'   7. Click any completed rule -> open its log in a modal
#'   8. Click the progress bar (after run) -> open main .run.log in a modal
#'
#' The module uses a shared global reactiveVal `pipeline_running` to prevent
#' concurrent runs across all tabs.
#'
#' @param id module namespace id
#' @param tab_name tab identifier (used to derive the Snakemake mode)
#' @param mode_override character — override the mode derived from tab_name.
#'   Used for the haplotype tab where two modes are available.
#' @param btn_label label shown on the Run button (default: "Run")
#' @noRd
mod_pipeline_runner_ui <- function(id, btn_label = NULL, mode_override = NULL) {
    ns <- shiny::NS(id)

    htmltools::div(
        class = "pipeline-runner",

        # Row 1: Run / Stop buttons + status badge
        htmltools::div(
            class = "d-flex align-items-center gap-2 flex-wrap",

            shiny::actionButton(
                ns("run_btn"),
                label = htmltools::tagList(
                    bsicons::bs_icon("play-circle", size = "0.9em"),
                    " ", if (!is.null(btn_label)) btn_label else "Run"
                ),
                class = "btn btn-primary btn-sm pipeline-run-btn"
            ),

            shinyjs::hidden(
                shiny::actionButton(
                    ns("stop_btn"),
                    label = htmltools::tagList(
                        bsicons::bs_icon("stop-circle", size = "0.9em"),
                        " Stop"
                    ),
                    class = "btn btn-danger btn-sm pipeline-stop-btn",
                    id    = ns("stop_btn")
                )
            ),

            shiny::uiOutput(ns("status_badge"))
        ),

        # Row 2: Progress bar + label + expand toggle (hidden until first run)
        shinyjs::hidden(
            htmltools::div(
                id    = ns("progress_container"),
                class = "pipeline-progress-container w-100 mt-2",

                # Progress bar (clickable after completion to open .run.log)
                htmltools::div(
                    class   = "progress pipeline-progress-bar",
                    style   = "height: 8px; cursor: pointer;",
                    onclick = sprintf(
                        "Shiny.setInputValue('%s', Date.now(), {priority: 'event'});",
                        ns("progress_bar_click")
                    ),
                    htmltools::div(
                        id    = ns("progress_bar_inner"),
                        class = "progress-bar progress-bar-striped progress-bar-animated bg-success",
                        role  = "progressbar",
                        style = "width: 0%"
                    )
                ),

                # Label row: "32/50 rules -- rule_name" + chevron toggle
                htmltools::div(
                    class = "d-flex align-items-center justify-content-between mt-1",
                    shiny::uiOutput(ns("progress_label")),
                    shiny::actionLink(
                        ns("toggle_rules"),
                        label = bsicons::bs_icon("chevron-down", size = "0.8em"),
                        class = "pipeline-toggle-rules text-muted"
                    )
                )
            )
        ),

        # Row 3: Warning summary (hidden until warnings found after completion)
        shinyjs::hidden(
            htmltools::div(
                id    = ns("warning_container"),
                class = "pipeline-warnings mt-2",
                shiny::uiOutput(ns("warning_summary"))
            )
        ),

        # Row 4: Collapsible rule list (hidden by default)
        shinyjs::hidden(
            htmltools::div(
                id    = ns("rule_list_container"),
                class = "pipeline-rule-list mt-1",
                shiny::uiOutput(ns("rule_list"))
            )
        )
    )
}

#' @param id module namespace id
#' @param config_state reactiveValues from app_server ($saved, $working, $project)
#' @param pipeline_running reactiveVal (logical) — global run lock
#' @param project_data_trigger reactiveVal — increment to force project_data refresh
#' @param tab_name tab name for mode derivation
#' @param mode_override character — override derived mode (for haplotype)
#' @noRd
mod_pipeline_runner_server <- function(id, config_state, pipeline_running,
                                        project_data_trigger,
                                        tab_name, mode_override = NULL,
                                        snp_sets_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Derive mode
        mode <- if (!is.null(mode_override)) mode_override else tab_to_mode(tab_name)

        # ── Maladaptation: gate Run button on presence of saved SNP sets ─────
        if (identical(tab_name, "maladaptation") && !is.null(snp_sets_trigger)) {
            shiny::observe({
                snp_sets_trigger()
                proj   <- config_state$project
                n_sets <- if (is.null(proj)) 0L else nrow(list_snp_sets(proj))
                if (n_sets == 0L) shinyjs::disable("run_btn")
                else              shinyjs::enable("run_btn")
            })
        }

        # Internal state
        proc             <- shiny::reactiveVal(NULL)
        run_status       <- shiny::reactiveVal("idle")  # idle | running | success | failed
        poll_timer       <- shiny::reactiveTimer(1000)
        rules_expanded   <- shiny::reactiveVal(FALSE)
        selected_log     <- shiny::reactiveVal(NULL)
        cached_progress  <- shiny::reactiveVal(list()) # last read progress; updated only while running
        shown_warnings   <- shiny::reactiveVal(integer(0)) # line numbers of warnings already toasted

        # ── Run button ──────────────────────────────────────────────────────────
        shiny::observeEvent(input$run_btn, {
            shiny::req(!pipeline_running())
            shiny::req(!is.null(config_state$project))
            if (is.null(mode)) {
                shiny::showNotification(
                    "Select a specific mode for this tab (Scan or Haplotype).",
                    type = "warning"
                )
                return()
            }

            project  <- config_state$project
            cpu      <- config_get_by_path(config_state$working, "cpu") %||% max(1L, parallel::detectCores() - 2L)
            pip_path <- get_pipeline_path()

            # ── Maladaptation: hard guard + "all" injection ─────────────────
            if (identical(tab_name, "maladaptation")) {
                if (nrow(list_snp_sets(project)) == 0L) {
                    shiny::showNotification(
                        "No saved SNP sets. Go to the GEA tab and save at least one set first.",
                        type = "warning", duration = 6)
                    return()
                }
                # Inject "all" sentinel if snp_sets is absent in working config
                if (is.null(config_get_by_path(config_state$working, "Maladaptation.snp_sets"))) {
                    config_state$working <- config_set_by_path(
                        config_state$working, "Maladaptation.snp_sets", "all")
                }
            }

            # Write working config to YAML before launching
            ok <- write_project_config(config_state$working, project, pip_path)
            if (is.null(ok)) {
                shiny::showNotification(
                    "Failed to write config file. Check file permissions.",
                    type = "error", duration = 8
                )
                return()
            }

            config_state$saved <- config_state$working

            # Notify if a stale lock will be auto-cleared (pipeline_run handles the actual unlock)
            if (pipeline_is_locked(project, pip_path)) {
                shiny::showNotification(
                    "Cleared a stale lock from a previously interrupted run. Resuming now.",
                    type = "warning", duration = 6
                )
            }

            p <- tryCatch(
                pipeline_run(mode, project, as.integer(cpu), pip_path),
                error = function(e) {
                    shiny::showNotification(
                        paste("Could not start pipeline:", conditionMessage(e)),
                        type = "error", duration = 10
                    )
                    NULL
                }
            )
            if (is.null(p)) return()

            proc(p)
            pipeline_running(TRUE)
            run_status("running")

            # Reset rule list + warning state for new run
            rules_expanded(FALSE)
            shown_warnings(integer(0))
            shinyjs::hide("rule_list_container")
            shinyjs::hide("warning_container")

            # Reset progress bar to initial animated state
            shinyjs::runjs(sprintf(
                'var b = document.getElementById("%s");
                 b.className = "progress-bar progress-bar-striped progress-bar-animated bg-success";
                 b.style.width = "0%%";',
                ns("progress_bar_inner")
            ))

            shinyjs::show("stop_btn")
            shinyjs::show("progress_container")
            shinyjs::disable("run_btn")
        })

        # ── Stop button ─────────────────────────────────────────────────────────
        shiny::observeEvent(input$stop_btn, {
            p <- proc()
            if (!is.null(p)) {
                pipeline_kill(p, config_state$project, mode)
            }
            proc(NULL)
            pipeline_running(FALSE)
            run_status("idle")
            shinyjs::hide("stop_btn")
            shinyjs::hide("progress_container")
            shinyjs::hide("rule_list_container")
            rules_expanded(FALSE)
            shinyjs::enable("run_btn")
            shiny::showNotification("Pipeline stopped.", type = "warning")
        })

        # ── Poll for progress ───────────────────────────────────────────────────
        shiny::observe({
            poll_timer()
            p <- proc()
            shiny::req(pipeline_is_running(p))

            prog <- pipeline_read_progress(config_state$project)
            cached_progress(prog)  # update cache; triggers label + rule_list renders

            pct <- if (!is.na(prog$total) && prog$total > 0)
                round(100 * (prog$done + prog$failed) / prog$total)
            else 0L

            bar_class <- if (prog$failed > 0)
                "progress-bar progress-bar-striped progress-bar-animated bg-danger"
            else
                "progress-bar progress-bar-striped progress-bar-animated bg-success"

            shinyjs::runjs(sprintf(
                'var b = document.getElementById("%s");
                 b.style.width = "%d%%";
                 b.className = "%s";',
                ns("progress_bar_inner"), pct, bar_class
            ))

            # Toast new warnings as they appear in the log
            warns <- pipeline_read_warnings(config_state$project)
            if (nrow(warns) > 0) {
                new_w <- warns[!warns$line %in% shown_warnings(), , drop = FALSE]
                for (i in seq_len(nrow(new_w))) {
                    shiny::showNotification(
                        htmltools::tagList(
                            bsicons::bs_icon("exclamation-triangle", size = "0.9em"), " ",
                            new_w$message[i]
                        ),
                        type = "warning", duration = 8
                    )
                }
                if (nrow(new_w) > 0) shown_warnings(c(shown_warnings(), new_w$line))
            }
        })

        # ── Detect completion ───────────────────────────────────────────────────
        shiny::observe({
            poll_timer()
            p <- proc()
            if (is.null(p)) return()
            if (pipeline_is_running(p)) return()

            exit_code <- tryCatch(p$get_exit_status(), error = function(e) -1L)
            proc(NULL)
            pipeline_running(FALSE)
            # Final progress read — updates label + rule list one last time, then stops
            cached_progress(pipeline_read_progress(config_state$project))
            shinyjs::hide("stop_btn")
            shinyjs::enable("run_btn")

            if (identical(exit_code, 0L)) {
                run_status("success")
                # Solid green bar at 100%
                shinyjs::runjs(sprintf(
                    'var b = document.getElementById("%s");
                     b.style.width = "100%%";
                     b.className = "progress-bar bg-success";',
                    ns("progress_bar_inner")
                ))
                shiny::showNotification(
                    paste0("Pipeline completed: mode=", mode, "."),
                    type = "message", duration = 5
                )
                project_data_trigger(project_data_trigger() + 1L)
            } else {
                run_status("failed")
                # Solid red bar
                shinyjs::runjs(sprintf(
                    'var b = document.getElementById("%s");
                     b.className = "progress-bar bg-danger";',
                    ns("progress_bar_inner")
                ))
                shiny::showNotification(
                    paste0("Pipeline failed (exit ", exit_code, "). ",
                           "Click the progress bar to view the log."),
                    type = "error", duration = 10
                )
            }
            # Keep progress_container visible after completion so user can review

            # Show warning summary if any warnings were found
            warns <- pipeline_read_warnings(config_state$project)
            if (nrow(warns) > 0) {
                shinyjs::show("warning_container")
            }
        })

        # ── Disable run when another tab is running ─────────────────────────────
        shiny::observe({
            if (pipeline_running() && run_status() != "running") {
                shinyjs::disable("run_btn")
            } else if (!pipeline_running() && run_status() != "running") {
                shinyjs::enable("run_btn")
            }
        })

        # ── Hide progress when project changes ──────────────────────────────────
        shiny::observeEvent(config_state$project, {
            run_status("idle")
            shinyjs::hide("progress_container")
            shinyjs::hide("rule_list_container")
            shinyjs::hide("warning_container")
            rules_expanded(FALSE)
            shown_warnings(integer(0))
        }, ignoreInit = TRUE)

        # ── Toggle rule list ────────────────────────────────────────────────────
        shiny::observeEvent(input$toggle_rules, {
            expanded <- !rules_expanded()
            rules_expanded(expanded)
            if (expanded) {
                shinyjs::show("rule_list_container")
                shinyjs::runjs(sprintf(
                    'var el = document.getElementById("%s");
                     if (el) { var icon = el.querySelector("svg, i");
                     el.innerHTML = "";
                     el.appendChild(document.createElementNS("http://www.w3.org/2000/svg","svg")); }',
                    ns("toggle_rules")
                ))
            } else {
                shinyjs::hide("rule_list_container")
            }
            # Swap chevron icon direction via class toggle
            shinyjs::runjs(sprintf(
                '(function() {
                    var btn = document.getElementById("%s");
                    if (!btn) return;
                    var svg = btn.querySelector("svg");
                    if (!svg) return;
                    var use = svg.querySelector("use");
                    if (!use) return;
                    var href = use.getAttribute("href") || use.getAttribute("xlink:href") || "";
                    if (href.indexOf("chevron-down") !== -1) {
                        use.setAttribute("href", href.replace("chevron-down","chevron-up"));
                    } else {
                        use.setAttribute("href", href.replace("chevron-up","chevron-down"));
                    }
                })();',
                ns("toggle_rules")
            ))
        })

        # ── Status badge ────────────────────────────────────────────────────────
        output$status_badge <- shiny::renderUI({
            st <- run_status()
            switch(st,
                idle    = NULL,
                running = htmltools::span(
                    class = "badge bg-warning text-dark pipeline-status-badge",
                    bsicons::bs_icon("arrow-repeat", size = "0.8em"),
                    " Running..."
                ),
                success = htmltools::span(
                    class = "badge bg-success pipeline-status-badge",
                    bsicons::bs_icon("check-circle", size = "0.8em"),
                    " Done"
                ),
                failed  = htmltools::span(
                    class = "badge bg-danger pipeline-status-badge",
                    bsicons::bs_icon("x-circle", size = "0.8em"),
                    " Failed"
                )
            )
        })

        # ── Progress label ──────────────────────────────────────────────────────
        output$progress_label <- shiny::renderUI({
            prog <- cached_progress()
            st   <- run_status()
            if (st == "idle" || length(prog) == 0) return(NULL)
            total_str <- if (!is.na(prog$total)) as.character(prog$total) else "?"
            n_done <- prog$done + prog$failed
            rule_str <- if (st == "running" && nzchar(prog$current_rule))
                paste0(" \u2014 ", prog$current_rule) else ""
            htmltools::p(
                class = "pipeline-progress-label small text-muted mb-0",
                paste0(n_done, "/", total_str, " rules", rule_str)
            )
        })

        # ── Warning summary (shown after run if warnings exist) ─────────────────
        output$warning_summary <- shiny::renderUI({
            st <- run_status()
            if (!st %in% c("success", "failed")) return(NULL)

            warns <- pipeline_read_warnings(config_state$project)
            if (nrow(warns) == 0) return(NULL)

            htmltools::tagList(
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-1",
                    bsicons::bs_icon("exclamation-triangle-fill",
                                     class = "text-warning", size = "0.85em"),
                    htmltools::strong(
                        class = "small",
                        paste0(nrow(warns), " warning", if (nrow(warns) > 1) "s" else "",
                               " during run")
                    )
                ),
                htmltools::tags$ul(
                    class = "mb-0 ps-3",
                    lapply(warns$message, function(m) {
                        htmltools::tags$li(class = "small", m)
                    })
                )
            )
        })

        # ── Rule list ───────────────────────────────────────────────────────────
        output$rule_list <- shiny::renderUI({
            prog <- cached_progress()
            st   <- run_status()
            if (st == "idle" || length(prog) == 0) return(NULL)

            if (length(prog$rules) == 0) {
                return(htmltools::p(class = "small text-muted px-2 mb-0",
                                    "Waiting for first rule..."))
            }

            items <- lapply(prog$rules, function(r) {
                icon <- switch(r$status,
                    done    = bsicons::bs_icon("check-circle-fill",
                                               class = "text-success flex-shrink-0",
                                               size  = "0.85em"),
                    failed  = bsicons::bs_icon("x-circle-fill",
                                               class = "text-danger flex-shrink-0",
                                               size  = "0.85em"),
                    running = bsicons::bs_icon("arrow-repeat",
                                               class = "text-warning spin-icon flex-shrink-0",
                                               size  = "0.85em")
                )

                has_log <- !is.null(r$log_file) && !is.na(r$log_file) &&
                           nzchar(r$log_file)

                name_el <- if (has_log) {
                    safe_path <- gsub("'", "\\\\'", r$log_file)
                    safe_rule <- gsub("'", "\\\\'", r$rule)
                    htmltools::span(
                        class   = "pipeline-rule-link small",
                        style   = "cursor:pointer;",
                        onclick = sprintf(
                            "Shiny.setInputValue('%s', {path:'%s',rule:'%s',ts:Date.now()},{priority:'event'});",
                            ns("rule_click"), safe_path, safe_rule
                        ),
                        r$rule
                    )
                } else {
                    htmltools::span(class = "small text-muted", r$rule)
                }

                htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2",
                    icon,
                    name_el
                )
            })

            htmltools::div(class = "pipeline-rule-items", items)
        })

        # ── Rule log click -> modal ─────────────────────────────────────────────
        shiny::observeEvent(input$rule_click, {
            info <- input$rule_click
            shiny::req(!is.null(info$path) && nzchar(info$path))

            selected_log(info$path)
            log_content <- pipeline_read_log_tail(info$path, n = 200L)

            shiny::showModal(shiny::modalDialog(
                title     = htmltools::tagList(
                    bsicons::bs_icon("file-text", size = "1em"), " ",
                    paste("Log:", info$rule)
                ),
                size      = "l",
                easyClose = TRUE,
                htmltools::pre(
                    class = "pipeline-log-content",
                    paste(log_content, collapse = "\n")
                ),
                footer = htmltools::tagList(
                    shiny::downloadButton(
                        ns("download_log"),
                        label = htmltools::tagList(
                            bsicons::bs_icon("download", size = "0.85em"), " Download"
                        ),
                        class = "btn btn-outline-secondary btn-sm"
                    ),
                    shiny::modalButton("Close")
                )
            ))
        })

        # ── Progress bar click -> main .run.log modal ───────────────────────────
        shiny::observeEvent(input$progress_bar_click, {
            st <- run_status()
            shiny::req(st %in% c("success", "failed"))

            log_path <- pipeline_run_log_path(config_state$project)
            selected_log(log_path)
            log_content <- pipeline_read_log_tail(log_path, n = 300L)

            title_icon <- if (st == "success")
                bsicons::bs_icon("check-circle", class = "text-success")
            else
                bsicons::bs_icon("x-circle", class = "text-danger")

            shiny::showModal(shiny::modalDialog(
                title     = htmltools::tagList(
                    title_icon, " Pipeline Log \u2014 ",
                    if (st == "success") "Completed" else "Failed"
                ),
                size      = "l",
                easyClose = TRUE,
                htmltools::pre(
                    class = "pipeline-log-content",
                    paste(log_content, collapse = "\n")
                ),
                footer = htmltools::tagList(
                    shiny::downloadButton(
                        ns("download_log"),
                        label = htmltools::tagList(
                            bsicons::bs_icon("download", size = "0.85em"), " Download"
                        ),
                        class = "btn btn-outline-secondary btn-sm"
                    ),
                    shiny::modalButton("Close")
                )
            ))
        })

        # ── Download log handler ────────────────────────────────────────────────
        output$download_log <- shiny::downloadHandler(
            filename = function() {
                p <- selected_log()
                if (!is.null(p) && nzchar(p)) basename(p) else "pipeline.log"
            },
            content = function(file) {
                p <- selected_log()
                if (!is.null(p) && nzchar(p) && file.exists(p)) {
                    file.copy(p, file)
                } else {
                    writeLines("Log file not available.", file)
                }
            }
        )

        invisible(NULL)
    })
}
