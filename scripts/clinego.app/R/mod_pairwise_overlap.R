#' Pairwise Trait Overlap UI
#'
#' On-demand async computation (background Rscript) of pairwise trait overlaps
#' at the thresholds and clumping distances currently set in the GEA/GWAS filter
#' bars. A "Compute" button triggers the run; a spinner shows while it runs;
#' results are cached to disk keyed by a parameter hash so they survive restart.
#'
#' @param id module namespace id
#' @noRd
mod_pairwise_overlap_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        bslib::card(
            full_screen = FALSE,
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span(
                    bsicons::bs_icon("intersect"),
                    " Pairwise Trait Overlap"
                ),
                shiny::downloadButton(ns("dl_pairwise"), "CSV",
                                      class = "btn-sm btn-outline-secondary")
            ),
            bslib::card_body(
                htmltools::div(
                    class = "alert alert-info py-2 mb-2 small",
                    bsicons::bs_icon("info-circle"), " ",
                    "Uses the thresholds and clumping distances from the GEA and GWAS filter bars above.",
                    " Click ", htmltools::strong("Compute"), " to run (or re-run after changing parameters).",
                    " Results are cached to disk and reload automatically on app restart."
                ),
                htmltools::div(
                    class = "control-bar mb-2",
                    bslib::layout_columns(
                        col_widths = c(4, 4, 4),
                        shiny::selectInput(
                            ns("ctype_filter"), "Comparison type",
                            choices  = c("All", "GEA-GEA", "GWAS-GWAS", "GEA-GWAS"),
                            selected = "All"
                        ),
                        htmltools::div(
                            class = "mt-4",
                            shiny::actionButton(
                                ns("run_pairwise"),
                                label = htmltools::tagList(
                                    bsicons::bs_icon("play-fill"), " Compute pairwise overlaps"
                                ),
                                class = "btn-primary"
                            )
                        ),
                        shiny::uiOutput(ns("status_badge"))
                    )
                ),
                shiny::uiOutput(ns("pairwise_ui"))
            )
        ),

        # Per-pair Miami — shown when a row is selected
        shiny::uiOutput(ns("miami_panel_ui"))
    )
}

#' Pairwise Trait Overlap Server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param coords reactive coords list (chr_offsets + chr_lengths)
#' @param gea_sigsnps reactive — named list of per-method sig-SNP data.tables (GEA side)
#' @param gwas_sigsnps reactive — named list of per-method sig-SNP data.tables (GWAS side)
#' @param pairwise_params reactive — list with threshold/regime/window params + fingerprint
#' @noRd
mod_pairwise_overlap_server <- function(id, project_data, coords,
                                         gea_sigsnps, gwas_sigsnps,
                                         pairwise_params) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Parameter hash (cache key) ─────────────────────────────────────────
        pairwise_hash <- shiny::reactive({
            p <- pairwise_params()
            if (is.null(p)) return(NULL)
            substr(digest::digest(p, algo = "md5"), 1, 8)
        })

        # ── Async state ────────────────────────────────────────────────────────
        subprocess_handles <- shiny::reactiveValues()  # hash -> handle
        pairwise_cache     <- shiny::reactiveValues()  # hash -> result list

        is_running <- function(hash) {
            if (is.null(hash)) return(FALSE)
            h <- subprocess_handles[[hash]]
            !is.null(h) && h$process$is_alive()
        }

        # ── Disk hydration: auto-load cached results on hash change ────────────
        shiny::observe({
            hash <- pairwise_hash()
            if (is.null(hash)) return()
            if (!is.null(pairwise_cache[[hash]]) || is_running(hash)) return()
            pd <- project_data()
            tbl_path <- file.path(ondemand_pairwise_dir(pd$name, hash),
                                   "pairwise_overlap_table.tsv")
            if (file_ok(tbl_path)) {
                handle_stub <- list(
                    out_dir    = ondemand_pairwise_dir(pd$name, hash),
                    gea_window  = pairwise_params()$gea_window  %||% 1000000L,
                    gwas_window = pairwise_params()$gwas_window %||% 1000000L
                )
                result <- tryCatch(.collect_pairwise_result(handle_stub),
                                   error = function(e) list(error = conditionMessage(e)))
                pairwise_cache[[hash]] <- result
            }
        })

        # ── Poller ─────────────────────────────────────────────────────────────
        shiny::observe({
            handles_list <- shiny::reactiveValuesToList(subprocess_handles)
            active <- Filter(Negate(is.null), handles_list)
            if (length(active) == 0) return()
            shiny::invalidateLater(1000)
            for (hash in names(active)) {
                h <- active[[hash]]
                if (h$process$is_alive()) next
                status <- tryCatch(h$process$get_exit_status(), error = function(e) 1L)
                result <- tryCatch(.collect_pairwise_result(h),
                                   error = function(e) list(error = conditionMessage(e)))
                if (!is.null(result$error) && status != 0L) {
                    log_lines <- tryCatch(readLines(h$log_file, warn = FALSE),
                                          error = function(e) character())
                    result$log_tail <- paste(utils::tail(log_lines, 20), collapse = "\n")
                }
                pairwise_cache[[hash]] <- result
                subprocess_handles[[hash]] <- NULL
            }
        })

        # ── Run button ─────────────────────────────────────────────────────────
        shiny::observeEvent(input$run_pairwise, {
            hash <- pairwise_hash()
            if (is.null(hash) || is_running(hash)) return()
            # Allow re-run even if cached (user may want fresh run after pipeline re-run)
            p  <- pairwise_params()
            pd <- project_data()
            handle <- tryCatch(
                launch_pairwise_subprocess(
                    gea_sig_list  = gea_sigsnps(),
                    gwas_sig_list = gwas_sigsnps(),
                    gea_window    = p$gea_window  %||% 1000000L,
                    gwas_window   = p$gwas_window %||% 1000000L,
                    min_snps      = p$min_snps    %||% 2L,
                    hash          = hash,
                    project_data  = pd
                ),
                error = function(e)
                    list(error = paste("Failed to launch pairwise computation:", conditionMessage(e)))
            )
            if (!is.null(handle$error)) {
                pairwise_cache[[hash]] <- handle  # store error immediately
            } else {
                # Clear old cached result for this hash so the spinner shows
                pairwise_cache[[hash]] <- NULL
                subprocess_handles[[hash]] <- handle
            }
        })

        # ── Status badge (stale indicator) ────────────────────────────────────
        output$status_badge <- shiny::renderUI({
            hash <- pairwise_hash()
            if (is.null(hash)) return(NULL)
            cached <- pairwise_cache[[hash]]
            running <- is_running(hash)
            if (running) {
                htmltools::span(class = "badge bg-warning mt-4 ms-1",
                                bsicons::bs_icon("arrow-repeat", class = "spin-icon"),
                                " Computing...")
            } else if (!is.null(cached) && is.null(cached$error)) {
                n <- if (!is.null(cached$table)) nrow(cached$table) else 0L
                htmltools::span(class = "badge bg-success mt-4 ms-1",
                                bsicons::bs_icon("check-circle"), sprintf(" %d pairs", n))
            } else if (!is.null(cached) && !is.null(cached$error)) {
                htmltools::span(class = "badge bg-danger mt-4 ms-1",
                                bsicons::bs_icon("x-circle"), " Failed")
            } else {
                htmltools::span(class = "badge bg-secondary mt-4 ms-1",
                                bsicons::bs_icon("dash-circle"), " Not computed")
            }
        })

        # ── Main table UI ──────────────────────────────────────────────────────
        output$pairwise_ui <- shiny::renderUI({
            hash <- pairwise_hash()
            running <- is_running(hash)
            cached  <- if (!is.null(hash)) pairwise_cache[[hash]] else NULL

            if (running) {
                return(htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-3",
                    bsicons::bs_icon("arrow-repeat", class = "text-warning spin-icon flex-shrink-0"),
                    htmltools::span(class = "small",
                        "Computing pairwise trait overlaps in the background...")
                ))
            }

            if (!is.null(cached) && !is.null(cached$error)) {
                return(htmltools::div(
                    htmltools::div(
                        class = "alert alert-danger small",
                        bsicons::bs_icon("x-circle"), " ",
                        "Computation failed: ", cached$error,
                        if (!is.null(cached$log_tail))
                            htmltools::pre(class = "pipeline-log-content mt-2 small",
                                           cached$log_tail)
                    ),
                    shiny::actionButton(ns("run_pairwise"), "Retry",
                                         class = "btn-sm btn-warning")
                ))
            }

            if (is.null(cached)) {
                return(htmltools::div(
                    class = "text-muted small fst-italic py-3 text-center",
                    "Click “Compute pairwise overlaps” to compute for the current thresholds and clumping distances."
                ))
            }

            # Result available
            DT::DTOutput(ns("pairwise_table"))
        })

        # ── Result data reactives ──────────────────────────────────────────────
        current_result <- shiny::reactive({
            hash <- pairwise_hash()
            if (is.null(hash)) return(NULL)
            pairwise_cache[[hash]]
        })

        pairwise_tbl <- shiny::reactive({
            r <- current_result()
            if (is.null(r) || !is.null(r$error)) return(data.table::data.table())
            r$table %||% data.table::data.table()
        })

        collapsed_snps <- shiny::reactive({
            r <- current_result()
            if (is.null(r) || !is.null(r$error)) return(data.table::data.table())
            r$collapsed %||% data.table::data.table()
        })

        pairwise_window <- shiny::reactive({
            r <- current_result()
            if (is.null(r)) return(list(gea = 500000L, gwas = 500000L))
            r$windows %||% list(gea = 500000L, gwas = 500000L)
        })

        # ── Filtered table for display ─────────────────────────────────────────
        filtered_tbl <- shiny::reactive({
            tbl   <- pairwise_tbl()
            if (nrow(tbl) == 0) return(tbl)
            ctype <- input$ctype_filter
            if (!is.null(ctype) && ctype != "All")
                tbl <- tbl[tbl$comparison_type == ctype, ]
            tbl
        })

        output$pairwise_table <- DT::renderDataTable({
            tbl <- filtered_tbl()
            if (nrow(tbl) == 0) {
                return(safe_datatable(data.frame(
                    message = "No overlapping traits found at the current thresholds and clumping distances."
                )))
            }
            display <- as.data.frame(tbl[, c("trait_a","source_a","n_snps_a",
                                              "trait_b","source_b","n_snps_b",
                                              "exact_matches","window_overlaps",
                                              "overlap_score","comparison_type")])
            names(display) <- c("Trait A","Source A","SNPs A",
                                 "Trait B","Source B","SNPs B",
                                 "Exact","Window","Score","Type")
            safe_datatable(
                display,
                selection = "single",
                options   = list(
                    dom        = "frtip",
                    pageLength = 15,
                    scrollX    = TRUE,
                    order      = list(list(8L, "desc"))
                ),
                rownames = FALSE
            )
        })

        output$dl_pairwise <- shiny::downloadHandler(
            filename = function() paste0("pairwise_overlap_", project_data()$name, ".csv"),
            content  = function(file) {
                utils::write.csv(as.data.frame(pairwise_tbl()), file, row.names = FALSE)
            }
        )

        # ── Selected pair ──────────────────────────────────────────────────────
        selected_pair <- shiny::reactive({
            rows <- input$pairwise_table_rows_selected
            if (length(rows) == 0) return(NULL)
            tbl <- filtered_tbl()
            if (nrow(tbl) < rows) return(NULL)
            as.list(tbl[rows, ])
        })

        # ── Per-pair Miami panel UI ────────────────────────────────────────────
        output$miami_panel_ui <- shiny::renderUI({
            pair <- selected_pair()
            if (is.null(pair)) return(NULL)
            bslib::card(
                bslib::card_header(
                    htmltools::div(
                        class = "d-flex align-items-center gap-2",
                        htmltools::span(
                            bsicons::bs_icon("bar-chart-steps"),
                            sprintf(" %s  vs  %s", pair$trait_a, pair$trait_b)
                        ),
                        help_note("pairwise_miami")
                    )
                ),
                bslib::card_body(
                    plotly::plotlyOutput(ns("pairwise_miami"), height = "420px")
                )
            )
        })

        # ── Per-pair Miami plot ────────────────────────────────────────────────
        output$pairwise_miami <- plotly::renderPlotly({
            pair       <- selected_pair()
            coords_val <- coords()
            if (is.null(pair) || is.null(coords_val)) return(plotly::plot_ly())

            snps <- collapsed_snps()
            if (nrow(snps) == 0) return(plotly::plot_ly())

            dt_a <- snps[snps$trait == pair$trait_a, ]
            dt_b <- snps[snps$trait == pair$trait_b, ]
            if (nrow(dt_a) == 0 || nrow(dt_b) == 0) return(plotly::plot_ly())

            dt_a <- data.table::copy(dt_a)
            dt_b <- data.table::copy(dt_b)
            dt_a[, pvalue := min_pvalue]
            dt_b[, pvalue := min_pvalue]
            dt_a <- add_cum_pos(dt_a, coords_val)
            dt_b <- add_cum_pos(dt_b, coords_val)

            # Window for highlight: max of the two filter-bar clumping distances
            wins <- pairwise_window()
            highlight_window <- max(wins$gea %||% 500000L, wins$gwas %||% 500000L)

            exact_ids   <- intersect(dt_a$SNPID, dt_b$SNPID)
            window_snps <- .pairwise_window_snps(dt_a, dt_b, highlight_window, 1L)
            overlap_a   <- union(exact_ids, window_snps$a_ids)
            overlap_b   <- union(exact_ids, window_snps$b_ids)

            dt_a[, overlapping := SNPID %in% overlap_a]
            dt_b[, overlapping := SNPID %in% overlap_b]
            dt_a[, exact       := SNPID %in% exact_ids]
            dt_b[, exact       := SNPID %in% exact_ids]

            .build_pairwise_plotly(dt_a, dt_b, pair, exact_ids, window_snps, coords_val)
        })
    })
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Single-linkage clustering (returns cluster vector for sorted positions)
#' @noRd
.cluster_snps_pw <- function(pos, window_size) {
    n <- length(pos)
    if (n == 0L) return(integer(0))
    cluster <- integer(n)
    k <- 1L
    cluster[1L] <- k
    if (n > 1L) {
        for (i in 2L:n) {
            if (is.na(pos[i]) || is.na(pos[i - 1L]) || (pos[i] - pos[i - 1L]) > window_size) {
                k <- k + 1L
            }
            cluster[i] <- k
        }
    }
    cluster
}

#' Find SNP IDs that participate in window overlaps between two trait SNP sets
#' Returns list(a_ids, b_ids, windows) where windows is data.table(chr, x0, x1)
#' @noRd
.pairwise_window_snps <- function(dt_a, dt_b, window_size, min_snps) {
    all_pos <- data.table::rbindlist(list(
        dt_a[, .(chr, pos, SNPID, grp = "A", cum_pos)],
        dt_b[, .(chr, pos, SNPID, grp = "B", cum_pos)]
    ))
    data.table::setorder(all_pos, chr, pos)
    a_ids <- character(0)
    b_ids <- character(0)
    win_rows <- list()
    for (ch in unique(all_pos$chr)) {
        sub  <- all_pos[chr == ch]
        clus <- .cluster_snps_pw(sub$pos, window_size)
        sub[, cluster := clus]
        cl_summ <- sub[, .(
            n_a = sum(grp == "A"), n_b = sum(grp == "B"),
            start = min(cum_pos), end = max(cum_pos),
            start_pos = min(pos),  end_pos = max(pos)
        ), by = cluster]
        keep <- cl_summ[n_a >= min_snps & n_b >= min_snps]
        if (nrow(keep) > 0) {
            for (cl_id in keep$cluster) {
                cl_snps <- sub[cluster == cl_id]
                a_ids <- union(a_ids, cl_snps[grp == "A", SNPID])
                b_ids <- union(b_ids, cl_snps[grp == "B", SNPID])
                row   <- keep[cluster == cl_id]
                win_rows[[length(win_rows) + 1L]] <- data.table::data.table(
                    chr = ch, x0 = row$start, x1 = row$end
                )
            }
        }
    }
    list(
        a_ids   = a_ids,
        b_ids   = b_ids,
        windows = if (length(win_rows) > 0) data.table::rbindlist(win_rows) else
                  data.table::data.table(chr = character(), x0 = numeric(), x1 = numeric())
    )
}

#' Build the pairwise Miami plotly
#' @noRd
.build_pairwise_plotly <- function(dt_a, dt_b, pair, exact_ids, window_snps, coords) {
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")
    col_a <- okabe_ito[1L]
    col_b <- okabe_ito[2L]

    y_max <- max(c(dt_a$log10p, dt_b$log10p), na.rm = TRUE) * 1.05

    shapes <- list()

    wins <- window_snps$windows
    if (nrow(wins) > 0) {
        for (i in seq_len(nrow(wins))) {
            shapes[[length(shapes) + 1L]] <- list(
                type    = "rect", xref = "x", yref = "paper",
                x0      = wins$x0[i], x1 = wins$x1[i],
                y0      = 0, y1 = 1,
                fillcolor = "rgba(128, 0, 128, 0.08)",
                line    = list(color = "rgba(128, 0, 128, 0.25)", width = 1)
            )
        }
    }

    if (length(exact_ids) > 0) {
        shared_a <- dt_a[SNPID %in% exact_ids, .(cum_pos_a = cum_pos[1L], log10p_a = log10p[1L]), by = SNPID]
        shared_b <- dt_b[SNPID %in% exact_ids, .(log10p_b = log10p[1L]),                           by = SNPID]
        shared_m <- shared_a[shared_b, on = "SNPID", nomatch = 0L]
        for (i in seq_len(nrow(shared_m))) {
            shapes[[length(shapes) + 1L]] <- list(
                type = "line", xref = "x", yref = "y",
                x0   = shared_m$cum_pos_a[i], x1 = shared_m$cum_pos_a[i],
                y0   = shared_m$log10p_a[i],  y1 = -shared_m$log10p_b[i],
                line = list(color = "rgba(214, 158, 46, 0.5)", width = 1.5, dash = "dot")
            )
        }
    }

    chr_offsets <- unlist(coords$chr_offsets)
    for (off in chr_offsets[-1L]) {
        shapes[[length(shapes) + 1L]] <- list(
            type = "line", xref = "x", yref = "paper",
            x0 = off, x1 = off, y0 = 0, y1 = 1,
            line = list(color = "rgba(150,150,150,0.2)", width = 0.5)
        )
    }

    midpoints <- chr_midpoints(coords)

    opacity_a <- ifelse(dt_a$overlapping, 0.90, 0.20)
    size_a    <- ifelse(dt_a$exact, 10, 7)
    border_a  <- ifelse(dt_a$exact, "rgba(214,158,46,1)", col_a)

    opacity_b <- ifelse(dt_b$overlapping, 0.90, 0.20)
    size_b    <- ifelse(dt_b$exact, 10, 7)
    border_b  <- ifelse(dt_b$exact, "rgba(214,158,46,1)", col_b)

    p <- plotly::plot_ly(source = "pairwise_miami") |>
        plotly::add_trace(
            data      = dt_a, x = ~cum_pos, y = ~log10p,
            type      = "scatter", mode = "markers",
            name      = pair$trait_a,
            text      = ~paste0(SNPID, "<br>-log10p: ", round(log10p, 2)),
            hoverinfo = "text",
            marker    = list(color = col_a, opacity = opacity_a, size = size_a,
                             line = list(color = border_a, width = 1.5))
        ) |>
        plotly::add_trace(
            data      = dt_b, x = ~cum_pos, y = ~(-log10p),
            type      = "scatter", mode = "markers",
            name      = pair$trait_b,
            text      = ~paste0(SNPID, "<br>-log10p: ", round(log10p, 2)),
            hoverinfo = "text",
            marker    = list(color = col_b, opacity = opacity_b, size = size_b,
                             line = list(color = border_b, width = 1.5))
        ) |>
        plotly::add_segments(
            x    = 0, xend = max(c(dt_a$cum_pos, dt_b$cum_pos), na.rm = TRUE),
            y    = 0, yend = 0,
            line = list(color = "rgba(100,100,100,0.4)", width = 1),
            showlegend = FALSE, hoverinfo = "none"
        ) |>
        plotly::layout(
            shapes  = shapes,
            xaxis   = list(title = "", tickvals = midpoints, ticktext = names(midpoints),
                           tickangle = 0, showgrid = FALSE, zeroline = FALSE),
            yaxis   = list(title = "-log10(p)",
                           tickvals = pretty(c(-y_max, y_max)),
                           ticktext = as.character(abs(pretty(c(-y_max, y_max)))),
                           zeroline = FALSE),
            legend  = list(orientation = "h", x = 0, y = 1.05),
            margin  = list(l = 60, r = 20, t = 10, b = 40),
            hovermode = "closest"
        )

    plotly::config(p, displayModeBar = FALSE, responsive = TRUE)
}
