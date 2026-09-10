

#' PreStructure tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_prestructure_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        k_values <- shiny::reactive({
            find_k_values(project_data()$name)
        })

        output$k_selector <- shiny::renderUI({
            ks <- k_values()
            if (length(ks) == 0) return(shiny::p("No K values found.", class = "text-muted small"))
            pd <- project_data()
            selected <- if (!is.na(pd$k_best) && pd$k_best %in% ks) pd$k_best else ks[1]
            shiny::selectInput(ns("k"), "K value", choices = ks, selected = selected)
        })

        selected_k <- shiny::reactive({
            k <- input$k
            if (is.null(k)) {
                ks <- k_values()
                if (length(ks) == 0) return(NULL)
                pd <- project_data()
                if (!is.na(pd$k_best) && pd$k_best %in% ks) pd$k_best else ks[1]
            } else {
                as.integer(k)
            }
        })

        # ── K-independent images ───────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            kr <- find_k_range(pd$name, pd$config)
            k_best <- pd$k_best

            n_samples <- {
                dt <- if (!is.null(k_best) && !is.na(k_best)) load_clusters(pd$name, k_best)
                      else data.table::data.table()
                if (nrow(dt) == 0) NULL else nrow(dt)
            }

            mod_image_card_server("pca",
                path    = shiny::reactive(pca_path(pd$name)),
                title   = shiny::reactive("PCA"),
                dl_name = shiny::reactive("pca"),
                note    = shiny::reactive(help_note("pca",
                    results = if (!is.null(n_samples)) paste0("N = ", n_samples, " samples.") else NULL
                ))
            )
            mod_image_card_server("tracy_widom",
                path    = shiny::reactive(tracy_widom_path(pd$name)),
                title   = shiny::reactive("Tracy-Widom Test"),
                dl_name = shiny::reactive("tracy_widom"),
                note    = shiny::reactive(help_note("tracy_widom"))
            )
            mod_image_card_server("cross_entropy",
                path    = shiny::reactive(kr$path),
                title   = shiny::reactive("Cross-Entropy"),
                dl_name = shiny::reactive("cross_entropy"),
                note    = shiny::reactive({
                    # Changing sNMF.k_start/k_end leaves the previous run's PNG in
                    # the same directory. find_k_range() now picks the one matching
                    # the config, but say the others are there — otherwise the folder
                    # silently holds two curves for one .snmfProject. Built BEFORE the
                    # k_best guard: a project with sNMF.k_best unset is exactly the one
                    # where the user is about to read k_best off this curve.
                    stale_txt <- if (length(kr$stale) > 0) paste0(
                        " Superseded cross-entropy plot",
                        if (length(kr$stale) == 1) "" else "s", " also present (",
                        paste(kr$stale, collapse = ", "),
                        ") — from an earlier K range, not this one.") else ""
                    if (is.null(k_best) || is.na(k_best)) {
                        if (!nzchar(stale_txt)) return(NULL)
                        return(help_note("cross_entropy", results = trimws(stale_txt)))
                    }
                    range_txt <- if (!is.na(kr$k_start) && !is.na(kr$k_end))
                        paste0(" Tested K ", kr$k_start, "–", kr$k_end, ".") else ""
                    # config_k_best() is NA when sNMF.k_best is unset — resolve_k_best()
                    # then falls back to the highest K found on disk, which is NOT the
                    # lowest-cross-entropy K. Word the note accordingly so it never
                    # overclaims a "best" that wasn't actually computed from this plot.
                    results_txt <- if (!is.na(config_k_best(pd$config))) {
                        paste0("Selected K = ", k_best, " (sNMF.k_best).", range_txt)
                    } else {
                        paste0("sNMF.k_best is not set — showing the highest K found on disk (",
                              k_best, "), not the lowest cross-entropy value.", range_txt,
                              " Read the curve above and set sNMF.k_best to the K where it flattens.")
                    }
                    results_txt <- paste0(results_txt, stale_txt)
                    help_note("cross_entropy", results = results_txt, label = paste0("K=", k_best))
                })
            )
        })

        # ── K-dependent images ─────────────────────────────────────────────────
        shiny::observe({
            k  <- shiny::req(selected_k())
            pd <- project_data()

            pop_summary <- cluster_pop_summary(load_clusters(pd$name, k), k)
            pop_results <- if (!is.null(pop_summary)) {
                rng <- if (pop_summary$min_n == pop_summary$max_n) as.character(pop_summary$min_n)
                       else paste0(pop_summary$min_n, "–", pop_summary$max_n)
                paste0(pop_summary$n_pops, " populations at K=", k, " — ",
                      rng, " samples/pop (N=", pop_summary$n_samples, ").")
            } else NULL

            mod_image_card_server("structure_bar",
                path    = shiny::reactive(structure_k_path(pd$name, k)),
                title   = shiny::reactive("Structure"),
                dl_name = shiny::reactive(paste0("structure_K", k)),
                note    = shiny::reactive(help_note("structure_bar", results = pop_results,
                                                    label = paste0("K=", k)))
            )
            mod_image_card_server("pca_structure",
                path    = shiny::reactive(pca_structure_path(pd$name, k)),
                title   = shiny::reactive("PCA (Clustered)"),
                dl_name = shiny::reactive(paste0("pca_structure_K", k)),
                note    = shiny::reactive(help_note("pca_structure", results = pop_results,
                                                    label = paste0("K=", k)))
            )
            mod_image_card_server("pop_diff",
                path    = shiny::reactive(pop_diff_path(pd$name, k)),
                title   = shiny::reactive("Pop Differentiation"),
                dl_name = shiny::reactive(paste0("pop_diff_K", k)),
                note    = shiny::reactive(help_note("pop_diff", results = pop_results,
                                                    label = paste0("K=", k)))
            )
        })
    })
}

#' prestructure tab UI — dashboard layout.
#'
#' Promoted from scripts/layout_lab/ on 2026-08-29; the comments inside
#' carry the measurements behind each sizing decision.
#' @param id module namespace id
#' @noRd
mod_prestructure_ui <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "prestructure",

        # PreStructure has no badge output, hence alert_id = NULL
        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: the anchor ──────────────────────────────────────────────
            htmltools::div(
                class = "lab-hero-col",
                lab_hero(mod_image_card_ui(ns("pca"))),
                bslib::card(
                    class = "lab-dock-card lab-dock-qmatrix",
                    bslib::card_header(
                        bsicons::bs_icon("table"),
                        " Ancestry coefficients (Q-matrix)"
                    ),
                    bslib::card_body(DT::DTOutput(ns("clusters_table")))
                )
            ),

            # ── RIGHT: plenty of plots ────────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",

                # The K selector drives structure_bar / pca_structure / pop_diff,
                # all of which live in this column. PCA on the left is
                # K-independent, so the control belongs here, not page-wide --
                # and it rides the Ancestry divider rather than claiming a
                # `.control-bar` row, which cost 64px to hold a 128px select.
                lab_section_header("Ancestry", icon = "diagram-3",
                                   aside = shiny::uiOutput(ns("k_selector"))),
                # PCA (Clustered) then the barplot -- the same K-dependent story
                # told two ways -- STACKED, not side by side.
                #
                # They cannot share a row at any useful size. pca_structure is
                # 1:1 and structure_bar is 4:1, so filling one row of height H
                # needs H + 4H of width: in this 904px column that pins H at
                # ~180px. Side by side at 1+3 columns they measured 202px and
                # 166px inside 466px tiles -- both width-bound, each stranded in
                # ~270px of vertical whitespace. Stacked, each tile is sized to
                # its own aspect and the figures roughly double.
                htmltools::div(
                    class = "lab-prestr-ancestry",
                    htmltools::div(class = "lab-prestr-pca",
                                   mod_image_card_ui(ns("pca_structure"))),
                    htmltools::div(class = "lab-prestr-bar",
                                   mod_image_card_ui(ns("structure_bar")))
                ),

                # Diagnostics are read once, when choosing K -- not on every
                # visit. Collapsed, they cost one 35px header instead of a full
                # tile row, and that reclaimed height goes to the Ancestry pair
                # above, which is what this column is FOR. (A closed panel is
                # display:none, so the three PNGs are not transmitted until it
                # is opened -- see lab_depth_below_fold()'s note.)
                bslib::accordion(
                    class = "lab-below-fold lab-prestr-diag",
                    open = FALSE,
                    bslib::accordion_panel(
                        value = "diag",
                        title = htmltools::tagList(bsicons::bs_icon("graph-up"),
                                                   " Diagnostics (Tracy-Widom, cross-entropy, population differentiation)"),
                        # tracy_widom and cross_entropy are 1:1 -> one column each.
                        # pop_diff is 1.33 -> two columns, filling the row exactly.
                        lab_thumb_grid(
                            cols = 4,
                            lab_thumb(mod_image_card_ui(ns("tracy_widom"))),
                            lab_thumb(mod_image_card_ui(ns("cross_entropy"))),
                            lab_thumb(mod_image_card_ui(ns("pop_diff")), span = 2)
                        )
                    )
                )
            )
        )
    )
}

