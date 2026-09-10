#' Manhattan overlay module UI
#'
#' Displays a plotly Manhattan with a static background PNG (non-sig SNPs)
#' and interactive sig-SNP markers overlaid.
#'
#' @param id module namespace id
#' @param height plotly output height (default "500px")
#' @noRd
mod_manhattan_overlay_ui <- function(id, height = "400px", filter_ui = NULL) {
    ns <- shiny::NS(id)
    bslib::card(
        full_screen = TRUE,
        bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            htmltools::div(
                class = "d-flex align-items-center gap-2",
                shiny::textOutput(ns("title"), inline = TRUE),
                shiny::uiOutput(ns("note"), inline = TRUE)
            ),
            htmltools::span(
                class = "d-flex gap-2 align-items-center",
                shiny::uiOutput(ns("show_regions_toggle")),
                bslib::popover(
                    trigger = bsicons::bs_icon("download", title = "Download plot"),
                    title = "Download",
                    shiny::downloadButton(ns("dl_svg"), "SVG",
                                          class = "btn-sm btn-outline-secondary"),
                    shiny::downloadButton(ns("dl_png"), "PNG",
                                          class = "btn-sm btn-outline-secondary")
                )
            )
        ),
        # Optional filter controls injected between header and plot
        if (!is.null(filter_ui)) filter_ui,
        bslib::card_body(
            class = "p-1",
            htmltools::div(
                class = "plot-loading-wrap",
                htmltools::div(
                    class = "plot-loading-overlay",
                    bsicons::bs_icon("arrow-repeat", class = "spin-icon"),
                    htmltools::span("Loading plot…")
                ),
                plotly::plotlyOutput(ns("overlay"), height = height)
            )
        )
    )
}

#' Manhattan overlay module server
#'
#' @param id module namespace id
#' @param project_data reactive returning project data bundle from make_project_data()
#' @param module character: MOD_GEA, MOD_GWAS, or MOD_GEAXGWAS
#' @param method character reactive: method name (or NULL for combined)
#' @param trait character reactive: trait name (or NULL for combined)
#' @param combined logical: TRUE = combined multi-trait/method view
#' @param is_miami logical: TRUE = Miami plot mode
#' @param title_label reactive character for card title
#' @param note reactive returning an htmltools tag (e.g. help_note()) rendered
#'   next to the title, or NULL for no badge
#' @param regions reactive data.table of regions (optional, for overlay shapes)
#' @param show_regions_control logical: show toggle for region highlights
#' @param sig_snps_override reactive returning a long-format data.table to use
#'   instead of loading from files. When non-NULL, the file-based loading is
#'   skipped entirely. Must have columns SNPID, chr, pos, pvalue, method, trait.
#'   Pass \code{shiny::reactive(NULL)} (default) to use the standard loaders.
#' @return reactive of clicked region_id (character or NULL)
#' @noRd
mod_manhattan_overlay_server <- function(id, project_data,
                                          module             = MOD_GEA,
                                          method             = shiny::reactive(NULL),
                                          trait              = shiny::reactive(NULL),
                                          combined           = FALSE,
                                          is_miami           = FALSE,
                                          title_label        = shiny::reactive("Manhattan"),
                                          note               = shiny::reactive(NULL),
                                          regions            = shiny::reactive(NULL),
                                          current_region_id  = shiny::reactive(NULL),
                                          show_regions_control = TRUE,
                                          sig_snps_override    = shiny::reactive(NULL),
                                          trait_colors         = shiny::reactive(NULL),
                                          method_shapes        = shiny::reactive(NULL),
                                          # Optional path overrides (WZA regime, etc.)
                                          bg_path_override     = shiny::reactive(NULL),
                                          coords_path_override = shiny::reactive(NULL),
                                          # Live threshold line: -log10(threshold) scalar.
                                          # When non-NULL, drawn instead of coords$bonferroni_y.
                                          threshold_y          = shiny::reactive(NULL),
                                          # Overlapping-tab Miami extensions
                                          gwas_regions      = shiny::reactive(NULL),
                                          overlap_regions   = shiny::reactive(NULL)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Latching loading overlay ───────────────────────────────────────────
        # Bakes a debounced `plotly_afterplot` listener into every plotly widget
        # returned by output$overlay.  Once the plot settles after cold-start,
        # the `.plot-revealed` class is added to the wrapper and the CSS overlay
        # fades out — it is never removed by JS so param toggles stay smooth.
        # An 8 s max-timeout fallback ensures the overlay cannot stick forever.
        JS_REVEAL <- htmlwidgets::JS("function(el, x) {
            var wrap = el.closest('.plot-loading-wrap');
            if (!wrap) return;
            var debounce = null;
            if (!wrap.dataset.revealFallback) {
                wrap.dataset.revealFallback = '1';
                setTimeout(function() { wrap.classList.add('plot-revealed'); }, 8000);
            }
            function doReveal() {
                // Wait for SVG <image> elements (background PNG) to finish browser
                // decode before fading overlay — plotly_afterplot fires after Plotly's
                // JS cycle but before the large base64 PNG is decoded and painted.
                var imgs = Array.from(el.querySelectorAll('image'));
                if (!imgs.length) { wrap.classList.add('plot-revealed'); return; }
                var n = imgs.length;
                function done() { if (--n === 0) wrap.classList.add('plot-revealed'); }
                imgs.forEach(function(si) {
                    var src = si.getAttribute('href') || si.getAttribute('xlink:href') || '';
                    if (!src) { done(); return; }
                    var tmp = new Image();
                    tmp.onload = tmp.onerror = done;
                    tmp.src = src;
                });
            }
            function revealDebounced() {
                // INSTRUMENTATION — remove after blink fix is verified.
                console.log('[overlay] afterplot', performance.now());
                clearTimeout(debounce);
                debounce = setTimeout(doReveal, 350);
            }
            el.on('plotly_afterplot', revealDebounced);
            // No initial revealDebounced() — premature 350ms timer caused early reveal
            // on slow cold-start renders before plotly_afterplot had a chance to fire.
        }")
        reveal_when_settled <- function(p) htmlwidgets::onRender(p, JS_REVEAL)

        output$title <- shiny::renderText(title_label())
        output$note  <- shiny::renderUI(note())

        # Toggle for showing region rectangles (optional)
        output$show_regions_toggle <- shiny::renderUI({
            if (!show_regions_control) return(NULL)
            bslib::input_switch(ns("show_regions"), "Show regions", value = TRUE)
        })

        show_regions <- shiny::reactive({
            if (!show_regions_control) return(TRUE)
            isTRUE(input$show_regions)
        })

        # ── Resolve file paths ─────────────────────────────────────────────────
        bg_path <- shiny::reactive({
            # Caller-supplied override (e.g. WZA regime) takes priority
            ov <- bg_path_override()
            if (!is.null(ov)) return(ov)

            pd  <- project_data()
            k   <- pd$k_best
            if (is.na(k)) return(NULL)

            if (is_miami) {
                miami_bg_path(pd$name, k)
            } else if (combined) {
                combined_manhattan_bg_path(pd$name, module, k)
            } else {
                m <- method(); t <- trait()
                shiny::req(m, t)
                cfg <- pd$config
                adj <- resolve_adjust(cfg, m, if (module == MOD_GWAS) MOD_GWAS else MOD_GEA)
                if (is.null(adj)) return(NULL)
                manhattan_bg_path(pd$name, module, m, t, k, adj)
            }
        })

        coords_path <- shiny::reactive({
            ov <- coords_path_override()
            if (!is.null(ov)) return(ov)

            pd <- project_data()
            k  <- pd$k_best
            if (is.na(k)) return(NULL)

            if (is_miami)       miami_coords_path(pd$name, k)
            else if (combined)  combined_manhattan_coords_path(pd$name, module, k)
            else {
                m <- method(); t <- trait()
                shiny::req(m, t)
                adj <- resolve_adjust(pd$config, m,
                    if (module == MOD_GWAS) MOD_GWAS else MOD_GEA)
                if (is.null(adj)) return(NULL)
                manhattan_coords_path(pd$name, module, m, t, k, adj)
            }
        })

        # ── Load data ──────────────────────────────────────────────────────────
        coords <- shiny::reactive({
            load_coords(coords_path())
        })

        sig_snps <- shiny::reactive({
            pd <- project_data()

            # If a caller-supplied reactive override is provided, use it directly.
            # This enables the interactive filter/strategy switcher in the association tab.
            override_val <- sig_snps_override()
            if (!is.null(override_val)) return(override_val)

            if (!combined && !is_miami) {
                # Per-method view: load method-specific file with actual per-method p-values.
                # The combined selected_snps.tsv uses min_pvalue across all methods, which
                # causes sig SNP dots to be clipped above the method-specific y_range.
                m   <- method(); t <- trait()
                k   <- pd$k_best
                adj <- resolve_adjust(pd$config, m,
                    if (module == MOD_GWAS) MOD_GWAS else MOD_GEA)
                shiny::req(m, k, adj)
                cache_key <- paste0("method_sigsnps_", pd$name, "_", module,
                                    "_", m, "_", k, "_", adj)
                # Include file mtime so cache is bypassed when the pipeline regenerates
                # the sig-SNP file (e.g., after adding/removing bioclimatic factors).
                sig_path <- method_sigsnps_direct_path(pd$name, module, m, k, adj)
                fp <- if (file.exists(sig_path)) as.character(file.info(sig_path)$mtime) else "missing"
                load_cached(cache_key, function() {
                    load_method_sigsnps(pd$name, module, m, k, adj)
                }, fingerprint = fp)
            } else {
                # Combined / Miami: use selected_snps.tsv (min_pvalue is correct here
                # because the combined background PNG spans all methods' y-range)
                path <- selected_snps_path(pd$name, module)
                fp   <- if (file.exists(path)) as.character(file.info(path)$mtime) else "missing"
                load_cached(paste0("sig_snps_", path), function() {
                    load_selected_snps(pd$name, module, k_best = pd$k_best)
                }, fingerprint = fp)
            }
        })

        # ── Encode background ─────────────────────────────────────────────────
        # Hoisted outside renderPlotly so the PNG is read and base64-encoded only
        # when bg_path() changes (project / K / method switch), not on every
        # reactive invalidation (e.g., region click, threshold change).
        bg_uri <- shiny::reactive({
            b <- bg_path()
            if (file_ok(b)) encode_background_png(b) else NULL
        })

        # ── Render plotly ──────────────────────────────────────────────────────
        # Informative blank plot used for the "no coords yet" case and as a fail-soft
        # fallback when rendering errors — a single bad output must never crash the app.
        overlay_placeholder <- function(msg) {
            plotly::plot_ly(type = "scatter", mode = "markers") |>
                plotly::layout(
                    xaxis = list(visible = FALSE),
                    yaxis = list(visible = FALSE),
                    annotations = list(list(
                        text      = msg,
                        x         = 0.5, y = 0.5,
                        xref      = "paper", yref = "paper",
                        showarrow = FALSE,
                        font      = list(size = 14, color = "#888888")
                    )),
                    paper_bgcolor = "rgba(0,0,0,0)",
                    plot_bgcolor  = "rgba(0,0,0,0)"
                )
        }

        output$overlay <- plotly::renderPlotly({
          tryCatch({
            co <- coords()
            # coords.json absent: pipeline hasn't run yet or was interrupted mid-run.
            # Show an informative placeholder instead of a silent blank.
            if (is.null(co)) {
                return(reveal_when_settled(overlay_placeholder(
                    "Manhattan plot not available — run this mode from the sidebar to generate it.")))
            }
            bg_uri_val <- bg_uri()

            # INSTRUMENTATION — remove after blink fix is verified.
            # Counts how many times renderPlotly actually re-executes so we can
            # confirm the fix collapsed cold-load to a single render.
            message("[overlay render] id=", id, " bg_ok=", !is.null(bg_uri_val),
                    " snps_rows=", if (!is.null(sig_snps())) nrow(sig_snps()) else "NULL",
                    " thr_y=", shiny::isolate(threshold_y()))

            snps <- sig_snps()
            # For per-method views the loaded data is already method-specific;
            # only filter by trait. For combined/Miami, no filtering needed.
            if (!combined && !is_miami && nrow(snps) > 0) {
                t <- trait()
                if (!is.null(t) && "trait" %in% names(snps))
                    snps <- snps[snps$trait == t, ]
            }

            # Add cumulative positions for sig SNPs
            if (nrow(snps) > 0 && !"cum_pos" %in% names(snps)) {
                snps <- add_cum_pos(snps, co)
            }

            # Miami: GWAS traits need negative log10p (bottom half of plot)
            if (is_miami && nrow(snps) > 0 && !is.null(co$gwas_traits)) {
                gwas_tr <- unlist(co$gwas_traits)
                snps <- data.table::copy(snps)
                snps[trait %in% gwas_tr, log10p := -log10p]
            }

            # ── Shape inputs: isolate so changes don't trigger a full re-render ──
            #
            # threshold_y(), regions(), gwas_regions(), overlap_regions(), show_regions()
            # are shape-only: they don't affect the plotly widget skeleton (data traces,
            # axis ranges, background image) — only the shapes[] overlay.  Isolating them
            # here means late-resolving reactives (threshold_y is debounced 500 ms in
            # parent modules) and project-data-dependent region tables don't cause extra
            # cold-load renders.  The shape_trigger observer below picks up all
            # subsequent changes and pushes them via plotlyProxy relayout instead.
            gwas_reg  <- shiny::isolate(gwas_regions())
            ov_reg    <- shiny::isolate(overlap_regions())
            use_miami_shapes <- is_miami &&
                !is.null(gwas_reg) && !is.null(ov_reg)

            thr_y <- shiny::isolate(threshold_y())

            if (use_miami_shapes) {
                # 4-category Miami shapes (GEA-only gray / GWAS-only gray / overlap orange / selected blue)
                custom_shapes <- if (shiny::isolate(show_regions())) {
                    build_miami_region_shapes(
                        gea_regions                = shiny::isolate(regions()),
                        gwas_regions               = gwas_reg,
                        overlap_regions            = ov_reg,
                        selected_overlap_region_id = shiny::isolate(current_region_id()),
                        coords                     = co
                    )
                } else {
                    list()
                }
                reveal_when_settled(build_manhattan_plotly(
                    bg_uri            = bg_uri_val,
                    coords            = co,
                    sig_snps          = if (nrow(snps) > 0) snps else NULL,
                    regions           = NULL,
                    current_region_id = NULL,
                    trait_colors      = trait_colors(),
                    method_shapes     = method_shapes(),
                    is_miami          = is_miami,
                    source            = ns("overlay"),
                    extra_shapes      = custom_shapes,
                    threshold_y       = thr_y
                ))
            } else {
                reg_data <- if (shiny::isolate(show_regions())) shiny::isolate(regions()) else NULL
                reveal_when_settled(build_manhattan_plotly(
                    bg_uri            = bg_uri_val,
                    coords            = co,
                    sig_snps          = if (nrow(snps) > 0) snps else NULL,
                    regions           = reg_data,
                    current_region_id = shiny::isolate(current_region_id()),
                    trait_colors      = trait_colors(),
                    method_shapes     = method_shapes(),
                    is_miami          = is_miami,
                    source            = ns("overlay"),
                    threshold_y       = thr_y
                ))
            }
          }, error = function(e) {
              # Fail soft: one bad output degrades to a placeholder instead of
              # terminating the whole Shiny session.
              warning("Manhattan overlay render failed: ", conditionMessage(e))
              reveal_when_settled(overlay_placeholder("Could not render Manhattan overlay."))
          })
        })

        # ── Shape trigger → incremental shape relayout (no full re-render) ──
        #
        # All shape-only inputs (current_region_id, threshold_y, regions,
        # gwas_regions, overlap_regions, show_regions) are isolated inside
        # renderPlotly above so they don't cause extra cold-load re-renders.
        # This observer fires on ANY of their changes and pushes the updated
        # shapes array via plotlyProxy relayout, which is instant and blink-free.
        #
        # IMPORTANT: every relayout must include BOTH region rects AND threshold
        # lines in one combined array because Plotly's relayout replaces the entire
        # shapes list — sending only region shapes would erase the threshold line.
        #
        # ignoreInit = TRUE: first establishment is skipped because renderPlotly
        # already drew correct initial shapes from the isolated reads above.
        shape_trigger <- shiny::reactive({
            list(current_region_id(), threshold_y(), regions(),
                 gwas_regions(), overlap_regions(), show_regions())
        })
        shiny::observeEvent(shape_trigger(), {
            co     <- shiny::isolate(coords())
            if (is.null(co)) return()

            thr_y    <- shiny::isolate(threshold_y())
            reg_data <- shiny::isolate(regions())
            gwas_reg <- shiny::isolate(gwas_regions())
            ov_reg   <- shiny::isolate(overlap_regions())
            show_reg <- shiny::isolate(show_regions())
            cur_rid  <- current_region_id()

            use_miami_shapes <- is_miami &&
                !is.null(gwas_reg) && !is.null(ov_reg)

            region_shapes <- if (!show_reg) {
                list()
            } else if (use_miami_shapes) {
                build_miami_region_shapes(
                    gea_regions                = reg_data,
                    gwas_regions               = gwas_reg,
                    overlap_regions            = ov_reg,
                    selected_overlap_region_id = cur_rid,
                    coords                     = co
                )
            } else if (!is.null(reg_data) && nrow(reg_data) > 0) {
                build_dual_region_shapes(reg_data, cur_rid, co,
                    y_lo = co$y_range[1], y_hi = co$y_range[2])
            } else {
                list()
            }

            # Replicate threshold lines using shared helper (fct_manhattan.R)
            combined_shapes <- c(region_shapes,
                                 build_threshold_lines(is_miami, thr_y, co))

            plotly::plotlyProxy("overlay", session) |>
                plotly::plotlyProxyInvoke("relayout", list(shapes = combined_shapes))
        }, ignoreInit = TRUE, ignoreNULL = FALSE)

        # ── Re-arm overlay on project switch ─────────────────────────────────
        # Switching projects reloads data cold — the overlay must cover that
        # blink too.  Remove plot-revealed so the overlay re-shows; the next
        # plotly_afterplot burst re-reveals it.  ignoreInit = TRUE so we don't
        # re-arm on the very first load (overlay starts visible by default).
        shiny::observeEvent(project_data(), {
            session$sendCustomMessage("plot_rearm", list(id = ns("overlay")))
        }, ignoreInit = TRUE)

        # ── Sig SNP click → region_id reactive ────────────────────────────────
        selected_region <- shiny::reactive({
            click <- plotly::event_data("plotly_click", source = ns("overlay"))
            if (is.null(click) || is.null(click$customdata)) return(NULL)
            click$customdata[[1]]
        })

        # ── Download handlers (static files) ──────────────────────────────────
        output$dl_png <- shiny::downloadHandler(
            filename = function() {
                paste0("manhattan_", project_data()$name, "_K",
                       project_data()$k_best, ".png")
            },
            content = function(file) {
                p <- bg_path()
                if (file_ok(p)) file.copy(p, file) else writeLines("", file)
            }
        )

        output$dl_svg <- shiny::downloadHandler(
            filename = function() {
                paste0("manhattan_", project_data()$name, "_K",
                       project_data()$k_best, ".svg")
            },
            content = function(file) {
                p <- bg_path()
                svg_path <- sub("_background\\.png$", ".svg", p)
                src <- if (file_ok(svg_path)) svg_path else p
                if (file_ok(src)) file.copy(src, file) else writeLines("", file)
            }
        )

        # Return selected region for parent module
        selected_region
    })
}
