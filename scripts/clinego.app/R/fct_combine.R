# On-the-fly combine strategies for interactive Shiny exploration.
# Ports combine logic from scripts/combine_selected_snps.R using data.table
# (avoids GenomicRanges dependency in Shiny).

#' Okabe-Ito colour palette (8 colours, colour-blind safe)
#' @noRd
OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Map traits to Okabe-Ito colours (sorted order, consistent with build_manhattan_plotly)
#'
#' @param traits character vector of trait names
#' @return named character vector: trait → hex colour
#' @noRd
trait_color_map <- function(traits) {
    traits_sorted <- sort(unique(traits))
    n <- length(traits_sorted)
    colors <- rep(OKABE_ITO, ceiling(n / length(OKABE_ITO)))[seq_len(n)]
    stats::setNames(colors, traits_sorted)
}

#' Plotly symbol names for methods (8 symbols, reused via recycling)
#' @noRd
METHOD_SHAPES <- c("circle", "triangle-up", "square", "diamond",
                   "cross", "x", "star", "triangle-down")

#' Map methods to plotly marker symbols (sorted order, stable across filtering)
#'
#' @param methods character vector of method names
#' @return named character vector: method → plotly symbol name
#' @noRd
method_shape_map <- function(methods) {
    methods_sorted <- sort(unique(methods))
    n <- length(methods_sorted)
    shapes <- rep(METHOD_SHAPES, ceiling(n / length(METHOD_SHAPES)))[seq_len(n)]
    stats::setNames(shapes, methods_sorted)
}

#' Normalize legacy and alias strategy names to current canonical names
#' @noRd
.normalize_strategy <- function(s) {
    switch(s,
        # Legacy pipeline names
        All          = "Union",
        Sum          = "Union",
        Overlap      = "Cross-method",
        MethodOverlap = "Cross-method per-trait",
        PairOverlap  = "Cross-method per-trait",
        s  # pass through "Union", "Cross-method", "Cross-method per-trait" unchanged
    )
}

#' Combine per-method sig SNPs using a specified strategy (on-the-fly in Shiny)
#'
#' @param sigsnps_list Named list of data.tables, one per method.
#'   Each must have columns: SNPID, chr, pos, pvalue, method, trait.
#' @param strategy One of "All", "Overlap", "MethodOverlap", or a method name.
#'   Legacy values "Sum" and "PairOverlap" are accepted and normalized.
#' @param gap Integer. Max bp distance for spatial overlap (Overlap/MethodOverlap).
#' @return Long-format data.table with SNPID, chr, pos, pvalue, method, trait.
#'   Empty data.table if nothing passes.
#' @noRd
combine_sigsnps <- function(sigsnps_list, strategy = "All", gap = 200000L) {
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    strategy <- .normalize_strategy(strategy)

    # Drop empty methods
    sigsnps_list <- sigsnps_list[vapply(sigsnps_list, nrow, integer(1)) > 0]
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    methods_vec <- names(sigsnps_list)
    gap <- as.integer(gap)

    result <- if (strategy %in% methods_vec) {
        # ── Single-method ───────────────────────────────────────────────────────
        data.table::copy(sigsnps_list[[strategy]])

    } else if (strategy == "Union") {
        # ── Union: all sig SNPs from all methods ────────────────────────────────
        data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)

    } else if (strategy == "Cross-method") {
        # ── Cross-method spatial overlap (trait-agnostic) ───────────────────────
        # A SNP passes if any partner exists in another method within gap bp.
        selected_ids <- character(0)
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                ids <- .spatial_overlap_ids(sigsnps_list[[n1]], sigsnps_list[[n2]], gap)
                selected_ids <- c(selected_ids, ids)
            }
        }
        selected_ids <- unique(selected_ids)
        if (length(selected_ids) == 0) return(.empty_sigsnps())
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_snps[SNPID %in% selected_ids]

    } else if (strategy == "Cross-method per-trait") {
        # ── Cross-method spatial overlap, same-trait only ───────────────────────
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_traits <- unique(all_snps$trait)
        selected_rows <- list()
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                for (bio in all_traits) {
                    m1_bio <- sigsnps_list[[n1]][trait == bio]
                    m2_bio <- sigsnps_list[[n2]][trait == bio]
                    if (nrow(m1_bio) == 0 || nrow(m2_bio) == 0) next
                    ids1 <- .spatial_overlap_ids(m1_bio, m2_bio, gap)
                    ids2 <- .spatial_overlap_ids(m2_bio, m1_bio, gap)
                    if (length(ids1) > 0)
                        selected_rows <- c(selected_rows, list(m1_bio[SNPID %in% ids1]))
                    if (length(ids2) > 0)
                        selected_rows <- c(selected_rows, list(m2_bio[SNPID %in% ids2]))
                }
            }
        }
        if (length(selected_rows) == 0) return(.empty_sigsnps())
        unique(data.table::rbindlist(selected_rows, use.names = TRUE, fill = TRUE))

    } else {
        stop(paste0("Unknown combine strategy: '", strategy,
                    "'. Valid values: Union, Cross-method, Cross-method per-trait, or a method name."))
    }

    if (is.null(result) || nrow(result) == 0) return(.empty_sigsnps())

    # Deduplicate — same SNP can appear in multiple methods/traits
    result <- unique(result, by = c("SNPID", "method", "trait"))

    # Add min_pvalue column for combined Manhattan y-placement
    # (background PNG y-axis spans the minimum p-value per SNP across all methods)
    result[, min_pvalue := min(pvalue, na.rm = TRUE), by = "SNPID"]

    result
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Find SNP IDs from dt1 that have a partner in dt2 within gap bp on same chr
#' @noRd
.spatial_overlap_ids <- function(dt1, dt2, gap) {
    if (nrow(dt1) == 0 || nrow(dt2) == 0) return(character(0))

    # Expand dt1 to windows [pos-gap, pos+gap]; dt2 stays as points [pos, pos]
    win <- data.table::data.table(
        chr   = as.character(dt1$chr),
        start = pmax(1L, as.integer(dt1$pos) - gap),
        end   = as.integer(dt1$pos) + gap,
        SNPID = dt1$SNPID
    )
    pts <- data.table::data.table(
        chr   = as.character(dt2$chr),
        start = as.integer(dt2$pos),
        end   = as.integer(dt2$pos),
        SNPID = dt2$SNPID
    )
    data.table::setkey(win, chr, start, end)
    data.table::setkey(pts, chr, start, end)
    ov <- data.table::foverlaps(pts, win, type = "within", nomatch = NULL)
    if (nrow(ov) == 0) return(character(0))
    # Return IDs from dt1 (window side)
    unique(ov$SNPID)
}

#' Empty sig SNPs skeleton
#' @noRd
.empty_sigsnps <- function() {
    data.table::data.table(
        SNPID  = character(),
        chr    = character(),
        pos    = integer(),
        pvalue = numeric(),
        method = character(),
        trait  = character()
    )
}

#' Empty sig SNPs skeleton with region_id column (for interactive combine output)
#' @noRd
.empty_sigsnps_assoc <- function() {
    data.table::data.table(
        SNPID     = character(),
        chr       = character(),
        pos       = integer(),
        pvalue    = numeric(),
        method    = character(),
        trait     = character(),
        region_id = character()
    )
}

# ── Shared UI / server helpers ────────────────────────────────────────────────

#' Extract default threshold type and value from project config.
#'
#' Reads the first entry in GEA.configs or GWAS.configs.
#' Returns list(type="bonf", value=0.05) as safe fallback.
#' @noRd
default_threshold <- function(config, module = MOD_GEA) {
    config_key <- if (module == MOD_GWAS) "GWAS" else "GEA"
    configs <- config_get(config, config_key, "configs", default = list())
    if (length(configs) == 0) return(list(type = "bonf", value = 0.05))

    first <- configs[[1]]
    type  <- if (is.list(first)) (first$adjust  %||% "bonf") else "bonf"
    value <- if (is.list(first)) suppressWarnings(as.numeric(first$threshold %||% 0.05))
             else 0.05
    if (is.na(value) || value <= 0) value <- 0.05
    list(type = type, value = value)
}

#' Sensible default threshold value for a given type.
#' @noRd
threshold_value_default_for_type <- function(type) {
    switch(type,
        bonf   = 0.05,
        qval   = 0.05,
        top    = 100,
        custom = 1e-5,
        0.05
    )
}

#' Is the threshold value valid for the given type?
#' @noRd
threshold_value_valid_for_type <- function(type, value) {
    if (is.null(value) || is.na(value)) return(FALSE)
    v <- suppressWarnings(as.numeric(value))
    if (is.na(v)) return(FALSE)
    switch(type,
        bonf   = v > 0 && v <= 1,
        qval   = v > 0 && v <= 1,
        top    = v >= 1,
        custom = v > 0,
        FALSE
    )
}

#' Build the always-present threshold bar (regime switch + threshold type/value + hint).
#'
#' Separate from the trait matrix so it is never gated on traits being present.
#' Per-method/per-cell rules no longer live here -- they're set directly in the
#' trait x method matrix (build_filter_bar_ui()'s cell/row/column popups) via
#' R/fct_threshold_rules.R. This bar only owns the MASTER rule.
#'
#' @param ns           Shiny namespace function
#' @param input_prefix Character prefix for input IDs (e.g. "" or "gea_")
#' @param regime_value Logical: current WZA regime state
#' @param threshold_type_value  Character: "bonf"/"qval"/"top"/"custom"
#' @param threshold_value_value Numeric: current threshold value
#' @param show_apply_to_config Logical: render the "Apply rules to config"
#'   button, which collapses the current master/registry/override rules into
#'   config_state$working's GEA.configs/GWAS.configs (one adjust/threshold
#'   per method -- the pipeline has no per-trait granularity). Omitted by
#'   GEAxGWAS's two bars, which don't drive a pipeline run.
#' @return tagList
#' @noRd
#' Significance-threshold modes offered by the threshold bar.
#'
#' Hoisted out of build_threshold_bar_ui() so the selectInput and the collapsed
#' summary badge (mod_gea.R::config_badges) name the same mode the same way --
#' they sit far apart in the UI and would otherwise drift.
#' @noRd
THRESHOLD_TYPE_CHOICES <- c(
    "Bonferroni"     = "bonf",
    "FDR (qval)"     = "qval",
    "Top N SNPs"     = "top",
    "Custom (raw p)" = "custom"
)

#' One-line rendering of the active significance rule, for a summary badge.
#'
#' The threshold VALUE means something different per mode -- an alpha, an FDR
#' target, a count of SNPs per trait, a raw p cutoff -- so the badge spells the
#' unit out rather than printing a bare number next to a mode name. Mirrors the
#' per-mode wording of output$threshold_hint in mod_gea.R.
#'
#' @param type one of THRESHOLD_TYPE_CHOICES' values
#' @param value the numeric threshold currently in force
#' @noRd
format_threshold_rule <- function(type, value) {
    type <- type %||% "bonf"
    if (is.null(value) || is.na(value)) return(names(which(THRESHOLD_TYPE_CHOICES == type))[1] %||% type)
    num <- function(v) format(v, scientific = (v > 0 && v < 1e-3), trim = TRUE)
    switch(type,
        bonf   = paste0("Bonferroni \u03b1 ", num(value)),
        qval   = paste0("FDR q \u2264 ", num(value)),
        top    = paste0("Top ", format(round(value), trim = TRUE), " SNPs/trait"),
        custom = paste0("raw p < ", num(value)),
        paste0(type, " ", num(value))
    )
}

build_threshold_bar_ui <- function(ns, input_prefix = "",
                                   regime_value          = FALSE,
                                   threshold_type_value  = "bonf",
                                   threshold_value_value = 0.05,
                                   regime_context        = "gea",
                                   show_apply_to_config  = FALSE) {
    pid <- function(name) ns(paste0(input_prefix, name))

    htmltools::div(
        class = "threshold-bar mb-2",
        # Regime switch row
        htmltools::div(
            class = "filter-row align-items-center mb-2",
            bslib::input_switch(
                pid("regime"), "WZA regime",
                value = isTRUE(regime_value)
            ),
            htmltools::span(
                class = "text-muted small ms-2",
                "Toggle to use Weighted-Z Analysis windows instead of per-SNP p-values"
            ),
            wza_window_note(regime_context)
        ),
        # Threshold type + value
        htmltools::div(
            class = "d-flex align-items-end gap-2",
            htmltools::div(
                class = "d-flex flex-column",
                htmltools::span("Significance threshold", class = "filter-label mb-1"),
                shiny::selectInput(
                    pid("threshold_type"),
                    label = NULL,
                    choices = THRESHOLD_TYPE_CHOICES,
                    selected = threshold_type_value,
                    width = "150px"
                )
            ),
            htmltools::div(
                class = "d-flex flex-column",
                htmltools::span(" ", class = "filter-label mb-1"),
                shiny::numericInput(
                    pid("threshold_value"),
                    label = NULL,
                    value = threshold_value_value,
                    min = 0, step = NA,
                    width = "110px"
                )
            ),
            shiny::uiOutput(pid("threshold_hint")),
            if (show_apply_to_config) htmltools::div(
                class = "d-flex flex-column ms-2",
                htmltools::span(" ", class = "filter-label mb-1"),
                shiny::actionButton(
                    pid("apply_to_config"),
                    label = htmltools::tagList(bsicons::bs_icon("cloud-arrow-up"),
                                               " Apply rules to config"),
                    class = "btn btn-outline-secondary btn-sm",
                    title = paste(
                        "Collapse the master/registry/override rules into one",
                        "adjust+threshold per method and write it into the",
                        "project config (Run or Save Project Files to persist to YAML)."
                    )
                )
            )
        )
    )
}

#' Build the trait x method matrix + strategy + clumping filter bar.
#'
#' Pure function -- call inside renderUI. Handles its own namespacing via `ns`.
#' The threshold/regime controls live in a separate always-present build_threshold_bar_ui().
#'
#' Per-cell significance rules are set directly here: click a cell's threshold
#' text to open a rule popup for that (trait, method); click the gear icon
#' revealed on hovering a row/column header to open a bulk popup for that
#' whole trait (all methods) or method (all traits). Row/column headers keep
#' their existing click-to-toggle-all behaviour -- the gear is a distinct
#' target so the two gestures never collide (see fct_threshold_rules.R for
#' the resulting override map and its precedence: cell override > registry
#' default (RDA) > master).
#'
#' @param ns                        Shiny namespace function from session$ns
#' @param traits                    Character vector of ALL trait names (full list, not just sig ones)
#' @param methods                   Character vector of method names
#' @param trait_colors              Named character vector: trait -> hex colour (from all_trait_names for stability)
#' @param combo_counts              Named list: "trait::method" -> integer SNP/window count
#' @param combo_thresholds          Named list: "trait::method" -> threshold p-value (from compute_method_thresholds)
#' @param default_strategy_value    Character scalar: "Union", "Cross-method", or "Cross-method per-trait"
#' @param snp_clumping_distance_value Integer scalar: current clumping distance (bp)
#' @param input_prefix              Character scalar: prefix for all Shiny input IDs inside this bar.
#'   Use "" (default) for a single filter bar; use e.g. "gea_" or "gwas_" when two bars
#'   coexist in the same module to avoid input ID collisions.
#' @param selected_pairs   Character vector of "trait::method" pairs currently
#'   active, or NULL to default every non-empty cell to active (first render /
#'   no persisted selection yet). Seeds which cells render `tm-active` so a
#'   caller-side re-render (e.g. after a threshold edit) does not silently
#'   reset the user's on/off choices -- see mod_gea.R's tm_selection_rv.
#' @param overrides         Named list "trait::method" -> list(type=,value=),
#'   the cells that deviate from the registry-or-master fallback.
#' @param registry_defaults Named list method -> list(adjust=,threshold=,family=),
#'   from gea_method_significance_defaults() -- pins non-univariate methods
#'   (RDA) to their registry rule unless a cell override exists.
#' @param master_type / master_value  The threshold bar's current master rule,
#'   needed here only to resolve each cell's rule SOURCE for the visual marker
#'   (override / registry / master) -- the resolved cutoff itself already
#'   arrives via combo_thresholds.
#' @return tagList (or a muted placeholder div when no data available)
#' @noRd
build_filter_bar_ui <- function(ns, traits, methods, trait_colors,
                                combo_counts, combo_thresholds = list(),
                                default_strategy_value,
                                snp_clumping_distance_value = 100000L,
                                input_prefix = "",
                                selected_pairs = NULL,
                                overrides = list(),
                                registry_defaults = list(),
                                master_type = "bonf",
                                master_value = 0.05) {
    # Helper: produce an input id with optional prefix
    pid <- function(name) ns(paste0(input_prefix, name))
    # Guard only for truly empty project (no data at all)
    if (length(methods) == 0 || length(traits) == 0) {
        return(htmltools::div(
            class = "text-muted small fst-italic py-2",
            "No association results available."
        ))
    }

    # Build table header row: blank + method column headers (each with a
    # hover-revealed gear icon that opens the bulk "all traits, this method" popup)
    header_cells <- lapply(methods, function(m) {
        htmltools::tags$th(
            class = "tm-col-header",
            `data-method` = m,
            m,
            htmltools::tags$span(
                class = "tm-gear",
                title = paste0("Set threshold for every trait — ", m),
                bsicons::bs_icon("gear-fill", size = "0.65em")
            )
        )
    })
    header_row <- htmltools::tags$tr(
        htmltools::tags$th(),  # blank corner
        header_cells
    )

    # Build body rows: one per trait
    body_rows <- lapply(traits, function(t) {
        dot_html <- paste0(
            '<span class="filter-trait-dot" style="background:', trait_colors[t],
            ';display:inline-block;width:8px;height:8px;border-radius:50%;',
            'margin-right:4px;"></span>'
        )
        row_header <- htmltools::tags$th(
            class = "tm-row-header",
            `data-trait` = t,
            htmltools::HTML(paste0(dot_html, htmltools::htmlEscape(t))),
            htmltools::tags$span(
                class = "tm-gear",
                title = paste0("Set threshold for every method — ", t),
                bsicons::bs_icon("gear-fill", size = "0.65em")
            )
        )
        cells <- lapply(methods, function(m) {
            key      <- paste0(t, "::", m)
            n        <- combo_counts[[key]] %||% 0L
            thr      <- combo_thresholds[[key]]
            thr_txt  <- if (is.null(thr) || is.na(thr)) "n/a" else
                            formatC(thr, format = "e", digits = 1)
            rule_src <- effective_rule_for(method = m, trait = t, overrides = overrides,
                                           master_type = master_type, master_value = master_value,
                                           registry_defaults = registry_defaults)$source
            thr_class <- if (rule_src == "master") "tm-thr" else paste0("tm-thr tm-thr-", rule_src)
            is_active <- n > 0 && (is.null(selected_pairs) || key %in% selected_pairs)
            if (n > 0) {
                htmltools::tags$td(
                    htmltools::tags$button(
                        class = paste("tm-cell", if (is_active) "tm-active" else ""),
                        `data-trait` = t,
                        `data-method` = m,
                        htmltools::tags$span(class = "tm-count", as.character(n)),
                        htmltools::tags$span(class = thr_class, thr_txt)
                    )
                )
            } else {
                htmltools::tags$td(
                    htmltools::tags$button(
                        class = "tm-cell tm-empty",
                        `data-trait` = t,
                        `data-method` = m,
                        htmltools::tags$span(class = "tm-count", "—"),
                        htmltools::tags$span(class = thr_class, thr_txt)
                    )
                )
            }
        })
        htmltools::tags$tr(row_header, cells)
    })

    matrix_table <- htmltools::tags$table(
        class = "tm-matrix",
        htmltools::tags$thead(header_row),
        htmltools::tags$tbody(body_rows)
    )

    # Hidden input bridge updated by JS (selection) and the popup trigger bridge
    input_id <- pid("tm_selection")
    popup_id <- pid("thr_popup")
    hidden_input <- shiny::textInput(input_id, label = NULL, value = "")

    # JS: delegated click on matrix container. Popup branches (cell threshold
    # text, row/column gear) are checked BEFORE the existing toggle-all
    # branches and each `return`s -- a single listener, not two competing
    # ones, so there is no ordering ambiguity between "open popup" and
    # "toggle selection" on the same click.
    container_id <- pid("tm_container")
    js_code <- sprintf('
(function() {
    var container = document.getElementById("%s");
    if (!container) return;
    // window-keyed last-sent cache so syncSelection() is idempotent across renderUI
    // re-creations. Each re-run of the renderUI builds a fresh IIFE (new closure), so a
    // closure-local lastSent would reset to undefined — window persists across re-renders.
    if (!window.__tmSelLast) window.__tmSelLast = {};
    function syncSelection() {
        var active = container.querySelectorAll(".tm-cell.tm-active");
        var pairs = Array.from(active).map(function(el) {
            return el.dataset.trait + "::" + el.dataset.method;
        });
        var json = JSON.stringify(pairs);
        // Skip setInputValue when selection unchanged — prevents reactive storms on
        // renderUI re-runs caused by cold-load param restoration (no {priority:"event"};
        // default priority adds a second dedup layer for genuine user interactions).
        if (window.__tmSelLast["%s"] === json) return;
        window.__tmSelLast["%s"] = json;
        Shiny.setInputValue("%s", json);
    }
    container.addEventListener("click", function(e) {
        // ── Popup triggers (checked first) ──────────────────────────────
        var thr = e.target.closest(".tm-thr");
        if (thr) {
            var cell = thr.closest(".tm-cell");
            if (cell) {
                Shiny.setInputValue("%s", {
                    scope: "cell", trait: cell.dataset.trait, method: cell.dataset.method,
                    ts: Date.now()
                }, {priority: "event"});
            }
            return;
        }
        var gear = e.target.closest(".tm-gear");
        if (gear) {
            var gth = gear.closest("th");
            if (gth.classList.contains("tm-row-header")) {
                Shiny.setInputValue("%s", {
                    scope: "trait", trait: gth.dataset.trait, method: null, ts: Date.now()
                }, {priority: "event"});
            } else if (gth.classList.contains("tm-col-header")) {
                Shiny.setInputValue("%s", {
                    scope: "method", trait: null, method: gth.dataset.method, ts: Date.now()
                }, {priority: "event"});
            }
            return;
        }
        // ── Existing toggle-selection behaviour ─────────────────────────
        var cellBtn = e.target.closest(".tm-cell");
        if (cellBtn && !cellBtn.classList.contains("tm-empty")) {
            cellBtn.classList.toggle("tm-active");
            syncSelection();
            return;
        }
        var rh = e.target.closest(".tm-row-header");
        if (rh) {
            var trait = rh.dataset.trait;
            var cells = container.querySelectorAll(".tm-cell[data-trait=\\"" + trait + "\\"]:not(.tm-empty)");
            var allActive = Array.from(cells).every(function(c) { return c.classList.contains("tm-active"); });
            cells.forEach(function(c) { c.classList.toggle("tm-active", !allActive); });
            syncSelection();
            return;
        }
        var ch = e.target.closest(".tm-col-header");
        if (ch) {
            var method = ch.dataset.method;
            var cells = container.querySelectorAll(".tm-cell[data-method=\\"" + method + "\\"]:not(.tm-empty)");
            var allActive = Array.from(cells).every(function(c) { return c.classList.contains("tm-active"); });
            cells.forEach(function(c) { c.classList.toggle("tm-active", !allActive); });
            syncSelection();
        }
    });
    syncSelection();
})();
', container_id, input_id, input_id, input_id, popup_id, popup_id, popup_id)

    strategy_choices <- c("Union", "Cross-method", "Cross-method per-trait")
    strategy_labels  <- setNames(
        strategy_choices,
        c(
            "Union — all sig SNPs from all methods",
            "Cross-method — SNPs significant in ≥2 methods within clumping distance",
            "Cross-method per-trait — same as above, same trait only"
        )
    )

    default_strat_norm <- .normalize_strategy(default_strategy_value)
    if (!default_strat_norm %in% strategy_choices) default_strat_norm <- "Union"

    n_overridden <- length(overrides)

    htmltools::div(
        class = "manhattan-filter-bar",
        htmltools::div(
            class = "filter-row align-items-start",
            # Matrix
            htmltools::div(
                class = "d-flex flex-column me-4",
                if (n_overridden > 0) htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-1",
                    htmltools::tags$span(class = "badge bg-warning text-dark",
                                         sprintf("%d cell%s overridden", n_overridden,
                                                 if (n_overridden == 1) "" else "s")),
                    shiny::actionButton(
                        pid("reset_overrides"), "Reset overrides",
                        class = "btn btn-link btn-sm text-muted p-0"
                    )
                ),
                htmltools::div(
                    id    = container_id,
                    class = "tm-container",
                    matrix_table,
                    hidden_input
                )
            ),
            # Strategy
            htmltools::div(
                class = "d-flex flex-column me-4",
                htmltools::span("Strategy", class = "filter-label mb-1"),
                shiny::radioButtons(
                    pid("combine_strategy"), label = NULL,
                    choices  = strategy_labels,
                    selected = default_strat_norm,
                    inline   = FALSE
                ),
                htmltools::span(
                    class = "text-muted small mt-1",
                    "Cross-method shows only SNPs/windows significant in ≥2 methods ",
                    "within the clumping distance. A hit significant in a single method ",
                    "appears only under Union."
                )
            ),
            # Clumping distance
            htmltools::div(
                class = "d-flex flex-column me-4",
                htmltools::div(
                    class = "mb-1",
                    filter_note(
                        "Clumping distance (bp)",
                        paste0(
                            "Region merge distance. Pick a value informed by your data's LD decay ",
                            "(Structure tab → LD Decay panel) or set a custom number."
                        )
                    )
                ),
                shiny::numericInput(
                    pid("snp_clumping_distance"), label = NULL,
                    value = snp_clumping_distance_value, min = 1000L, step = 100000L,
                    width = "160px"
                ),
                htmltools::span(
                    class = "text-muted small mt-1",
                    "Merge distance for regions; also used for Cross-method overlap"
                )
            )
        ),
        htmltools::tags$script(htmltools::HTML(js_code))
    )
}

#' Compute interactive sig SNPs from per-method data + matrix selection
#'
#' Pure function — call inside reactive(). Filters, combines, and stamps region IDs.
#'
#' @param all_method_sigsnps Named list of data.tables (method -> dt with SNPID/chr/pos/pvalue/method/trait)
#' @param tm_selection_json  JSON string from tm_selection hidden input (list of "trait::method" pairs)
#' @param combo_counts       Named list: "trait::method" -> integer count (for fallback)
#' @param known_traits       Character vector of traits to restrict to (avoids cross-tab bleed)
#' @param strategy           Character: "Union", "Cross-method", "Cross-method per-trait", or method name
#' @param clumping_distance  Integer: SNP clumping distance in bp
#' @param project_name       Character: project name (for assign_region_ids)
#' @param module             Character: MOD_GEA, MOD_GWAS, or MOD_GEAXGWAS
#' @return data.table (SNPID/chr/pos/pvalue/method/trait/region_id) or .empty_sigsnps_assoc()
#' @noRd
compute_interactive_sigsnps <- function(all_method_sigsnps, tm_selection_json,
                                        combo_counts, known_traits,
                                        strategy, clumping_distance,
                                        project_name, module) {
    if (length(all_method_sigsnps) == 0) return(NULL)
    strategy <- .normalize_strategy(strategy)

    # All valid trait::method pairs (positive count, restricted to known_traits)
    all_valid_pairs <- names(combo_counts)[
        vapply(combo_counts, function(n) n > 0, logical(1))
    ]
    all_valid_pairs <- all_valid_pairs[
        sub("::.*", "", all_valid_pairs) %in% known_traits
    ]

    # Parse selected pairs from JSON; fall back to all valid.
    # An empty JSON array ("[]") from JS means the matrix rendered with no
    # tm-active cells (un-initialized state), not a deliberate user deselect.
    # Treat length-0 parsed result the same as NULL / "" to show all by default.
    selected_pairs <- if (is.null(tm_selection_json) || !nzchar(tm_selection_json)) {
        all_valid_pairs
    } else {
        parsed <- tryCatch(
            jsonlite::fromJSON(tm_selection_json),
            error = function(e) NULL
        )
        if (is.null(parsed) || length(parsed) == 0) all_valid_pairs else parsed
    }

    # Filter each method to selected traits for that method
    filtered <- lapply(names(all_method_sigsnps), function(m) {
        paired_traits <- sub("::.*", "",
            selected_pairs[grepl(paste0("::", m, "$"), selected_pairs, fixed = FALSE)])
        if (length(paired_traits) == 0) return(NULL)
        dt <- all_method_sigsnps[[m]]
        if (nrow(dt) == 0) return(NULL)
        dt[trait %in% paired_traits]
    })
    names(filtered) <- names(all_method_sigsnps)
    filtered <- Filter(function(x) !is.null(x) && nrow(x) > 0, filtered)

    if (length(filtered) == 0) return(.empty_sigsnps_assoc())

    combined_dt <- combine_sigsnps(filtered, strategy = strategy,
                                   gap = as.integer(clumping_distance))

    if (nrow(combined_dt) == 0) return(.empty_sigsnps_assoc())

    # region_id is stamped downstream from the live-computed regions (see
    # assign_region_ids_from_regions() / plotted_sigsnps in mod_gea.R / mod_gwas.R),
    # not from the static pipeline regions_combined.tsv — keeps the clicked SNP's
    # region in sync with the interactive region table/rectangles.
    combined_dt[, region_id := NA_character_]
    combined_dt
}

#' Wire the trait×method matrix's rule popup + reset-overrides + apply-to-
#' config observers for ONE filter bar instance.
#'
#' Shared by mod_gea.R (GEA tab), mod_gwas.R (GWAS tab), and mod_gea_x_gwas.R
#' (called TWICE, once per side — hence `input_prefix` distinguishing "gea_"
#' from "gwas_" so the two popups/reset/apply observers never collide within
#' one moduleServer). See R/fct_threshold_rules.R for the precedence the
#' popup edits (cell override > registry default > master).
#'
#' @param input,output,session Standard Shiny module server args (the CALLING
#'   module's — this is not itself a moduleServer, just a helper that
#'   registers observers into the caller's reactive graph)
#' @param ns              The calling module's session$ns
#' @param input_prefix    "" for a single filter bar, "gea_"/"gwas_" for two
#'   sharing one moduleServer (GEAxGWAS)
#' @param project_data    Reactive project data bundle
#' @param module          MOD_GEA / MOD_GWAS / MOD_GEAXGWAS — region_params.json
#'   partition (GEAxGWAS's two sides still need input_prefix-scoped KEYS within
#'   that one module's params, since both sides share module = MOD_GEAXGWAS)
#' @param methods,all_traits  Reactives: character vectors for this side
#' @param threshold_type,threshold_value  Reactives: this side's master rule
#' @param registry_defaults  Reactive: method -> list(adjust=,threshold=,family=)
#' @param config_rules       Reactive: method -> "bonf_0.05" (pipeline rule string)
#' @param config_state       reactiveValues from app_server.R, or NULL to omit
#'   "Apply rules to config" (GEAxGWAS's two bars don't drive a pipeline run)
#' @param config_module_key  "GEA" or "GWAS" — which config_state$working key
#'   Apply-to-config collapses rules into
#' @return list(overrides = debounced reactive "trait::method"->rule,
#'   selected_pairs = reactive character vector or NULL, reload = function()
#'   that re-reads both from region_params.json into this instance's
#'   reactiveVals — call after writing into the SAME JSON keys from outside
#'   this instance, e.g. mod_gea_x_gwas.R's "Fill from …" buttons) — the
#'   first two feed straight into build_filter_bar_ui()
#' @noRd
setup_matrix_rules_server <- function(input, output, session, ns, input_prefix = "",
                                      project_data, module, methods, all_traits,
                                      threshold_type, threshold_value,
                                      registry_defaults, config_rules,
                                      config_state = NULL, config_module_key = "GEA") {
    pid <- function(name) paste0(input_prefix, name)

    # ── Per-cell overrides ("trait::method" -> rule) ────────────────────────
    threshold_overrides_rv <- shiny::reactiveVal(list())

    # ── Matrix selection (which cells are included in the combined view) ───
    # Separate from threshold_overrides — this is WHICH cells, not their rule.
    # NULL means "not yet set"; compute_interactive_sigsnps() then defaults to
    # every non-empty cell active (pre-rework behaviour).
    tm_selection_rv <- shiny::reactiveVal(NULL)

    # Re-read both from region_params.json — cold load on project switch, and
    # exposed as `reload()` so a caller who wrote INTO these same JSON keys
    # from outside this instance (mod_gea_x_gwas.R's "Fill from …" buttons,
    # which copy the OTHER side's persisted rules) can force this instance's
    # reactiveVals to pick up the change immediately, not just on next reload.
    reload_from_region_params <- function() {
        pd <- project_data(); if (is.null(pd)) return()
        rp <- read_region_params(pd$name)
        threshold_overrides_rv(normalize_threshold_overrides(
            get_global_param(rp, module, paste0(input_prefix, "threshold_overrides")),
            traits = all_traits(), methods = methods()
        ))
        saved <- get_global_param(rp, module, paste0(input_prefix, "snp_matrix_selection"))
        tm_selection_rv(if (is.null(saved)) NULL else as.character(unlist(saved)))
    }
    shiny::observe(reload_from_region_params())

    # Debounced — an override edit shouldn't fire a full-vector scan per keystroke.
    threshold_overrides <- shiny::debounce(
        shiny::reactive(threshold_overrides_rv()), 500)

    persist_overrides <- function(ov) {
        pd <- project_data(); if (is.null(pd)) return()
        threshold_overrides_rv(ov)
        rp <- read_region_params(pd$name)
        # set_global_param(..., NULL) deletes the key — a full reset leaves the
        # file byte-identical to the legacy (pre-override) shape.
        rp <- set_global_param(rp, module, paste0(input_prefix, "threshold_overrides"),
                               if (length(ov) == 0) NULL else ov)
        save_region_params(pd$name, rp)
    }

    shiny::observeEvent(input[[pid("tm_selection")]], {
        pd <- project_data(); if (is.null(pd)) return()
        raw    <- input[[pid("tm_selection")]]
        parsed <- tryCatch(jsonlite::fromJSON(raw), error = function(e) NULL)
        pairs  <- if (is.null(parsed) || length(parsed) == 0) NULL else as.character(parsed)
        tm_selection_rv(pairs)
        rp <- read_region_params(pd$name)
        rp <- set_global_param(rp, module, paste0(input_prefix, "snp_matrix_selection"), pairs)
        save_region_params(pd$name, rp)
    }, ignoreInit = TRUE)

    # ── Rule popup: cell / row (trait) / column (method) ────────────────────
    # input[[pid("thr_popup")]] is an object payload {scope, trait, method, ts}
    # set by fct_combine.R's delegated click listener with {priority:"event"} —
    # opening the same cell's popup twice must not be deduped away.
    popup_ctx <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input[[pid("thr_popup")]], {
        ctx <- input[[pid("thr_popup")]]
        popup_ctx(ctx)
        ov <- threshold_overrides_rv(); rd <- registry_defaults()
        mt <- threshold_type(); mv <- threshold_value()

        title <- switch(ctx$scope,
            cell   = sprintf("%s × %s", ctx$trait, ctx$method),
            trait  = sprintf("%s — all methods", ctx$trait),
            method = sprintf("%s — all traits", ctx$method),
            "Significance rule"
        )
        # A cell popup resolves its own cell; a row/column popup seeds from a
        # representative (method[1]/no-trait or trait/no-method) rule — Apply
        # always writes an explicit rule to every cell in scope regardless of
        # what they currently show, so the seed is a starting point, not a summary.
        seed_method <- ctx$method %||% methods()[1]
        seed_trait  <- ctx$trait  %||% NULL
        rule <- effective_rule_for(method = seed_method, trait = seed_trait,
                                   overrides = ov, master_type = mt, master_value = mv,
                                   registry_defaults = rd)
        cfg_rule <- if (!is.null(ctx$method)) config_rules()[[ctx$method]] %||% "" else ""
        source_label <- switch(rule$source,
            override = "your override", registry = "method's registry default",
            "master rule")

        shiny::showModal(shiny::modalDialog(
            title = title,
            shiny::selectInput(ns(pid("thr_popup_type")), "Rule",
                choices = c("Bonferroni" = "bonf", "FDR (qval)" = "qval",
                           "Top N" = "top", "Custom (raw p)" = "custom"),
                selected = rule$type),
            shiny::numericInput(ns(pid("thr_popup_value")), "Value",
                                value = rule$value, min = 0, step = NA),
            htmltools::p(class = "small text-muted mb-1",
                sprintf("Currently: %s %s (%s)", rule$type, rule$value, source_label)),
            if (nzchar(cfg_rule)) htmltools::p(class = "small text-muted",
                sprintf("Pipeline files (QQ/Manhattan background) were built with %s.", cfg_rule)),
            footer = htmltools::tagList(
                shiny::actionButton(ns(pid("thr_popup_follow")), "Follow master/registry",
                                    class = "btn btn-link btn-sm text-muted"),
                shiny::modalButton("Cancel"),
                shiny::actionButton(ns(pid("thr_popup_apply")), "Apply", class = "btn btn-primary")
            ),
            easyClose = TRUE
        ))
    })

    scope_keys <- function(ctx) {
        if (is.null(ctx)) return(character(0))
        switch(ctx$scope,
            cell   = paste0(ctx$trait, "::", ctx$method),
            trait  = paste0(ctx$trait, "::", methods()),
            method = paste0(all_traits(), "::", ctx$method),
            character(0)
        )
    }

    shiny::observeEvent(input[[pid("thr_popup_apply")]], {
        ctx <- popup_ctx(); if (is.null(ctx)) return()
        keys <- scope_keys(ctx)
        ty <- input[[pid("thr_popup_type")]]; va <- input[[pid("thr_popup_value")]]
        if (!threshold_value_valid_for_type(ty, va)) {
            shiny::showNotification("Invalid threshold value for this rule.", type = "error")
            return()
        }
        ov <- threshold_overrides_rv()
        for (k in keys) ov[[k]] <- list(type = ty, value = as.numeric(va))
        persist_overrides(ov)
        shiny::removeModal()
    })

    shiny::observeEvent(input[[pid("thr_popup_follow")]], {
        ctx <- popup_ctx(); if (is.null(ctx)) return()
        keys <- scope_keys(ctx)
        ov <- threshold_overrides_rv()
        ov[keys] <- NULL
        persist_overrides(ov)
        shiny::removeModal()
    })

    shiny::observeEvent(input[[pid("reset_overrides")]], {
        persist_overrides(list())
        shiny::showNotification("All threshold overrides reset to master/registry.",
                                type = "message", duration = 3)
    })

    # ── Apply rules to config ────────────────────────────────────────────────
    # Collapses master/registry/override into one adjust+threshold per method
    # and writes it into config_state$working. The pipeline has no per-trait
    # granularity (GEA.configs/GWAS.configs is per-method — see
    # workflow/rules/common.smk's parse_association_configs()), so a method
    # whose cells disagree needs an explicit choice, never a silent pick.
    if (!is.null(config_state)) {
        shiny::observeEvent(input[[pid("apply_to_config")]], {
            ov <- threshold_overrides_rv(); mt <- threshold_type(); mv <- threshold_value()
            rd <- registry_defaults(); ms <- methods(); tr <- all_traits()

            per_method <- lapply(ms, function(m) {
                rules <- lapply(tr, function(t) effective_rule_for(
                    method = m, trait = t, overrides = ov,
                    master_type = mt, master_value = mv, registry_defaults = rd))
                keys <- unique(vapply(rules, function(r) paste0(r$type, "_", r$value), character(1)))
                list(method = m, keys = keys, rule = rules[[1]])
            })
            conflicts <- Filter(function(x) length(x$keys) > 1, per_method)

            if (length(conflicts) == 0) {
                # Carry each method's existing `params` (K, n_pcs, axes, …)
                # forward untouched — this button only ever sets adjust/
                # threshold. Losing hyperparameters here would be exactly the
                # silent-discard mistake the sidebar's own reconcile (see
                # mod_config_sidebar.R) was written to avoid on the other side.
                cur_val <- config_get_by_path(config_state$working, paste0(config_module_key, ".configs"))
                cur_by_method <- stats::setNames(
                    cur_val %||% list(),
                    vapply(cur_val %||% list(), function(x) x$method %||% "", character(1))
                )
                new_configs <- lapply(per_method, function(x) {
                    existing <- cur_by_method[[x$method]]
                    list(method = x$method, adjust = x$rule$type,
                        threshold = as.character(x$rule$value),
                        params = existing$params %||% list())
                })
                config_state$working <- config_set_by_path(
                    config_state$working, paste0(config_module_key, ".configs"), new_configs)
                shiny::showNotification(
                    "Rules applied to config (working copy) — Run or Save Project Files to persist to YAML.",
                    type = "message", duration = 5)
            } else {
                names_conflict <- paste(vapply(conflicts, `[[`, character(1), "method"), collapse = ", ")
                shiny::showNotification(
                    sprintf(paste(
                        "%s %s per-trait rules that disagree — the pipeline config is",
                        "per-method only. Give %s a single rule (via its column popup)",
                        "before applying."
                    ), names_conflict, if (length(conflicts) == 1) "has" else "have",
                       if (length(conflicts) == 1) "it" else "them"),
                    type = "warning", duration = 8)
            }
        })
    }

    list(
        overrides      = threshold_overrides,
        selected_pairs = shiny::reactive(tm_selection_rv()),
        reload         = reload_from_region_params
    )
}
