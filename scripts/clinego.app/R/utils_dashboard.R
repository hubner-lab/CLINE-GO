# utils_dashboard.R — shared dashboard primitives.
#
# Promoted from scripts/layout_lab/R/lab_10_kit.R on 2026-08-29.
#
# Only helpers that TWO OR MORE of the three independent designs asked for.
# The kit is what the designs converged on, not a vocabulary imposed up front.

#' The eight image-card ids the untouched mod_processing_server() binds.
#' Single source of truth so three hand-written UIs cannot drift apart.
PROCESSING_QC_CARD_IDS <- c(
    "attrition", "sample_miss", "het_miss", "relatedness",
    "relatedness_mds", "snp_miss", "maf", "snp_density"
)

#' Every non-image output id in the server contract.
PROCESSING_OTHER_IDS <- c(
    "summary_boxes", "lea_losses_alert", "filtering_table",
    "het_table", "dl_het", "depth_section"
)

#' Variant scoping wrapper -> <div class="lab-v-a lab-mod-processing"> etc.
#'
#' `module` matters: without it, a height rule written for one module (e.g.
#' climate's two stacked heroes) can out-specify another module's rule and
#' silently resize it. That regression happened once -- Processing's hero went
#' 453px -> 276px when climate added `.lab-hero-col > .lab-hero > .card`.
lab_root <- function(variant, ..., module = NULL) {
    cls <- paste0("lab-v-", variant)
    if (!is.null(module)) cls <- paste0(cls, " lab-mod-", module)
    htmltools::div(class = cls, ...)
}

#' KPI ribbon: summary_boxes grows, the LEA badge sits beside it instead of
#' claiming its own row (~28px reclaimed). The server already right-aligns the
#' badge with its own justify-content-end wrapper -- do not add another.
#' @param alert_id id of a right-aligned badge output, or NULL when the module has
#'   none (Processing has the LEA-losses badge; PreStructure has nothing there).
lab_kpi_row <- function(ns, boxes_id = "summary_boxes", alert_id = "lea_losses_alert") {
    htmltools::div(
        class = "lab-kpi-row d-flex align-items-center gap-2",
        htmltools::div(class = "flex-grow-1 lab-kpi-boxes",
                       shiny::uiOutput(ns(boxes_id))),
        if (!is.null(alert_id))
            htmltools::div(class = "lab-kpi-alert",
                           shiny::uiOutput(ns(alert_id)))
    )
}

#' A KPI value box in the house style. Wraps bslib::value_box so every module's
#' strip looks identical without repeating showcase/theme boilerplate.
lab_value_box <- function(title, value, theme = "primary", icon = "dot") {
    bslib::value_box(
        title    = title,
        value    = value,
        theme    = theme,
        showcase = bsicons::bs_icon(icon)
    )
}

#' Uppercase domain divider (Sample QC / SNP QC).
#'
#' `aside` is an optional right-aligned slot for a control that governs the
#' section it labels -- PreStructure's K selector, for instance. A selector that
#' drives exactly one section does not need a `.control-bar` row of its own
#' (measured: a 128px select sitting in a 64px band), and putting it on the
#' divider keeps the control and the plots it changes in the same visual unit.
#' Omitted by every other caller, so the header renders exactly as before.
lab_section_header <- function(label, icon = NULL, aside = NULL) {
    htmltools::div(
        class = if (is.null(aside)) "lab-section-header"
                else "lab-section-header lab-section-header--aside",
        if (!is.null(icon)) bsicons::bs_icon(icon),
        htmltools::span(label),
        if (!is.null(aside))
            htmltools::div(class = "lab-section-aside ms-auto", aside)
    )
}

#' CSS-grid tiling. `cols` becomes a custom property so column count is data.
lab_thumb_grid <- function(..., cols = 4) {
    htmltools::div(
        class = "lab-thumb-grid",
        style = paste0("--lab-cols: ", cols, ";"),
        ...
    )
}

#' One grid item. `span` matters: snp_density is a chromosome strip chart and is
#' unreadable in a square tile.
lab_thumb <- function(..., span = 1) {
    htmltools::div(
        class = if (span > 1) "lab-thumb lab-thumb--wide" else "lab-thumb",
        style = paste0("grid-column: span ", span, ";"),
        ...
    )
}

#' Filtering Summary card — the SNP/sample rail table. Ids spelled once.
lab_filtering_summary_card <- function(ns) {
    bslib::card(
        full_screen = TRUE,
        class = "lab-dock-card lab-dock-filtering",
        bslib::card_header(bsicons::bs_icon("table"), " Filtering Summary"),
        bslib::card_body(DT::DTOutput(ns("filtering_table")))
    )
}

#' Sample Heterozygosity card — bundles het_table + dl_het so neither id is
#' ever retyped by hand.
lab_het_table_card <- function(ns) {
    bslib::card(
        full_screen = TRUE,
        class = "lab-dock-card lab-dock-het",
        bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            htmltools::span(bsicons::bs_icon("table"), " Sample Heterozygosity"),
            shiny::downloadButton(ns("dl_het"), "CSV",
                                  class = "btn-sm btn-outline-secondary")
        ),
        bslib::card_body(DT::DTOutput(ns("het_table")))
    )
}

#' depth_section in an accordion below the fold.
#'
#' `open` is a PARAMETER, not a constant: A and C want it closed (a collapsed
#' panel is display:none, so Shiny suspends the output and the PNG is never
#' transmitted -- verified by the Step 0 smoke test). B wants it open, because
#' B's whole premise is that nothing is ever hidden.
lab_depth_below_fold <- function(ns, open = FALSE) {
    bslib::accordion(
        class = "lab-below-fold",
        open  = open,
        bslib::accordion_panel(
            value = "depth",
            title = htmltools::tagList(
                bsicons::bs_icon("layers"), " Depth Distribution ",
                htmltools::span(class = "text-muted small",
                                "(context — often empty for GBS)")
            ),
            shiny::uiOutput(ns("depth_section"))
        )
    )
}

#' Anchor-figure wrapper: <div class="lab-hero">.
#' @noRd
lab_hero <- function(...) htmltools::div(class = "lab-hero", ...)
