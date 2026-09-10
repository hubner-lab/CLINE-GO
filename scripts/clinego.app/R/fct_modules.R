# fct_modules.R — module registry and the two-row module bar
#
# The registry is the single source of truth for top-level navigation: it drives
# the module bar, and (from Phase 2) which modules a project's regime exposes.
# The bslib nav_panels still exist for every module; the bar decides which ones
# a user can reach.

#' Module registry
#'
#' One row per top-level module, in pipeline order.
#'
#' * `id`        stable identifier — matches the `nav_panel` `value`, the
#'               config-schema tab name, and the runner's `tab_name`. Never
#'               rename: everything downstream keys off it.
#' * `label`     display text only. Safe to change without touching anything
#'               else, which is why nav panels carry an explicit `value`.
#' * `group`     bar grouping. A group with several members renders its name on
#'               the top row and its members as columns underneath; a group with
#'               a single member spans both rows in larger type.
#' * `gwas_only` whether the module survives the `gwas_only` regime.
#'
#' @noRd
module_registry <- function() {
    data.frame(
        id = c(
            "home", "processing",
            "prestructure", "structure",
            "climate", "traits",
            "pregea", "gea", "gwas", "gea_x_gwas",
            "maladaptation"
        ),
        label = c(
            "Home", "Processing",
            "PreStructure", "Structure",
            "Environmental", "Phenotypic",
            "PreGEA", "GEA", "GWAS", "GEAxGWAS",
            "Maladaptation"
        ),
        group = c(
            "Setup", "Setup",
            "Structure", "Structure",
            "Factors", "Factors",
            "Association", "Association", "Association", "Association",
            "Maladaptation"
        ),
        # Phenotypic survives gwas_only — it is the ONLY factor-characterization
        # module such a project has (mode=climate raises without climate), and
        # mode=traits needs neither coordinates nor climate.
        gwas_only = c(
            TRUE, TRUE,
            TRUE, TRUE,
            FALSE, TRUE,
            FALSE, FALSE, TRUE, FALSE,
            FALSE
        ),
        stringsAsFactors = FALSE
    )
}

#' Input id of a module bar button
#' @noRd
module_input_id <- function(id) paste0("mb_", id)

#' Module ids reachable under a regime
#' @noRd
module_ids_for_regime <- function(regime = "standard") {
    reg <- module_registry()
    if (identical(regime, "gwas_only")) reg <- reg[reg$gwas_only, , drop = FALSE]
    reg$id
}

#' Config-sidebar sections that do not apply under a regime
#'
#' Keyed by the schema entry's `section` label (fct_config_schema.R). LD Decay is
#' deliberately NOT hidden in gwas_only: with a single site it degrades to the
#' genome-wide "All" curve, which is still worth having.
#' @noRd
hidden_sidebar_sections <- function(regime = "standard") {
    if (identical(regime, "gwas_only")) c("Map", "Piemap", "Pop Stats") else character(0)
}

#' Registry rows visible under a regime, split into groups (bar order preserved)
#' @noRd
module_groups <- function(regime = "standard") {
    reg <- module_registry()
    if (identical(regime, "gwas_only")) reg <- reg[reg$gwas_only, , drop = FALSE]
    split(reg, factor(reg$group, levels = unique(reg$group)))
}

#' One clickable module cell
#' @noRd
module_bar_item <- function(id, label, active, solo = FALSE) {
    shiny::actionLink(
        inputId = module_input_id(id),
        label   = label,
        class   = paste(c(
            "mb-item",
            if (solo) "mb-item-solo",
            if (identical(id, active)) "mb-active"
        ), collapse = " ")
    )
}

#' Build the two-row module bar
#'
#' @param active id of the module currently shown
#' @param regime "standard" or "gwas_only"
#' @noRd
module_bar_ui <- function(active = "home", regime = "standard") {
    groups <- module_groups(regime)

    htmltools::div(
        class = "module-bar",
        lapply(names(groups), function(g_name) {
            members <- groups[[g_name]]

            # Single-member group: one cell spanning both rows, showing the
            # MODULE name (a lone "Factors" cell would not say which factors).
            if (nrow(members) == 1) {
                return(htmltools::div(
                    class = "mb-group mb-group-solo",
                    module_bar_item(members$id[1], members$label[1], active, solo = TRUE)
                ))
            }

            htmltools::div(
                class = "mb-group",
                htmltools::div(class = "mb-group-label", g_name),
                htmltools::div(
                    class = "mb-members",
                    Map(function(i, l) module_bar_item(i, l, active),
                        members$id, members$label)
                )
            )
        })
    )
}
