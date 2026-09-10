#' Build the application UI
#' @noRd
app_ui <- function(request) {
    pipeline_path <- get_pipeline_path()
    projects <- find_selectable_projects(pipeline_path)

    # Serve pipeline output files as static resources
    shiny::addResourcePath("pipeline", pipeline_path)

    bslib::page_navbar(
        title = NULL,
        theme    = app_theme(),
        id       = "main_tabs",
        # Matches the nav_panel VALUE, not its title (bslib resolves `fillable`
        # against values, which default to the title only when none is given —
        # every panel below now sets one explicitly).
        fillable = c("home"),

        # Enable shinyjs, load config-dirty.js, and show busy spinners on
        # recalculating outputs (plotlyOutput, tableOutput, etc.) — shiny >= 1.8.1
        #
        # `header` renders ONCE, between the navbar and the tab content — not per
        # panel — which is what makes it the right slot for the module bar. The
        # default navbar tab links are visually hidden in custom.scss; the bar
        # replaces them and drives them (see module-bar.js).
        header = htmltools::tagList(
            shinyjs::useShinyjs(),
            shiny::useBusyIndicators(spinners = TRUE),
            htmltools::tags$script(src = "www/config-dirty.js"),
            htmltools::tags$script(src = "www/module-bar.js"),

            # ONE bar at the top: brand | project | module tabs | dark toggle.
            # The app's own navbar is visually hidden by dashboard.scss (NOT
            # display:none -- module-bar.js switches tabs by firing synthetic
            # clicks on its anchors, so they must stay in the tree), and this
            # replaces both it and the separate module bar. Two stacked bars
            # cost ~104px of every module's vertical budget.
            htmltools::div(
                class = "lab-chrome-bar",
                htmltools::div(
                    class = "lab-bar-left",
                    brand_logo(),
                    shiny::selectInput(
                        "project_selector", label = NULL,
                        choices = c(
                            "+ New Project" = "__new__",
                            if (length(projects) > 0) setNames(projects, projects)
                            else character(0)
                        ),
                        selected = if (length(projects) > 0) projects[1] else "__new__",
                        width = "200px"
                    )
                ),
                htmltools::div(class = "lab-bar-modules",
                               shiny::uiOutput("module_bar")),
                # Display preference, not identity and not navigation -- so it
                # sits apart from both, at the far right.
                htmltools::div(class = "lab-bar-right",
                               bslib::input_dark_mode(id = "dark_mode", mode = "dark"))
            )
        ),

        # ── Tabs ───────────────────────────────────────────────────────────────

        bslib::nav_panel("Home", value = "home",
            class = "bslib-page-dashboard",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_home", "Project Files",
                    runner_ui = shiny::actionButton(
                        "save_project_files",
                        label = htmltools::tagList(
                            bsicons::bs_icon("floppy", size = "0.9em"),
                            " Save Project Files"
                        ),
                        class = "btn btn-primary btn-sm w-100 mt-1"
                    )
                ),
                fill = TRUE, fillable = TRUE,
                mod_home_ui("home")
            )
        ),

        bslib::nav_panel("Processing", value = "processing",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_processing", "Processing Config",
                    runner_ui = mod_pipeline_runner_ui("runner_processing",
                                                       btn_label = "Run Processing")
                ),
                fill = FALSE, fillable = FALSE,
                mod_processing_ui("processing")
            )
        ),

        bslib::nav_panel("PreStructure", value = "prestructure",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_prestructure", "PreStructure Config",
                    runner_ui = mod_pipeline_runner_ui("runner_prestructure",
                                                       btn_label = "Run PreStructure")
                ),
                fill = FALSE, fillable = FALSE,
                mod_prestructure_ui("prestructure")
            )
        ),

        bslib::nav_panel("Structure", value = "structure",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_structure", "Structure Config",
                    runner_ui = mod_pipeline_runner_ui("runner_structure",
                                                       btn_label = "Run Structure")
                ),
                fill = FALSE, fillable = FALSE,
                mod_structure_ui("structure")
            )
        ),

        bslib::nav_panel("Environmental", value = "climate",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_climate", "Climate Config",
                    runner_ui = mod_pipeline_runner_ui("runner_climate",
                                                       btn_label = "Run Climate")
                ),
                fill = FALSE, fillable = FALSE,
                mod_climate_ui("climate")
            )
        ),

        bslib::nav_panel("Phenotypic", value = "traits",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_traits", "Traits Config",
                    runner_ui = mod_pipeline_runner_ui("runner_traits",
                                                       btn_label = "Run Traits")
                ),
                fill = FALSE, fillable = FALSE,
                mod_traits_ui("traits")
            )
        ),

        bslib::nav_panel("PreGEA", value = "pregea",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_pregea", "PreGEA Config",
                    runner_ui = mod_pipeline_runner_ui("runner_pregea",
                                                       btn_label = "Run PreGEA")
                ),
                fill = FALSE, fillable = FALSE,
                mod_pregea_ui("pregea")
            )
        ),

        bslib::nav_panel("GEA", value = "gea",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_gea", "GEA Config",
                    runner_ui = mod_pipeline_runner_ui("runner_gea",
                                                       btn_label = "Run GEA")
                ),
                fill = FALSE, fillable = FALSE,
                mod_gea_ui("gea")
            )
        ),

        bslib::nav_panel("GWAS", value = "gwas",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_gwas", "GWAS Config",
                    runner_ui = mod_pipeline_runner_ui("runner_gwas",
                                                       btn_label = "Run GWAS")
                ),
                fill = FALSE, fillable = FALSE,
                mod_gwas_ui("gwas")
            )
        ),

        bslib::nav_panel("GEAxGWAS", value = "gea_x_gwas",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_gea_x_gwas", "GEAxGWAS Config",
                    runner_ui = mod_pipeline_runner_ui("runner_gea_x_gwas",
                                                       btn_label = "Run GEAxGWAS")
                ),
                fill = FALSE, fillable = FALSE,
                mod_gea_x_gwas_ui("gea_x_gwas")
            )
        ),

        bslib::nav_panel("Maladaptation", value = "maladaptation",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_maladaptation", "Maladaptation Config",
                    runner_ui = mod_pipeline_runner_ui("runner_maladaptation",
                                                       btn_label = "Run Maladaptation")
                ),
                fill = FALSE, fillable = FALSE,
                mod_maladaptation_ui("maladaptation")
            )
        )
    )
}
