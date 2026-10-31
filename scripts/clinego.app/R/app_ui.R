#' Build the application UI
#' @noRd
app_ui <- function(request) {
    pipeline_path <- get_pipeline_path()
    projects <- find_selectable_projects(pipeline_path)

    # Serve pipeline output files as static resources
    shiny::addResourcePath("pipeline", pipeline_path)

    # ...and the package's own assets. The header below loads two scripts as
    # src = "www/...", which resolves only if a "www" resource prefix exists.
    # dev.R registers one (pointing at the bind mount); the PACKAGE path never
    # did, so clinego.app::run_app() served 404s for config-dirty.js and
    # module-bar.js — and module-bar.js is what drives tab switching, since
    # dashboard.scss hides the real navbar links. The SCSS is unaffected: it goes
    # through sass_file(app_sys(...)) at theme build time, not over HTTP.
    #
    # Registered only when absent, so dev.R's mount-backed prefix keeps winning
    # in dev mode — otherwise the installed (stale) copy of the JS would shadow
    # the file being edited.
    if (!"www" %in% names(shiny::resourcePaths())) {
        www_dir <- app_sys("app", "www")
        if (!nzchar(www_dir))
            www_dir <- file.path(pipeline_path, "scripts", "clinego.app",
                                 "inst", "app", "www")
        if (dir.exists(www_dir)) {
            shiny::addResourcePath("www", www_dir)
        } else {
            warning("app_ui: no www/ asset directory found; ",
                    "config-dirty.js and module-bar.js will not load")
        }
    }

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
            # Dark mode only (A4 I2): the attribute Bootstrap's dark colour mode keys on,
            # set before the body renders (htmltools hoists tags$head into <head>), so
            # the page never paints light first. There is no light theme and no toggle.
            htmltools::tags$head(htmltools::tags$script(htmltools::HTML(
                "document.documentElement.setAttribute('data-bs-theme', 'dark');"))),
            shinyjs::useShinyjs(),
            shiny::useBusyIndicators(spinners = TRUE),
            htmltools::tags$script(src = "www/config-dirty.js"),
            htmltools::tags$script(src = "www/module-bar.js"),

            # ONE bar at the top: brand | project | module tabs | Save.
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
                    lab_input_width(shiny::selectInput(
                        "project_selector", label = NULL,
                        choices = c(
                            "+ New Project" = "__new__",
                            if (length(projects) > 0) setNames(projects, projects)
                            else character(0)
                        ),
                        selected = if (length(projects) > 0) projects[1] else "__new__"
                    ), 200)
                ),
                htmltools::div(class = "lab-bar-modules",
                               shiny::uiOutput("module_bar")),
                # The one Save (I9): writes the working config to config_{project}.yaml,
                # as Run does, from any tab. config-dirty.js marks it while the
                # working config differs from the saved one.
                htmltools::div(class = "lab-bar-right",
                    shiny::actionButton(
                        "save_project_files",
                        label = htmltools::tagList(bsicons::bs_icon("floppy"), "Save"),
                        class = "btn-sm btn-outline-light lab-save-btn",
                        title = "Save the config of this project (config_<project>.yaml)"
                    ))
            )
        ),

        # ── Tabs ───────────────────────────────────────────────────────────────

        bslib::nav_panel("Home", value = "home",
            class = "bslib-page-dashboard",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui("config_home", "Project Files"),
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
                    "config_climate", "Environmental Config",
                    runner_ui = mod_pipeline_runner_ui("runner_climate",
                                                       btn_label = "Run Environmental")
                ),
                fill = FALSE, fillable = FALSE,
                mod_climate_ui("climate")
            )
        ),

        bslib::nav_panel("Phenotypic", value = "traits",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_traits", "Phenotypic Config",
                    runner_ui = mod_pipeline_runner_ui("runner_traits",
                                                       btn_label = "Run Phenotypic")
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
