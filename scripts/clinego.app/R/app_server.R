#' Main application server
#' @noRd
app_server <- function(input, output, session) {

    # ── Config state: working copy + saved (YAML) ──────────────────────────────
    # config_state$saved   = last-written YAML config (matches pipeline results)
    # config_state$working = user's edits (lost on reload unless Run is clicked)
    # config_state$project = current project name
    config_state <- shiny::reactiveValues(
        saved   = list(),
        working = list(),
        project = NULL
    )

    # Track the last valid project for reverting selector after modal cancel
    prev_project <- shiny::reactiveVal(NULL)

    # Benchmark mode is gone: the merged chrome has no such switch, and a
    # project declares whether it is simulated in its own config, so it was
    # never a session mode. find_selectable_projects() returns one list of both
    # pools -- see fct_discovery.R for why the old partition stranded them.

    # Populate config_state on project switch; intercept __new__ to open modal
    shiny::observeEvent(input$project_selector, {
        shiny::req(input$project_selector)

        if (input$project_selector == "__new__") {
            # Revert selector to the previous valid project
            prev <- prev_project()
            if (!is.null(prev)) {
                shiny::updateSelectInput(session, "project_selector", selected = prev)
            }
            # Show creation modal
            shiny::showModal(mod_create_project_ui("create_project"))
            return()
        }

        proj <- input$project_selector
        prev_project(proj)
        cfg  <- read_project_config(proj)
        config_state$project <- proj
        config_state$saved   <- cfg
        config_state$working <- cfg

        # The bar re-renders for the new regime, but the active pane does not move
        # on its own -- switching from a standard project while sitting on GEA
        # would leave the GEA pane visible with no button for it anywhere. Send the
        # user home whenever the current module is out of regime.
        if (!active_module() %in% module_ids_for_regime(config_regime(cfg))) {
            active_module("home")
            session$sendCustomMessage("module_bar_select", list(value = "home"))
        }
    }, ignoreInit = FALSE)

    # ── Module bar ────────────────────────────────────────────────────────────
    # Replaces the navbar's flat tab strip with modules grouped by analysis
    # stage. The nav links still exist (hidden in custom.scss) and the bar drives
    # them with a synthetic click — see inst/app/www/module-bar.js for why
    # nav_select() cannot be used against bslib's nav markup.
    #
    # `active_module` is tracked here rather than read from `input$main_tabs`
    # because that input is produced by the same binding whose selector does not
    # match bslib's markup, so it cannot be relied on either.
    active_module <- shiny::reactiveVal("home")

    # Regime comes from `saved`, not `working`: the bar reflects what is on disk,
    # and `working` mutates on every sidebar keystroke.
    output$module_bar <- shiny::renderUI({
        module_bar_ui(active = active_module(),
                      regime = config_regime(config_state$saved))
    })

    # One observer per module, created once at startup — the bar re-renders on
    # every click, so per-render observers would pile up.
    for (.module_id in module_registry()$id) {
        local({
            mid <- .module_id
            shiny::observeEvent(input[[module_input_id(mid)]], {
                active_module(mid)
                session$sendCustomMessage("module_bar_select", list(value = mid))
            }, ignoreInit = TRUE)
        })
    }

    # Create project module
    mod_create_project_server(
        "create_project",
        pipeline_path_rv     = shiny::reactive(get_pipeline_path()),
        # Duplicate-name check must see EVERY project, including simulated ones and
        # internal clone lanes the selector hides -- otherwise a new project could be
        # created onto an existing results tree.
        existing_projects_rv = shiny::reactive(find_projects(get_pipeline_path(), all = TRUE)),
        on_created = function(new_project) {
            pip      <- get_pipeline_path()
            # Same one-list rule as the initial render, or a project created with
            # Simulation.enabled would vanish from the selector right after it
            # was made.
            projects <- find_selectable_projects(pip)
            choices  <- c(
                "+ New Project" = "__new__",
                setNames(projects, projects)
            )
            shiny::updateSelectInput(session, "project_selector",
                                     choices = choices, selected = new_project)
        }
    )

    # ── Global run lock (prevents concurrent pipeline runs) ───────────────────
    pipeline_running <- shiny::reactiveVal(FALSE)

    # ── Trigger to invalidate project_data after a successful run ─────────────
    project_data_trigger <- shiny::reactiveVal(0L)

    # ── Trigger for SNP-set changes (save in GEA / delete in Maladaptation) ───
    # Kept separate from project_data_trigger to avoid reloading heavy p-value tables.
    snp_sets_trigger <- shiny::reactiveVal(0L)

    # ── Reactive project data bundle ───────────────────────────────────────────
    # Depends on project selector AND project_data_trigger (incremented post-run)
    project_data <- shiny::reactive({
        shiny::req(input$project_selector, input$project_selector != "__new__")
        project_data_trigger()  # declare dependency
        proj <- input$project_selector
        make_project_data(proj)
    }) |> shiny::bindCache(input$project_selector, project_data_trigger(),
                           cache = "session")

    # ── Config sidebar servers (one per tab) ───────────────────────────────────
    mod_config_sidebar_server("config_home",          config_state, "home")
    mod_config_sidebar_server("config_processing",    config_state, "processing")
    mod_config_sidebar_server("config_prestructure",  config_state, "prestructure")
    mod_config_sidebar_server("config_structure",     config_state, "structure")
    mod_config_sidebar_server("config_climate",       config_state, "climate")
    mod_config_sidebar_server("config_traits",        config_state, "traits")
    mod_config_sidebar_server("config_pregea",        config_state, "pregea")
    mod_config_sidebar_server("config_gea",           config_state, "gea")
    mod_config_sidebar_server("config_gwas",          config_state, "gwas")
    mod_config_sidebar_server("config_gea_x_gwas",    config_state, "gea_x_gwas")
    mod_config_sidebar_server("config_maladaptation", config_state, "maladaptation",
                              snp_sets_trigger = snp_sets_trigger)

    # ── Save project files button (Home tab) ──────────────────────────────────
    shiny::observeEvent(input$save_project_files, {
        shiny::req(config_state$project)
        tryCatch({
            write_project_config(config_state$working, config_state$project)
            config_state$saved <- config_state$working
            shiny::showNotification("Project files saved.", type = "message", duration = 3)
        }, error = function(e) {
            shiny::showNotification(paste("Save failed:", conditionMessage(e)),
                                    type = "error", duration = 8)
        })
    })

    # ── Pipeline runner servers (one per tab / mode) ───────────────────────────
    mod_pipeline_runner_server("runner_processing",    config_state, pipeline_running,
                                project_data_trigger, tab_name = "processing")
    mod_pipeline_runner_server("runner_prestructure",  config_state, pipeline_running,
                                project_data_trigger, tab_name = "prestructure")
    mod_pipeline_runner_server("runner_structure",     config_state, pipeline_running,
                                project_data_trigger, tab_name = "structure")
    mod_pipeline_runner_server("runner_climate",       config_state, pipeline_running,
                                project_data_trigger, tab_name = "climate")
    mod_pipeline_runner_server("runner_traits",        config_state, pipeline_running,
                                project_data_trigger, tab_name = "traits")
    mod_pipeline_runner_server("runner_pregea",        config_state, pipeline_running,
                                project_data_trigger, tab_name = "pregea")
    mod_pipeline_runner_server("runner_gea",           config_state, pipeline_running,
                                project_data_trigger, tab_name = "gea")
    mod_pipeline_runner_server("runner_gwas",          config_state, pipeline_running,
                                project_data_trigger, tab_name = "gwas")
    mod_pipeline_runner_server("runner_gea_x_gwas",    config_state, pipeline_running,
                                project_data_trigger, tab_name = "gea_x_gwas")
    mod_pipeline_runner_server("runner_maladaptation", config_state, pipeline_running,
                                project_data_trigger, tab_name = "maladaptation",
                                snp_sets_trigger = snp_sets_trigger)

    # ── Module servers ─────────────────────────────────────────────────────────
    mod_home_server("home",                    project_data = project_data)
    mod_processing_server("processing",        project_data = project_data)
    mod_prestructure_server("prestructure",    project_data = project_data)
    mod_structure_server("structure",          project_data = project_data)
    mod_climate_server("climate",              project_data = project_data)
    mod_traits_server("traits",                project_data = project_data)
    mod_pregea_server("pregea",                project_data = project_data)
    mod_gea_server("gea",                      project_data = project_data, run_trigger = project_data_trigger,
                   snp_sets_trigger = snp_sets_trigger, config_state = config_state)
    mod_gwas_server("gwas",                    project_data = project_data, run_trigger = project_data_trigger,
                    config_state = config_state)
    mod_gea_x_gwas_server("gea_x_gwas",        project_data = project_data, run_trigger = project_data_trigger)
    mod_maladaptation_server("maladaptation",  project_data = project_data,
                             snp_sets_trigger = snp_sets_trigger)

    # ── Dashboard-layout outputs ──────────────────────────────────────────────
    # A second moduleServer on the SAME id as its module above, so it shares the
    # namespace and can read that server's inputs. These ADD outputs the original
    # module servers never defined (the KPI strips, PreStructure's Q-matrix
    # table, GEA's suspend override) rather than replacing anything -- see
    # R/mod_dashboard_outputs.R. Order matters for GEA: outputOptions() there
    # needs the outputs already registered by mod_gea_server().
    mod_prestructure_dashboard_outputs("prestructure", project_data)
    mod_structure_dashboard_outputs("structure",       project_data)
    mod_climate_dashboard_outputs("climate",           project_data)
    mod_traits_dashboard_outputs("traits",             project_data)
    mod_pregea_dashboard_outputs("pregea",             project_data)
    mod_gea_dashboard_outputs("gea",                   project_data)
    mod_gwas_dashboard_outputs("gwas",                 project_data)
    mod_geaxgwas_dashboard_outputs("gea_x_gwas",       project_data)
}
