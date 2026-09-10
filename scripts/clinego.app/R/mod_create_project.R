#' Create Project modal module
#'
#' Triggered from the project selector "+ New Project" entry.
#' Guides the user through setting up a new project: project name,
#' input directory, VCF, metadata, and optional GFF file.
#' On creation, copies inst/config_default.yaml, patches user values,
#' and creates the _results/ and _logs/ directories.
#'
#' @noRd

# ── Helper ─────────────────────────────────────────────────────────────────────

#' Convert an absolute path returned by shinyFiles to a path relative to base
#'
#' common.smk prepends '/pipeline/' to input.dir, so config values must be
#' relative to the pipeline root (e.g., "data/" not "/pipeline/data/").
#'
#' @param abs_path character absolute path
#' @param base_path character pipeline root (e.g. "/pipeline")
#' @return relative path with no leading slash
#' @noRd
.make_relative <- function(abs_path, base_path) {
    base <- sub("/$", "", base_path)
    rel  <- sub(paste0("^", base, "/?"), "", abs_path)
    if (!nzchar(rel)) rel <- "."
    rel
}

# ── UI ─────────────────────────────────────────────────────────────────────────

#' UI for the Create Project modal
#'
#' Called from app_server.R via shiny::showModal().
#'
#' @param id module namespace id
#' @noRd
mod_create_project_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::modalDialog(
        title = htmltools::tagList(
            bsicons::bs_icon("folder-plus"), " Create New Project"
        ),
        size  = "m",
        easyClose = FALSE,

        # ── Project name ───────────────────────────────────────────────────────
        shiny::textInput(
            ns("project_name"),
            label       = "Project name",
            placeholder = "MYPROJECT",
            width       = "100%"
        ),
        htmltools::tags$small(
            class = "text-muted",
            "Letters, numbers, underscores only. Must start with a letter."
        ),

        htmltools::tags$hr(),

        # ── Input directory ────────────────────────────────────────────────────
        htmltools::tags$label("Input directory", class = "form-label fw-semibold"),
        htmltools::div(
            class = "input-group mb-1",
            shiny::textInput(
                ns("input_dir"),
                label       = NULL,
                placeholder = "data/",
                width       = "100%"
            ),
            shinyFiles::shinyDirButton(
                ns("input_dir_btn"),
                label = "Browse",
                title = "Select input directory",
                class = "btn btn-outline-secondary"
            )
        ),
        htmltools::tags$small(
            class = "text-muted d-block mb-3",
            "Directory inside /pipeline/ containing your input files."
        ),

        # ── VCF file ───────────────────────────────────────────────────────────
        htmltools::tags$label("VCF file", class = "form-label fw-semibold"),
        htmltools::div(
            class = "input-group mb-1",
            shiny::textInput(
                ns("input_vcf"),
                label       = NULL,
                placeholder = "samples.vcf.gz",
                width       = "100%"
            ),
            shinyFiles::shinyFilesButton(
                ns("vcf_btn"),
                label    = "Browse",
                title    = "Select VCF file",
                multiple = FALSE,
                class    = "btn btn-outline-secondary"
            )
        ),
        htmltools::tags$small(
            class = "text-muted d-block mb-3",
            "VCF filename inside the input directory (.vcf or .vcf.gz)."
        ),

        # ── Metadata file ──────────────────────────────────────────────────────
        htmltools::tags$label("Metadata file", class = "form-label fw-semibold"),
        htmltools::div(
            class = "input-group mb-1",
            shiny::textInput(
                ns("input_metadata"),
                label       = NULL,
                placeholder = "samples_metadata.tsv",
                width       = "100%"
            ),
            shinyFiles::shinyFilesButton(
                ns("metadata_btn"),
                label    = "Browse",
                title    = "Select metadata file",
                multiple = FALSE,
                class    = "btn btn-outline-secondary"
            )
        ),
        htmltools::tags$small(
            class = "text-muted d-block mb-3",
            "TSV with columns: site, sample, latitude, longitude (columns 1\u20134)."
        ),

        # \u2500\u2500 Regime \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        # Sits right after the metadata field: the regime is a statement about
        # whether that file's coordinates are meaningful. Fixed at creation.
        htmltools::tags$label("Project regime", class = "form-label fw-semibold"),
        shiny::radioButtons(
            ns("regime"),
            label    = NULL,
            selected = "standard",
            choiceNames = list(
                htmltools::tagList(
                    htmltools::tags$strong("Standard"),
                    htmltools::tags$span(
                        class = "text-muted d-block small",
                        "Geography and climate. Population structure, GEA, and maladaptation."
                    )
                ),
                htmltools::tagList(
                    htmltools::tags$strong("GWAS only"),
                    htmltools::tags$span(
                        class = "text-muted d-block small",
                        "No coordinates. Structure and GWAS only \u2014 skips WorldClim, GEA, and maladaptation."
                    )
                )
            ),
            choiceValues = list("standard", "gwas_only")
        ),
        htmltools::tags$small(
            class = "text-muted d-block mb-3",
            "Set once, at creation \u2014 it decides which modules this project has."
        ),

        # ── GFF file (optional) ────────────────────────────────────────────────
        htmltools::tags$label(
            "GFF file ",
            htmltools::tags$span(class = "text-muted fw-normal", "(optional)"),
            class = "form-label fw-semibold"
        ),
        htmltools::div(
            class = "input-group mb-1",
            shiny::textInput(
                ns("input_gff"),
                label       = NULL,
                placeholder = "genome.gff3",
                width       = "100%"
            ),
            shinyFiles::shinyFilesButton(
                ns("gff_btn"),
                label    = "Browse",
                title    = "Select GFF file",
                multiple = FALSE,
                class    = "btn btn-outline-secondary"
            )
        ),
        htmltools::tags$small(
            class = "text-muted d-block mb-3",
            "GFF3 annotation file for gene finding and GO enrichment. Leave empty to skip."
        ),

        # ── Validation messages ────────────────────────────────────────────────
        shiny::uiOutput(ns("validation_msg")),

        # ── Footer ─────────────────────────────────────────────────────────────
        footer = htmltools::tagList(
            shiny::modalButton("Cancel"),
            shinyjs::disabled(
                shiny::actionButton(
                    ns("create_btn"),
                    label = htmltools::tagList(
                        bsicons::bs_icon("folder-check"), " Create Project"
                    ),
                    class = "btn-primary"
                )
            )
        )
    )
}

# ── Server ─────────────────────────────────────────────────────────────────────

#' Server for the Create Project modal
#'
#' @param id module namespace id
#' @param pipeline_path_rv reactive returning the pipeline root path string
#' @param existing_projects_rv reactive returning character vector of existing project names
#' @param on_created function(project_name) called after successful project creation
#' @noRd
mod_create_project_server <- function(id, pipeline_path_rv,
                                      existing_projects_rv, on_created) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── shinyFiles roots ───────────────────────────────────────────────────
        roots <- shiny::reactive({
            c(pipeline = pipeline_path_rv())
        })

        # defaultPath for file pickers: input_dir when set, else empty
        file_default_path <- shiny::reactive({
            sub("/$", "", trimws(input$input_dir %||% ""))
        })

        # Register dir picker (always from pipeline root)
        shiny::observe({
            shinyFiles::shinyDirChoose(input, "input_dir_btn",
                                       roots = roots(), session = session)
        })

        # Register file pickers — re-run when defaultPath changes to update
        # the starting directory (browser opens at input_dir if set)
        shiny::observe({
            r   <- roots()
            def <- file_default_path()
            shinyFiles::shinyFileChoose(input, "vcf_btn",
                roots = r, session = session, defaultPath = def)
            shinyFiles::shinyFileChoose(input, "metadata_btn",
                roots = r, session = session, defaultPath = def)
            shinyFiles::shinyFileChoose(input, "gff_btn",
                roots = r, session = session, defaultPath = def)
        })

        # ── shinyFiles → textInput observers ──────────────────────────────────
        shiny::observeEvent(input$input_dir_btn, {
            req(!is.integer(input$input_dir_btn))
            d <- shinyFiles::parseDirPath(roots(), input$input_dir_btn)
            if (length(d) > 0 && nzchar(d)) {
                rel <- .make_relative(d, pipeline_path_rv())
                # Ensure trailing slash for directories
                if (!grepl("/$", rel)) rel <- paste0(rel, "/")
                shiny::updateTextInput(session, "input_dir", value = rel)
            }
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$vcf_btn, {
            req(!is.integer(input$vcf_btn))
            f <- shinyFiles::parseFilePaths(roots(),input$vcf_btn)
            if (nrow(f) > 0) {
                # VCF is stored as filename only (basename), not a full path
                shiny::updateTextInput(session, "input_vcf", value = basename(f$datapath[1]))
            }
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$metadata_btn, {
            req(!is.integer(input$metadata_btn))
            f <- shinyFiles::parseFilePaths(roots(),input$metadata_btn)
            if (nrow(f) > 0) {
                shiny::updateTextInput(session, "input_metadata", value = basename(f$datapath[1]))
            }
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$gff_btn, {
            req(!is.integer(input$gff_btn))
            f <- shinyFiles::parseFilePaths(roots(),input$gff_btn)
            if (nrow(f) > 0) {
                shiny::updateTextInput(session, "input_gff", value = basename(f$datapath[1]))
            }
        }, ignoreInit = TRUE)

        # ── Validation ─────────────────────────────────────────────────────────
        validation_errors <- shiny::reactive({
            errors <- character(0)

            name <- trimws(input$project_name %||% "")
            if (!nzchar(name)) {
                errors <- c(errors, "Project name is required.")
            } else if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", name)) {
                errors <- c(errors, "Project name must start with a letter and contain only letters, numbers, and underscores.")
            } else if (tolower(name) %in% tolower(existing_projects_rv())) {
                errors <- c(errors, paste0("A project named '", name, "' already exists."))
            }

            if (!nzchar(trimws(input$input_dir %||% "")))
                errors <- c(errors, "Input directory is required.")
            if (!nzchar(trimws(input$input_vcf %||% "")))
                errors <- c(errors, "VCF file is required.")
            if (!nzchar(trimws(input$input_metadata %||% "")))
                errors <- c(errors, "Metadata file is required.")

            errors
        })

        # Enable / disable Create button based on validation
        shiny::observe({
            shinyjs::toggleState(
                id        = "create_btn",
                condition = length(validation_errors()) == 0
            )
        })

        # Show validation messages
        output$validation_msg <- shiny::renderUI({
            errs <- validation_errors()
            if (length(errs) == 0) return(NULL)
            htmltools::div(
                class = "alert alert-warning py-2 mt-2",
                lapply(errs, function(e) {
                    htmltools::div(bsicons::bs_icon("exclamation-triangle"), " ", e)
                })
            )
        })

        # ── Create project handler ─────────────────────────────────────────────
        shiny::observeEvent(input$create_btn, {
            req(length(validation_errors()) == 0)

            pip  <- pipeline_path_rv()
            name <- trimws(input$project_name)
            dir_val      <- trimws(input$input_dir)
            vcf_val      <- trimws(input$input_vcf)
            metadata_val <- trimws(input$input_metadata)
            gff_val      <- trimws(input$input_gff %||% "")

            tryCatch({
                # 1. Copy template config
                template <- app_sys("config_default.yaml")
                dest     <- file.path(pip, paste0("config_", name, ".yaml"))
                file.copy(template, dest, overwrite = FALSE)

                # 2. Read, patch, write
                cfg <- yaml::read_yaml(dest)
                cfg$project_name    <- name
                cfg$Input$dir       <- dir_val
                cfg$Input$vcf       <- vcf_val
                cfg$Input$metadata  <- metadata_val
                if (nzchar(gff_val)) {
                    cfg$Input$gff <- gff_val
                }
                # Regime is declared here and nowhere else; Climate.enabled is the
                # derived value the pipeline reads (see fct_config_writer.R).
                regime_val <- input$regime %||% "standard"
                cfg$Regime$mode     <- regime_val
                cfg$Climate$enabled <- !identical(regime_val, "gwas_only")
                yaml::write_yaml(cfg, dest)

                # 3. Create result and log directories
                dir.create(file.path(pip, paste0(name, "_results")),
                           showWarnings = FALSE, recursive = TRUE)
                dir.create(file.path(pip, paste0(name, "_logs")),
                           showWarnings = FALSE, recursive = TRUE)

                # 4. Close modal and notify parent
                shiny::removeModal()
                on_created(name)

            }, error = function(e) {
                shiny::showNotification(
                    paste("Failed to create project:", conditionMessage(e)),
                    type = "error", duration = 10
                )
            })
        })
    })
}
