#' @import shiny
#' @importFrom golem get_golem_options
NULL

#' Access files in the inst/ directory
#' @param ... path components passed to file.path
#' @noRd
app_sys <- function(...) {
    # In dev mode the source tree is mounted at /pipeline — use it directly so
    # CSS/static changes take effect immediately without a Docker rebuild.
    src_path <- file.path("/pipeline/scripts/clinego.app/inst", ...)
    if (file.exists(src_path)) return(src_path)
    system.file(..., package = "clinego.app")
}

#' Read golem config values
#' @param value config key to read
#' @param config which config environment
#' @param use_parent whether to use parent env
#' @noRd
get_golem_config <- function(value, config = Sys.getenv("R_CONFIG_ACTIVE", "default"),
                              use_parent = TRUE) {
    config::get(
        value = value,
        config = config,
        file = app_sys("golem-config.yml"),
        use_parent = use_parent
    )
}

#' Get the project root: the folder that holds config_<p>.yaml, <p>_results/,
#' <p>_logs/ and the project's inputs
#'
#' Priority: option clinego.project_root > env CLINEGO_PROJECT_ROOT > the code root
#' (get_pipeline_path()), where code and projects sat together before v1.0. Never set
#' in golem-config.yml: there it would beat the env var `clinego app` passes. Code
#' paths (Snakefile, scripts/, workflow/, inst/) keep get_pipeline_path().
#' @noRd
get_project_root <- function() {
    opt <- getOption("clinego.project_root")
    if (!is.null(opt) && nzchar(opt)) return(opt)
    env <- Sys.getenv("CLINEGO_PROJECT_ROOT", unset = "")
    if (nzchar(env)) return(env)
    get_pipeline_path()
}

#' Get the code root (where the Snakefile, scripts/ and workflow/ live)
#' @noRd
get_pipeline_path <- function() {
    # Priority: runtime option > golem config > env var > /pipeline
    opt <- getOption("clinego.pipeline_path")
    if (!is.null(opt)) return(opt)
    cfg <- tryCatch(get_golem_config("pipeline_path"), error = function(e) NULL)
    if (!is.null(cfg) && nzchar(cfg)) return(cfg)
    env <- Sys.getenv("PIPELINE_PATH", unset = "")
    if (nzchar(env)) return(env)
    "/pipeline"
}
