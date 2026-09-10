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

#' Get the pipeline root path (where {PROJECT}_results/ directories live)
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
