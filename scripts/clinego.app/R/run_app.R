#' Run the CLINE-GO results viewer
#'
#' @param pipeline_path path to the pipeline root (default: /pipeline)
#' @param ... additional arguments passed to shiny::shinyApp() options
#' @export
run_app <- function(pipeline_path = NULL, ...) {
    if (!is.null(pipeline_path)) {
        options(clinego.pipeline_path = pipeline_path)
    }

    shiny::shinyApp(
        ui     = app_ui,
        server = app_server,
        options = list(...)
    )
}
