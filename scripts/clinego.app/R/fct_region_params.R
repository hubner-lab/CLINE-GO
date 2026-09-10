#' Per-region exploratory parameter persistence
#'
#' Stores parameters used for on-demand computations (hap scan/viz, region_distance)
#' in a project-local JSON file, separate from the pipeline config YAML.
#' This is Shiny-side "session memory" — not a pipeline parameter.
#'
#' File: {PROJECT}_results/_intermediate/region_params.json
#' Structure:
#'   global$module$key          — e.g., global$gea$region_distance
#'   regions$module$region_id$type — e.g., regions$gea$5_7100000-19350376$hap_scan
#'
#' @noRd

.empty_region_params <- function() {
    list(global = list(), regions = list())
}

#' Path to the region params JSON file for a project
#' @noRd
region_params_path <- function(project) {
    file.path(project_base(project), "_intermediate", "region_params.json")
}

#' Read per-region params from disk. Returns empty skeleton on missing/corrupt file.
#' @noRd
read_region_params <- function(project) {
    path <- region_params_path(project)
    if (!file_ok(path)) return(.empty_region_params())
    tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(e) .empty_region_params()
    )
}

#' Write per-region params atomically (tempfile + rename).
#' @noRd
save_region_params <- function(project, params_list) {
    path <- region_params_path(project)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".json")
    tryCatch({
        jsonlite::write_json(params_list, tmp, auto_unbox = TRUE, pretty = FALSE)
        file.rename(tmp, path)
    }, error = function(e) {
        if (file.exists(tmp)) unlink(tmp)
        message("WARNING: Could not save region params: ", e$message)
    })
    invisible(params_list)
}

#' Get a stored computation's params for a specific region.
#' @param params Result of read_region_params()
#' @param module Module constant string ("gea", "gwas", "gea_x_gwas")
#' @param region_id Region id string (e.g., "5_7100000-19350376")
#' @param type Computation type ("hap_scan" or "hap_viz")
#' @return Named list of params, or NULL if not found
#' @noRd
get_region_param <- function(params, module, region_id, type) {
    if (is.null(params) || is.null(module) || is.null(region_id) || is.null(type))
        return(NULL)
    tryCatch(
        params$regions[[module]][[region_id]][[type]],
        error = function(e) NULL
    )
}

#' Store a computation's params for a specific region.
#' Returns the modified params list (does not write to disk — caller must call save_region_params).
#' @noRd
set_region_param <- function(params, module, region_id, type, values) {
    if (is.null(params$regions)) params$regions <- list()
    if (is.null(params$regions[[module]])) params$regions[[module]] <- list()
    if (is.null(params$regions[[module]][[region_id]]))
        params$regions[[module]][[region_id]] <- list()
    params$regions[[module]][[region_id]][[type]] <- values
    params
}

#' Get a global (non-region-specific) param.
#' @param params Result of read_region_params()
#' @param module Module constant string
#' @param key Param key (e.g., "region_distance")
#' @return Scalar value or NULL
#' @noRd
get_global_param <- function(params, module, key) {
    if (is.null(params)) return(NULL)
    tryCatch(params$global[[module]][[key]], error = function(e) NULL)
}

#' Set a global (non-region-specific) param.
#' Returns the modified params list.
#' @noRd
set_global_param <- function(params, module, key, value) {
    if (is.null(params$global)) params$global <- list()
    if (is.null(params$global[[module]])) params$global[[module]] <- list()
    params$global[[module]][[key]] <- value
    params
}
