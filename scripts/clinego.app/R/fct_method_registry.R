#' GEA method registry bridge (Python -> R)
#'
#' The GEA method registry (engine, capability flags, per-method hyperparameter
#' specs) lives in workflow/methods/gea.py — the single source of truth shared
#' with the Snakemake pipeline (see workflow/rules/common.smk). Rather than
#' duplicating that spec in R (which would drift the moment a method's params
#' change), this shells out to python3 against the mounted pipeline directory
#' and reads it back as JSON. This is what makes "add a method to gea.py and
#' the Shiny sidebar picks it up with no code change" actually true.
#'
#' Memoised per pipeline_path + gea.py mtime so a pipeline edit (e.g. during
#' development, editing gea.py on the mounted volume) is picked up without an
#' app restart, without re-shelling out on every request.
#' @noRd

.gea_registry_cache <- new.env(parent = emptyenv())

#' Read the GEA_METHODS registry from workflow/methods/gea.py as a nested R list.
#' @return named list: method -> list(engine=, supports_phenotypes=, multivariate=,
#'   pseudo_trait=, supports_wza=, params=list(param_name -> list(type=,default=,...)))
#'   Falls back to a minimal name-only registry (empty params) if python3 or the
#'   import fails, with a one-time warning notification.
#' @noRd
gea_method_registry <- function(pipeline_path = get_pipeline_path()) {
    gea_py <- file.path(pipeline_path, "workflow", "methods", "gea.py")
    mtime  <- tryCatch(as.character(file.info(gea_py)$mtime), error = function(e) NA_character_)
    cache_key <- paste(pipeline_path, mtime)

    cached <- .gea_registry_cache[[cache_key]]
    if (!is.null(cached)) return(cached)

    workflow_dir <- file.path(pipeline_path, "workflow")
    py_code <- paste(
        "import sys, json",
        sprintf("sys.path.insert(0, %s)", shQuote(workflow_dir)),
        "from methods.gea import GEA_METHODS",
        "print(json.dumps(GEA_METHODS))",
        sep = "; "
    )

    result <- tryCatch({
        out <- system2("python3", c("-c", shQuote(py_code)), stdout = TRUE, stderr = TRUE)
        status <- attr(out, "status")
        if (!is.null(status) && status != 0) stop(paste(out, collapse = "\n"))
        jsonlite::fromJSON(paste(out, collapse = ""), simplifyVector = FALSE)
    }, error = function(e) {
        warning("gea_method_registry(): failed to read workflow/methods/gea.py via python3: ",
                conditionMessage(e))
        NULL
    })

    if (is.null(result)) {
        # Degraded-mode fallback: name-only list, no per-method params. The sidebar
        # still functions (method picker works), it just can't show hyperparameters.
        result <- stats::setNames(
            lapply(c("EMMAX", "LFMM", "GLM", "MLM", "CMLM", "ECMLM", "SUPER",
                     "MLMM", "FarmCPU", "BLINK", "RDA"),
                   function(m) list(params = list())),
            c("EMMAX", "LFMM", "GLM", "MLM", "CMLM", "ECMLM", "SUPER",
              "MLMM", "FarmCPU", "BLINK", "RDA")
        )
    }

    .gea_registry_cache[[cache_key]] <- result
    result
}

#' Read the GWAS_METHODS registry (workflow/methods/gwas.py — GEA_METHODS
#' filtered to supports_phenotypes=True; RDA is excluded there since
#' multivariate methods aren't wired into phenotype/GWAS mode).
#' @noRd
gwas_method_registry <- function(pipeline_path = get_pipeline_path()) {
    gwas_py <- file.path(pipeline_path, "workflow", "methods", "gwas.py")
    mtime   <- tryCatch(as.character(file.info(gwas_py)$mtime), error = function(e) NA_character_)
    cache_key <- paste("gwas", pipeline_path, mtime)

    cached <- .gea_registry_cache[[cache_key]]
    if (!is.null(cached)) return(cached)

    workflow_dir <- file.path(pipeline_path, "workflow")
    py_code <- paste(
        "import sys, json",
        sprintf("sys.path.insert(0, %s)", shQuote(workflow_dir)),
        "from methods.gwas import GWAS_METHODS",
        "print(json.dumps(GWAS_METHODS))",
        sep = "; "
    )

    result <- tryCatch({
        out <- system2("python3", c("-c", shQuote(py_code)), stdout = TRUE, stderr = TRUE)
        status <- attr(out, "status")
        if (!is.null(status) && status != 0) stop(paste(out, collapse = "\n"))
        jsonlite::fromJSON(paste(out, collapse = ""), simplifyVector = FALSE)
    }, error = function(e) {
        warning("gwas_method_registry(): failed to read workflow/methods/gwas.py via python3: ",
                conditionMessage(e))
        NULL
    })

    if (is.null(result)) {
        result <- stats::setNames(
            lapply(c("EMMAX", "GLM", "MLM", "CMLM", "ECMLM", "SUPER", "MLMM", "FarmCPU", "BLINK"),
                   function(m) list(params = list())),
            c("EMMAX", "GLM", "MLM", "CMLM", "ECMLM", "SUPER", "MLMM", "FarmCPU", "BLINK")
        )
    }

    .gea_registry_cache[[cache_key]] <- result
    result
}

#' Registry default value for one param, sentinel-resolved against k_best.
#' @noRd
gea_param_default <- function(spec, k_best = NULL) {
    d <- spec$default
    if (identical(d, "@k_best")) return(k_best %||% 3)
    d
}

#' method -> list(param_name -> resolved default) for every method in the registry.
#' Used both by the R-side row renderer and injected into the row-editor JS
#' (as PARAM_SPECS/PARAM_DEFAULTS) so JS-created rows and R-created rows use
#' identical defaults.
#' @noRd
gea_method_param_defaults <- function(registry, k_best = NULL) {
    stats::setNames(lapply(registry, function(cfg) {
        params <- cfg$params %||% list()
        stats::setNames(lapply(params, gea_param_default, k_best = k_best), names(params))
    }), names(registry))
}

#' Merge registry defaults for `method` with any user-supplied param overrides
#' from a config row (m$params, if present). Unknown param names are dropped
#' silently here (the pipeline-side resolve_method_params() in common.smk is
#' the authoritative validator; this is just for display).
#' @noRd
resolve_row_params <- function(method, row_params, registry, k_best = NULL) {
    spec_dict <- registry[[method]]$params %||% list()
    resolved <- stats::setNames(
        lapply(spec_dict, gea_param_default, k_best = k_best),
        names(spec_dict)
    )
    row_params <- row_params %||% list()
    for (name in intersect(names(row_params), names(spec_dict))) {
        resolved[[name]] <- row_params[[name]]
    }
    resolved
}

#' method -> list(adjust=, threshold=, family=) from the registry's
#' significance fields (adjust_default/threshold_default/significance_family,
#' workflow/methods/gea.py). Falls back to bonf/0.05/univariate_pvalue per
#' method for a degraded-mode registry (gea_method_registry() python failure,
#' :52-62 above) where those keys are absent entirely.
#' @noRd
gea_method_significance_defaults <- function(registry) {
    stats::setNames(lapply(registry, function(cfg) list(
        adjust    = cfg$adjust_default      %||% "bonf",
        threshold = as.character(cfg$threshold_default %||% "0.05"),
        family    = cfg$significance_family %||% "univariate_pvalue"
    )), names(registry))
}
