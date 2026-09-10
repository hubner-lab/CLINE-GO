#' Config writer — serialize in-app config state back to YAML
#'
#' Provides utilities for writing the working config (collected from sidebar
#' inputs) back to the project YAML file. Used by the pipeline runner when
#' the user clicks Run: working config is persisted before Snakemake starts.
#'
#' Note: R's yaml library does not preserve comments. The written YAML will
#' have the same structure and values but no inline comments.

#' Get the config file path for a project
#'
#' Returns config_{project}.yaml if it exists, otherwise config.yaml.
#' @noRd
config_file_path <- function(project, pipeline_path = get_pipeline_path()) {
    primary <- file.path(pipeline_path, paste0("config_", project, ".yaml"))
    if (file.exists(primary)) primary else file.path(pipeline_path, "config.yaml")
}

#' Collect all sidebar input values for a given set of schema entries
#'
#' Reads each input by its namespaced ID, coerces to the correct type, and
#' sets the value in a copy of the base config using config_set_by_path().
#'
#' @param input Shiny input object (from moduleServer context)
#' @param ns namespace function for the sidebar module
#' @param schema_entries list of schema entries for this tab
#' @param base_config list — starting config (typically config_state$working)
#' @return modified config list
#' @noRd
collect_tab_config <- function(input, ns, schema_entries, base_config) {
    cfg <- base_config
    for (e in schema_entries) {
        if (e$type == "method_table") next  # read-only
        input_id <- gsub("\\.", "_", e$key)
        raw <- input[[input_id]]
        if (is.null(raw)) next
        val <- input_to_config_value(raw, e$type)
        cfg <- config_set_by_path(cfg, e$key, val)
    }
    cfg
}

#' Write a config list to YAML atomically
#'
#' Writes to a temp file in the same directory, then renames into place.
#' This prevents partial writes from corrupting the config file.
#'
#' @param config named list (full project config)
#' @param project project name
#' @param pipeline_path pipeline root path
#' @return invisible path to the written file, or NULL on error
#' @noRd
write_project_config <- function(config, project,
                                  pipeline_path = get_pipeline_path()) {
    dest <- config_file_path(project, pipeline_path)

    # Materialize the regime and the value the pipeline actually reads. Every write
    # goes through here (Save button and the pipeline runner both hand over
    # config_state$working), so this is the one place the derivation belongs.
    # Regime.mode is written explicitly, not just the derived flag: that pins the
    # migration inference to a single run, so a later hand-edit of Climate.enabled
    # cannot silently re-classify the project.
    regime <- config_regime(config)
    config <- config_set_by_path(config, "Regime.mode",     regime)
    config <- config_set_by_path(config, "Climate.enabled", !identical(regime, "gwas_only"))

    # Prepare YAML handlers for special R values
    prepared <- prepare_config_for_yaml(config)

    # Write to temp file in same directory (atomic rename requires same fs)
    tmp <- tempfile(tmpdir = dirname(dest), fileext = ".yaml")
    result <- tryCatch({
        yaml::write_yaml(prepared, tmp)
        file.rename(tmp, dest)
        invisible(dest)
    }, error = function(e) {
        if (file.exists(tmp)) file.remove(tmp)
        warning("write_project_config: failed to write config: ", conditionMessage(e))
        NULL
    })
    result
}

#' Recursively prepare a config list for yaml::write_yaml()
#'
#' - Removes NULL entries (omit from YAML)
#' - Converts NA to NULL
#' - Ensures association.configs list structure is preserved
#' @noRd
prepare_config_for_yaml <- function(x) {
    if (is.list(x) && !is.null(names(x))) {
        # Named list: recurse and drop NULLs
        result <- lapply(x, prepare_config_for_yaml)
        result[vapply(result, is.null, logical(1))] <- NULL
        result
    } else if (is.list(x)) {
        # Unnamed list (e.g. association.configs): recurse but keep all entries
        lapply(x, prepare_config_for_yaml)
    } else if (length(x) == 1 && is.na(x)) {
        NULL
    } else {
        x
    }
}
