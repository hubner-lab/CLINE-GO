# Tell data.table that this package uses its syntax (avoids cedta() errors)
.datatable.aware <- TRUE

#' Keep only canonical merged per-method pvalue/wza files from a glob result.
#'
#' GWAS DROP-mode (Path B) writes per-trait files ({trait}_pvalues_K{k}.tsv)
#' directly into the same methods/{METHOD}/ directory as the merged
#' {METHOD}_pvalues_K{k}.tsv. The canonical merged file always starts with the
#' directory (method) name followed by the regime token ("_pvalues_K" or "_wza_K").
#' Per-trait files never start with the method name, so this filter excludes them.
#'
#' @param files Character vector of glob results (already filtered for sig_snps/wza/block).
#' @param regime "snp" (default) or "wza".
#' @noRd
keep_merged_method_files <- function(files, regime = "snp") {
    if (length(files) == 0) return(files)
    token <- if (regime == "wza") "_wza_K" else "_pvalues_K"
    keep <- vapply(files, function(f) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        mi <- which(parts == "methods")
        if (length(mi) == 0) return(FALSE)
        startsWith(basename(f), paste0(parts[mi + 1L], token))
    }, logical(1))
    files[keep]
}

#' Check if a file path exists and is non-empty
#' @noRd
file_ok <- function(path) {
    !is.null(path) && length(path) == 1 && !is.na(path) &&
        nzchar(path) && file.exists(path) && file.size(path) > 0
}

#' Safely read a JSON file, returning NULL on error
#' @noRd
safe_read_json <- function(path) {
    if (!file_ok(path)) return(NULL)
    tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

#' Resolve K_best: from config, then fall back to max available K
#' @noRd
resolve_k_best <- function(project, config = NULL) {
    k_cfg <- if (!is.null(config)) config_k_best(config) else NA_integer_
    if (!is.na(k_cfg) && k_cfg > 0L) return(k_cfg)
    ks <- find_k_values(project)
    if (length(ks) == 0) return(NA_integer_)
    # Flag the fallback: config k_best is absent/invalid so we guess max available K.
    # This is the condition that can pin a stale K (e.g. K6 while outputs are K5) — keep
    # the fallback (needed pre-sNMF) but make it visible in logs rather than silent.
    k_fallback <- max(ks)
    warning(sprintf(
        "resolve_k_best: config sNMF.k_best absent/invalid for project '%s'; falling back to max available K = %s",
        project, k_fallback))
    k_fallback
}

#' Resolve adjust string for a method from ASSOC_CONFIGS entries
#' Returns "bonf_0.05" style string or NULL
#' @noRd
resolve_adjust <- function(config, method, module = "GEA") {
    # Try module-specific configs first (e.g. GWAS.configs)
    configs <- config_get(config, module, "configs", default = list())
    # Fall back to GEA.configs (GWAS inherits from GEA)
    if (length(configs) == 0) configs <- config_assoc_configs(config)
    for (cfg in configs) {
        if (identical(cfg$method, method)) {
            return(paste0(cfg$adjust, "_", cfg$threshold))
        }
    }
    NULL
}

#' Create a reactive project data bundle for use across modules
#' @noRd
make_project_data <- function(project, pipeline_path = get_pipeline_path()) {
    config <- read_project_config(project, pipeline_path)
    k_best <- resolve_k_best(project, config)
    list(
        name    = project,
        path    = project_base(project, pipeline_path),
        config  = config,
        k_best  = k_best,
        status  = check_module_status(project)
    )
}

#' Map Shiny tab name to Snakemake mode string
#' Returns NULL for haplotype (handled by two separate runner instances)
#' @noRd
tab_to_mode <- function(tab) {
    switch(tab,
        processing    = "processing",
        prestructure  = "prestructure",
        climate       = "climate",
        traits        = "traits",
        pregea        = "pregea",
        structure     = "structure",
        gea           = "gea",
        gwas          = "gwas",
        gea_x_gwas    = "gea_x_gwas",
        maladaptation = "maladaptation",
        NULL
    )
}

#' Safe DT rendering with "no data" message fallback
#' @noRd
safe_datatable <- function(dt, ...) {
    if (is.null(dt) || nrow(dt) == 0) {
        dt <- data.frame(Message = "No data available")
    }
    DT::datatable(dt, ...)
}

#' Format p-value columns in a DT with signif formatting
#' @noRd
dt_format_pvals <- function(tbl, cols = c("pvalue", "p_adjust")) {
    cols_present <- intersect(cols, names(tbl))
    if (length(cols_present) == 0) return(tbl)
    DT::formatSignif(tbl, columns = cols_present, digits = 3)
}

#' Default interactive SNP clumping distance (bp) for the GEA/GWAS filter bar.
#'
#' Clumping distance is a purely interactive/exploratory parameter (not read from
#' the pipeline config) so the filter-bar control is always the single source of
#' truth for region merging/padding — no separate config-driven default to diverge
#' from. User-chosen values persist per-project via region_params.json
#' (get_global_param()/set_global_param()), independent of this constant.
#' @noRd
DEFAULT_CLUMPING_DISTANCE <- 100000L

#' Reactive wrapper that only propagates when the VALUE actually changes.
#'
#' Shiny's `reactive()` has no value memo: it re-invalidates every downstream
#' consumer whenever an upstream dependency fires, even when the recomputed
#' value is identical. That is normally harmless, but it is expensive for the
#' per-method Manhattan, whose driving inputs resolve in stages on a cold tab:
#' `input$method_tab` and `input$per_method_trait` live in `renderUI` outputs, so
#' each reactive first yields its `%||%` fallback ("EMMAX", first trait) and then
#' yields the SAME string again a beat later when the real input registers. Each
#' of those identical values forced a full `Plotly.newPlot`, which re-rasterizes
#' the large base64 background PNG -- measured as 3 renders in ~280 ms, seen as a
#' flicker.
#'
#' `reactiveVal` already refuses to invalidate on an `identical()` write, so
#' routing the value through one collapses the burst to a single render while
#' leaving genuine changes untouched (one extra flush of latency, no more).
#'
#' @param r a reactive expression
#' @return a reactive that fires only on real value changes
#' @noRd
dedupe_reactive <- function(r) {
    # Seed from the current value rather than starting NULL: reactiveVal() would
    # otherwise report NULL for the one flush before the observer first runs, and
    # a consumer reading in that window would see NULL where the plain reactive
    # gave it the `%||%` fallback. Observer priority happens to cover this today
    # -- seeding removes the dependency on that ordering.
    rv <- shiny::reactiveVal(shiny::isolate(r()))
    shiny::observe(rv(r()))
    shiny::reactive(rv())
}
