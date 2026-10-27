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

#' Stand-in clumping distance (bp) for the GEA/GWAS filter bar.
#'
#' Used only when default_clumping_distance() cannot read the project's own LD-decay
#' distance: no LD-decay table yet (mode=structure not run), a table that predates the
#' PLINK estimator, or a rule the curve never reaches.
#' @noRd
DEFAULT_CLUMPING_DISTANCE <- 100000L

#' The filter bar's starting clumping distance (bp): the project's own.
#'
#' Reads `snp_clumping_distance` of `source` as the pipeline's resolve_region_params()
#' does: GWAS falls back to GEA, and the legacy `region_distance` / `combine_gap` names
#' stand in for an unset key. A bp value is used as given. An auto_* value, or none at
#' all (the pipeline's default since 2026-10-26), is the genome-wide LD-decay distance
#' under the project's `clumping_ld_rule` and `ld_decay_group`, read with the pipeline's
#' own resolve_clumping_distance() (scripts/R/lib/region_distance.R). The bar holds one
#' number, so auto_per_chromosome is read genome-wide, as in declared SNP sets.
#' The filter bar stays the source of truth once moved: a saved value
#' (region_params.json) wins over this one.
#'
#' @param cfg      Project config (pd$config).
#' @param source   MOD_GEA or MOD_GWAS.
#' @param ld_table Path of the LD-decay table, ld_decay_table_path(project).
#' @return Integer bp; DEFAULT_CLUMPING_DISTANCE when the LD distance cannot be read
#'   or the spec is neither bp nor auto_*.
#' @noRd
default_clumping_distance <- function(cfg, source = MOD_GEA, ld_table = "NULL") {
    r <- .resolve_config_clumping(cfg, source, ld_table)
    if (is.null(r$bp)) DEFAULT_CLUMPING_DISTANCE else r$bp
}

#' Why the project's own clumping distance could not be read, or NULL when it was.
#'
#' The pipeline stops on these cases (e.g. an auto_* distance on an unreliable LD-decay
#' curve) while the filter bar starts at DEFAULT_CLUMPING_DISTANCE; this is the reason to
#' show beside it, so the stand-in is never silent.
#' @inheritParams default_clumping_distance
#' @return Character message or NULL.
#' @noRd
clumping_distance_issue <- function(cfg, source = MOD_GEA, ld_table = "NULL") {
    .resolve_config_clumping(cfg, source, ld_table)$error
}

# The config's clumping distance of `source` in bp (`bp`), or the reason it could not be
# read (`error`); both NULL for a spec that is neither bp nor auto_*.
.resolve_config_clumping <- function(cfg, source, ld_table) {
    group_keys <- function(group) {
        d <- config_get(cfg, group, default = list())
        if (!is.list(d)) d <- list()
        d$snp_clumping_distance <- d$snp_clumping_distance %||% d$region_distance %||% d$combine_gap
        d$ld_decay_group <- d$ld_decay_group %||% d$region_ld_decay_group
        d
    }
    own <- group_keys(if (identical(source, MOD_GWAS)) "GWAS" else "GEA")
    gea <- group_keys("GEA")
    pick <- function(key, default) own[[key]] %||% gea[[key]] %||% default

    spec <- pick("snp_clumping_distance", "auto_genome_wide")
    bp <- suppressWarnings(as.numeric(spec))
    if (length(bp) == 1 && !is.na(bp)) return(list(bp = as.integer(bp), error = NULL))
    if (!tolower(trimws(as.character(spec)))[1] %in%
            c("auto_genome_wide", "auto_per_chromosome", "auto")) {
        return(list(bp = NULL, error = NULL))
    }
    tryCatch(
        list(bp = suppressMessages(resolve_clumping_distance(
                 "auto_genome_wide", ld_table, rule = pick("clumping_ld_rule", "critical"),
                 group = pick("ld_decay_group", "All"))),
             error = NULL),
        error = function(e) list(bp = NULL, error = conditionMessage(e)))
}

#' A filter bar's clumping distance: one value per project, saved only when moved.
#'
#' Call from a module server; the box must be rendered from the returned reactive,
#' never from its own input. The start is the value saved under `module`/`rp_key`,
#' else default_clumping_distance() for `source`. A move is kept for the session and
#' saved under the project it was made in; nothing else is saved, so opening a tab
#' writes nothing.
#'
#' Every re-render binds a new box, which reports its own starting value, and on a cold
#' app that report can arrive after a project switch. Taken as a move, the previous
#' project's value was saved into the new project's region_params.json and the bar
#' re-rendered in a loop between the two values (measured 2026-10-26, once per-project
#' starting values made them differ). So after the app puts a value in the box, every
#' report is an echo until the box reports that value back.
#'
#' @param input        The module's input.
#' @param project_data Reactive project bundle (make_project_data()).
#' @param input_id     The box's input id (un-namespaced).
#' @param module       region_params.json module key (MOD_GEA, MOD_GWAS, MOD_GEAXGWAS).
#' @param rp_key       Key under `module`.
#' @param source       MOD_GEA or MOD_GWAS: whose config distance is the start.
#' @return Reactive integer bp.
#' @noRd
clumping_distance_state <- function(input, project_data, input_id, module, rp_key, source) {
    start <- shiny::reactive({
        pd <- project_data()
        if (is.null(pd)) return(DEFAULT_CLUMPING_DISTANCE)
        saved <- get_global_param(read_region_params(pd$name), module, rp_key)
        if (!is.null(saved)) as.integer(saved)
        else default_clumping_distance(pd$config, source, ld_decay_table_path(pd$name))
    })
    moved <- shiny::reactiveVal(list())   # this session's moves, by project name
    value <- shiny::reactive({
        pd <- project_data()
        (if (!is.null(pd)) moved()[[pd$name]]) %||% start()
    })
    # The value the app last put in the box that the box has not reported back yet;
    # NULL once it has, or when the box already held it.
    pending <- NULL
    shiny::observe({
        now <- value()
        cur <- shiny::isolate(input[[input_id]])
        held <- !is.null(cur) && !is.na(cur) && identical(as.integer(cur), now)
        pending <<- if (held) NULL else now
    })
    shiny::observeEvent(input[[input_id]], {
        v  <- input[[input_id]]
        pd <- project_data()
        if (is.null(v) || is.na(v) || v < 1L || is.null(pd)) return()
        v <- as.integer(v)
        if (!is.null(pending)) {
            if (identical(v, pending)) pending <<- NULL
            return()
        }
        if (identical(v, shiny::isolate(value()))) return()
        m <- shiny::isolate(moved())
        m[[pd$name]] <- v
        moved(m)
        rp <- read_region_params(pd$name)
        save_region_params(pd$name, set_global_param(rp, module, rp_key, v))
    }, ignoreInit = TRUE)
    value
}

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
