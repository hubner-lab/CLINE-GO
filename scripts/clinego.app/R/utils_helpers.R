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

#' A config list as a character vector: a comma-separated string (Climate.predictors,
#' Future.models), a YAML list or a vector; entries trimmed, empty ones dropped.
#' @noRd
parse_csv_list <- function(x) {
    if (is.null(x)) return(character(0))
    v <- trimws(unlist(strsplit(as.character(unlist(x)), ",")))
    v[nzchar(v)]
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

#' What make_project_data() reads from disk, as a cache key part.
#'
#' The modification times of the project's config file and pipeline_summary.tsv:
#' a project bundle cached under an older key is stale (a config saved, a run
#' finished in another session or on the command line), fb53fe.
#' @noRd
project_files_stamp <- function(project, pipeline_path = get_pipeline_path()) {
    if (is.null(project) || length(project) != 1 || !nzchar(project)) return("")
    paths <- c(config_file_path(project, pipeline_path), pipeline_summary_path(project))
    paste(format(file.mtime(paths), "%Y%m%d%H%M%OS6"), collapse = "|")
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
safe_datatable <- function(dt, ..., options = list()) {
    if (is.null(dt) || nrow(dt) == 0) {
        dt <- data.frame(Message = "No data available")
    }
    # a table's own CSV button exports the data, not the 3-digit display (N2)
    if (!is.null(options$buttons)) {
        options$buttons <- lapply(options$buttons, function(b) {
            if (identical(b, "csv")) list(extend = "csv", exportOptions = list(orthogonal = "export"))
            else b
        })
    }
    dt_round_doubles(DT::datatable(dt, ..., options = options), dt)
}

#' Numbers on screen at 3 significant digits (N2), per cell as A3's mock had it: each value
#' of a double column that holds non-integer values shows as
#' String(Number(Number(x).toPrecision(3))), so 0.185926 reads 0.186, 2.5e-09 reads 2.5e-9
#' and 1234.5 reads 1230, whatever else the column holds. The render runs for display only:
#' sorting, the card's CSV download and the table's own CSV button (orthogonal "export")
#' keep full precision. Integer-valued columns (positions, counts, ids) are left alone.
#' @param tbl A DT::datatable() of `df`.
#' @noRd
dt_round_doubles <- function(tbl, df) {
    frac <- names(df)[vapply(df, function(x) {
        is.double(x) && any(is.finite(x) & x != round(x))
    }, logical(1))]
    # by the data's own names: a rownames column shifts the index, colnames = does not
    targets <- stats::na.omit(match(frac, names(tbl$x$data))) - 1L
    if (length(targets) == 0) return(tbl)
    tbl$x$options$columnDefs <- c(tbl$x$options$columnDefs, list(list(
        targets = as.integer(targets),
        render  = DT::JS(
            "function(d, t) {",
            "  if (t !== 'display' || d === null || d === '' || isNaN(Number(d))) return d;",
            "  return String(Number(Number(d).toPrecision(3)));",
            "}"))))
    tbl
}

#' PreGEA's recommendations as the app shows them (D3, cc17b6): Method, Parameter, Value,
#' Rule (a "note" badge carries the row's note in its title) and Evidence as
#' "metric = value" at 3 significant digits. HTML-escaped here: the table renders with
#' escape = FALSE for the badge.
#' @param dt pregea_recommendations.tsv as read (any column may be missing).
#' @noRd
recommendations_display <- function(dt) {
    if (is.null(dt) || nrow(dt) == 0) return(data.frame())
    col <- function(name) {
        v <- if (name %in% names(dt)) as.character(dt[[name]]) else rep("", nrow(dt))
        v[is.na(v)] <- ""
        v
    }
    esc <- function(x) vapply(x, function(v) as.character(htmltools::htmlEscape(v)), "",
                              USE.NAMES = FALSE)
    value <- vapply(col("evidence_value"), function(v) {
        num <- suppressWarnings(as.numeric(v))
        if (is.na(num)) v else format(signif(num, 3))
    }, "", USE.NAMES = FALSE)
    metric <- col("evidence_metric")
    note   <- col("note")
    data.frame(
        Method    = esc(col("method")),
        Parameter = esc(col("param")),
        Value     = esc(col("recommended_value")),
        Rule      = paste0(esc(col("rule")),
                           ifelse(nzchar(note),
                                  paste0(' <span class="badge bg-secondary" title="', esc(note),
                                         '">note</span>'), "")),
        Evidence  = ifelse(!nzchar(metric), "",
                           ifelse(nzchar(value), paste0(esc(metric), " = ", esc(value)),
                                  esc(metric))),
        stringsAsFactors = FALSE
    )
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

#' The region merge's own distance from the config, when it sets one.
#'
#' The pipeline merges significant SNPs into regions at `region_clumping_distance`
#' when the config sets it (GWAS inheriting GEA's, as resolve_region_params() reads
#' it), else at the clumping distance (acff31). Read as default_clumping_distance()
#' reads the clumping distance, an auto_* value genome-wide from the LD-decay table.
#'
#' @inheritParams default_clumping_distance
#' @return list(bp, error): both NULL when the config sets no region distance;
#'   `error` the reason an auto_* distance could not be read.
#' @noRd
config_region_clumping_distance <- function(cfg, source = MOD_GEA, ld_table = "NULL") {
    .resolve_config_clumping(cfg, source, ld_table, key = "region_clumping_distance")
}

# The config's clumping distance of `source` in bp (`bp`), or the reason it could not be
# read (`error`); both NULL for a spec that is neither bp nor auto_*. `key` is
# snp_clumping_distance (unset = the pipeline's auto_genome_wide default) or
# region_clumping_distance (unset, null or "none" = not set: both NULL).
.resolve_config_clumping <- function(cfg, source, ld_table, key = "snp_clumping_distance") {
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

    spec <- pick(key, if (identical(key, "snp_clumping_distance")) "auto_genome_wide")
    if (is.null(spec) || tolower(trimws(as.character(spec)))[1] %in% c("", "none", "null"))
        return(list(bp = NULL, error = NULL))
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

#' A filter-bar control's value: one per project, saved only when the user moves it.
#'
#' Call from a module server; the control must be rendered (or updated) from the
#' returned reactive, never from its own input. The start is `start_fn(pd, saved)`,
#' where `saved` is the value region_params.json holds under `module`/`rp_key` (NULL
#' when none). A move is kept for the session and saved under the project it was made
#' in; nothing else is saved, so opening a tab writes nothing.
#'
#' Every re-render binds a new control, which reports its own starting value, and on a
#' cold app that report can arrive after a project switch. Taken as a move, the
#' previous project's value was saved into the new project's region_params.json and
#' the bar re-rendered in a loop between the two values (measured 2026-10-26 on the
#' clumping distance, once per-project starting values made them differ), and every
#' control wrote its starting value on first show (ac8228). So after the app puts a
#' value in the control, every report is an echo until the control reports that value
#' back.
#'
#' @param input        The module's input.
#' @param project_data Reactive project bundle (make_project_data()).
#' @param input_id     The control's input id (un-namespaced).
#' @param module       region_params.json module key (MOD_GEA, MOD_GWAS, MOD_GEAXGWAS).
#' @param rp_key       Key under `module`.
#' @param start_fn     function(pd, saved): the project's starting value; `pd` is NULL
#'   before a project is loaded.
#' @param coerce       function(v): the control's report as the stored value, or NULL
#'   for a report to ignore (an empty or invalid box).
#' @return Reactive value.
#' @noRd
filter_param_state <- function(input, project_data, input_id, module, rp_key,
                               start_fn, coerce = identity) {
    start <- shiny::reactive({
        pd <- project_data()
        if (is.null(pd)) return(start_fn(NULL, NULL))
        start_fn(pd, get_global_param(read_region_params(pd$name), module, rp_key))
    })
    moved <- shiny::reactiveVal(list())   # this session's moves, by project name
    value <- shiny::reactive({
        pd <- project_data()
        (if (!is.null(pd)) moved()[[pd$name]]) %||% start()
    })
    # The value the app last put in the control that the control has not reported back
    # yet; NULL once it has, or when the control already held it.
    pending <- NULL
    shiny::observe({
        now <- value()
        cur <- shiny::isolate(input[[input_id]])
        held <- !is.null(cur) && identical(coerce(cur), now)
        pending <<- if (held) NULL else now
    })
    shiny::observeEvent(input[[input_id]], {
        pd <- project_data()
        v  <- coerce(input[[input_id]])
        if (is.null(v) || is.null(pd)) return()
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

#' A filter bar's clumping distance (bp), saved only when moved (filter_param_state()).
#'
#' The start is the value saved under `module`/`rp_key`, else default_clumping_distance()
#' for `source`.
#'
#' @inheritParams filter_param_state
#' @param source MOD_GEA or MOD_GWAS: whose config distance is the start.
#' @return Reactive integer bp.
#' @noRd
clumping_distance_state <- function(input, project_data, input_id, module, rp_key, source) {
    filter_param_state(
        input, project_data, input_id, module, rp_key,
        start_fn = function(pd, saved) {
            if (is.null(pd)) return(DEFAULT_CLUMPING_DISTANCE)
            if (!is.null(saved)) as.integer(saved)
            else default_clumping_distance(pd$config, source, ld_decay_table_path(pd$name))
        },
        coerce = function(v) {
            if (is.null(v) || length(v) != 1 || is.na(v) || v < 1L) NULL else as.integer(v)
        })
}

#' The four controls of a GEA/GWAS threshold and filter bar, saved only when moved.
#'
#' The threshold rule (type and value), the combine strategy and the regime switch,
#' each a filter_param_state(). The ids are `prefix` + "threshold_type",
#' "threshold_value", "combine_strategy", "regime"; the region_params.json keys are the
#' same names unless `rp_keys` renames them (GEAxGWAS keeps its historic
#' "gea_strategy" / "gwas_strategy").
#'
#' @inheritParams filter_param_state
#' @param prefix            Input-id prefix ("" or "gea_" / "gwas_").
#' @param source            MOD_GEA or MOD_GWAS: whose config rule is the start.
#' @param registry_defaults Reactive method -> list(adjust, threshold, family).
#' @param default_strategy  Strategy when none is saved.
#' @param rp_keys           Named character: region_params.json key per control.
#' @return list(threshold_type, threshold_value, strategy, regime), each a reactive.
#' @noRd
filter_bar_state <- function(input, project_data, module, source, registry_defaults,
                             prefix = "", default_strategy = "Union",
                             rp_keys = c(threshold_type = "threshold_type",
                                         threshold_value = "threshold_value",
                                         combine_strategy = "combine_strategy",
                                         regime = "regime")) {
    key <- function(name) {
        k <- rp_keys[[name]]
        if (is.null(k)) paste0(prefix, name) else k
    }
    config_rule <- function(pd) {
        if (is.null(pd)) return(list(type = "bonf", value = 0.05))
        default_threshold(pd$config, source, registry_defaults())
    }
    list(
        threshold_type = filter_param_state(
            input, project_data, paste0(prefix, "threshold_type"), module, key("threshold_type"),
            start_fn = function(pd, saved) {
                if (!is.null(saved)) as.character(saved) else config_rule(pd)$type
            },
            coerce = function(v) {
                if (is.null(v) || length(v) != 1 || !nzchar(v)) NULL else as.character(v)
            }),
        threshold_value = filter_param_state(
            input, project_data, paste0(prefix, "threshold_value"), module, key("threshold_value"),
            start_fn = function(pd, saved) {
                if (!is.null(saved)) as.numeric(saved) else config_rule(pd)$value
            },
            coerce = function(v) {
                if (is.null(v) || length(v) != 1 || is.na(v) || v <= 0) NULL else as.numeric(v)
            }),
        strategy = filter_param_state(
            input, project_data, paste0(prefix, "combine_strategy"), module, key("combine_strategy"),
            start_fn = function(pd, saved) .normalize_strategy(saved %||% default_strategy),
            coerce = function(v) {
                if (is.null(v) || length(v) != 1 || !nzchar(v)) NULL else .normalize_strategy(v)
            }),
        regime = filter_param_state(
            input, project_data, paste0(prefix, "regime"), module, key("regime"),
            start_fn = function(pd, saved) isTRUE(saved),
            coerce = function(v) if (is.null(v)) NULL else isTRUE(v))
    )
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
