#' Discover all available projects from the pipeline path
#' @noRd
#' Discover every results directory, tagged
#'
#' Returns one row per `*_results/` directory with:
#'   name         project name (directory minus the `_results` suffix)
#'   is_simulated declared by `Simulation.enabled` in that project's own config
#'   sim_source   free-text provenance for the UI badge
#'   is_internal  scaffolding that should never be user-visible (clone lanes)
#'
#' `is_simulated` is read from a DECLARED FLAG, never guessed from the name. The
#' de-facto "starts with MVP" convention is exactly what let two different datasets
#' share the name TEST in 2026-07. `is_internal` likewise comes from an explicit
#' `.internal` marker file rather than a name pattern, because `grepl("p[0-9]+$")`
#' would also swallow a legitimately-named project.
#'
#' A project with no config of its own is NOT simulated — absence means false. Note
#' this deliberately does not use read_project_config(), which falls back to a shared
#' config.yaml; that fallback would let one project's Simulation block leak onto
#' every project that lacks its own config.
#' @noRd
discover_projects <- function(pipeline_path = get_pipeline_path()) {
    dirs  <- list.dirs(pipeline_path, full.names = FALSE, recursive = FALSE)
    dirs  <- dirs[grepl("_results$", dirs)]
    names <- sort(gsub("_results$", "", dirs))
    empty <- data.frame(name = character(0), is_simulated = logical(0),
                        sim_source = character(0), is_internal = logical(0),
                        stringsAsFactors = FALSE)
    if (!length(names)) return(empty)

    cfg_files <- file.path(pipeline_path, paste0("config_", names, ".yaml"))
    res_dirs  <- file.path(pipeline_path, paste0(names, "_results"))
    # Fingerprint over the project set, its configs AND its marker files, so the cache
    # turns over when a project appears, disappears, has its Simulation block edited, or
    # gets marked internal/simulated -- without re-parsing ~100 YAML files every render.
    fp <- paste0(length(names), "_",
                 sum(file.exists(cfg_files)), "_",
                 as.integer(max(c(0, file.mtime(cfg_files[file.exists(cfg_files)])))), "_",
                 sum(file.exists(file.path(res_dirs, ".internal"))), "_",
                 sum(file.exists(file.path(res_dirs, ".simulation"))))

    load_cached("discover_projects", function() {
        rows <- lapply(seq_along(names), function(i) {
            sim <- FALSE; src <- ""
            res_dir <- file.path(pipeline_path, paste0(names[i], "_results"))
            if (file.exists(cfg_files[i])) {
                cfg <- tryCatch(yaml::read_yaml(cfg_files[i]), error = function(e) list())
                blk <- cfg[["Simulation"]]
                if (is.list(blk)) {
                    sim <- isTRUE(blk[["enabled"]])
                    src <- if (!is.null(blk[["source"]])) as.character(blk[["source"]]) else ""
                }
            }
            # Fallback marker for results trees that have NO config_{project}.yaml of
            # their own -- sweep variants whose configs live under compound filenames
            # (config_MVP1232218M05_c1_p11.yaml and friends). The config flag is still
            # the primary declaration; this only covers trees it cannot reach, and it
            # is never consulted to override an existing config.
            if (!sim && file.exists(file.path(res_dir, ".simulation"))) {
                sim <- TRUE
                src <- tryCatch(trimws(readLines(file.path(res_dir, ".simulation"), n = 1)),
                                error = function(e) "")
            }
            data.frame(name = names[i], is_simulated = sim, sim_source = src,
                       is_internal = file.exists(file.path(res_dir, ".internal")),
                       stringsAsFactors = FALSE)
        })
        do.call(rbind, rows)
    }, fingerprint = fp)
}

#' Project names for the selector
#'
#' Analysis mode (the default) shows only real datasets; benchmark mode shows only
#' simulated ones. Internal scaffolding is hidden in both. `all = TRUE` returns every
#' name regardless -- used for duplicate-name checks when creating a project, which
#' must see names the selector hides.
#' @noRd
find_projects <- function(pipeline_path = get_pipeline_path(),
                          benchmark_mode = FALSE, all = FALSE) {
    p <- discover_projects(pipeline_path)
    if (!nrow(p)) return(character(0))
    if (all) return(sort(p$name))
    p <- p[!p$is_internal, , drop = FALSE]
    p <- p[if (isTRUE(benchmark_mode)) p$is_simulated else !p$is_simulated, , drop = FALSE]
    sort(p$name)
}

#' Every project a user may select, real and simulated alike.
#'
#' find_projects() partitions on `is_simulated`, which existed to serve the
#' navbar's Benchmark switch: one pool or the other, never both. The merged
#' chrome has no such switch (a project declares its own nature via
#' `Simulation.enabled` in its config, so it was never a session mode), and
#' without it `input$benchmark_mode` is NULL, `isTRUE(NULL)` is FALSE, and
#' simulated projects became unreachable from the selector AND from
#' app_server's repopulate observer.
#'
#' One list instead. `is_internal` scaffolding stays excluded -- that flag marks
#' things that should never be user-visible, which is a different question from
#' whether the data is simulated.
#' @noRd
find_selectable_projects <- function(pipeline_path = get_pipeline_path()) {
    p <- discover_projects(pipeline_path)
    if (!nrow(p)) return(character(0))
    sort(p$name[!p$is_internal])
}

#' Find K values for a project (from PreStructure/plots/K{k}/ directories)
#' @noRd
find_k_values <- function(project) {
    struct_plots <- mod_path(project, MOD_PRESTRUCT, "plots")
    if (!dir.exists(struct_plots)) return(integer(0))
    k_dirs <- list.dirs(struct_plots, full.names = FALSE, recursive = FALSE)
    k_nums <- as.integer(gsub("^K", "", k_dirs[grepl("^K\\d+$", k_dirs)]))
    sort(k_nums)
}

#' Find K_start and K_end from cross_entropy plot filenames
#' @noRd
find_k_range <- function(project) {
    struct_plots <- mod_path(project, MOD_PRESTRUCT, "plots")
    files <- list.files(struct_plots, pattern = "^cross_entropy_K.*\\.png$",
                        full.names = TRUE)
    if (length(files) == 0) return(list(k_start = NA, k_end = NA, path = NULL))
    # Match floats or integers: K2-10 or K2.0-10.0
    m <- regmatches(basename(files[1]),
                    regexpr("K([0-9.]+)-([0-9.]+)", basename(files[1])))
    if (length(m) == 0) return(list(k_start = NA, k_end = NA, path = files[1]))
    parts <- strsplit(gsub("^K", "", m), "-")[[1]]
    list(k_start = as.integer(as.numeric(parts[1])),
         k_end   = as.integer(as.numeric(parts[2])),
         path    = files[1])
}

#' Find all association methods available for a module
#' @noRd
find_assoc_methods <- function(project, module = MOD_GEA) {
    methods_dir <- mod_path(project, module, "tables", "methods")
    if (!dir.exists(methods_dir)) return(character(0))
    list.dirs(methods_dir, full.names = FALSE, recursive = FALSE)
}

#' Find traits from regions_per_trait.tsv (unique trait column values)
#' @noRd
find_assoc_traits <- function(project, module = MOD_GEA) {
    rpt <- regions_per_trait_path(project, module)
    if (!file.exists(rpt)) return(character(0))
    tryCatch({
        dt <- data.table::fread(rpt, select = "trait")
        sort(unique(dt$trait))
    }, error = function(e) character(0))
}

#' Find bio variable numbers from piemap filenames
#' @noRd
find_bio_values <- function(project) {
    piemap_dir <- mod_path(project, MOD_STRUCT, "plots", "piemap")
    if (!dir.exists(piemap_dir)) return(integer(0))
    files <- list.files(piemap_dir, pattern = "^piemap_bio_\\d+\\.png$")
    sort(as.integer(gsub("piemap_bio_(\\d+)\\.png", "\\1", files)))
}

#' Find GF run suffixes from maladaptation output directories
#' @noRd
find_gf_suffixes <- function(project, method = "gradient_forest") {
    method_dir <- mod_path(project, MOD_MALAD, "plots", method)
    if (!dir.exists(method_dir)) return(character(0))
    dirs <- list.dirs(method_dir, full.names = FALSE, recursive = FALSE)
    dirs[dirs != "zoom"]
}

#' Find all available maladaptation models across all methods.
#' Returns a named character vector: names = display labels, values = "method:::suffix".
#' @noRd
find_mala_models <- function(project) {
    known_methods <- c("gradient_forest", "geometric_offset", "rda_offset")
    out_vals  <- character(0)
    out_names <- character(0)
    for (m in known_methods) {
        suffs <- find_gf_suffixes(project, method = m)
        if (length(suffs) == 0) next
        method_label <- switch(m,
            gradient_forest  = "GF",
            geometric_offset = "GeoOff",
            rda_offset       = "RDAoff",
            m
        )
        vals  <- paste0(m, ":::", suffs)
        names <- paste0(method_label, " / ", suffs)
        out_vals  <- c(out_vals,  vals)
        out_names <- c(out_names, names)
    }
    if (length(out_vals) == 0) return(character(0))
    setNames(out_vals, out_names)
}

#' Parse a "method:::suffix" model key into its components.
#' Returns list(method = ..., suffix = ...) or NULL on bad input.
#' @noRd
parse_model_key <- function(key) {
    if (is.null(key) || !grepl(":::", key, fixed = TRUE)) return(NULL)
    parts <- strsplit(key, ":::", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NULL)
    list(method = parts[1], suffix = parts[2])
}

#' Find haplotype tags (validated: must have HapObject files).
#' Auto-creates the scan_done.flag if HapObjects exist but flag is missing
#' (self-healing for pipeline runs that predate the flag mechanism).
#' @noRd
find_haplotype_tags <- function(project) {
    hap_inter <- mod_path(project, MOD_INTER, "haplotype")
    if (!dir.exists(hap_inter)) return(character(0))
    tags <- list.dirs(hap_inter, full.names = FALSE, recursive = FALSE)
    valid <- vapply(tags, function(tag) {
        hap_files <- list.files(
            mod_path(project, MOD_INTER, "haplotype", tag),
            pattern = "HapObject\\.qs$"
        )
        if (length(hap_files) == 0L) return(FALSE)
        flag <- mod_path(project, MOD_INTER, "flags",
                         paste0("haplotype_", tag, "_scan_done.flag"))
        if (!file.exists(flag)) {
            flag_dir <- mod_path(project, MOD_INTER, "flags")
            dir.create(flag_dir, recursive = TRUE, showWarnings = FALSE)
            file.create(flag, showWarnings = FALSE)
        }
        TRUE
    }, logical(1))
    tags[valid]
}

#' Check which pipeline modules have completed (based on pipeline_summary.tsv steps)
#' @param project project name
#' @param summary_dt optional pre-loaded summary data.table; loaded from disk if NULL
#' @noRd
check_module_status <- function(project, summary_dt = NULL) {
    if (is.null(summary_dt)) summary_dt <- load_pipeline_summary(project)
    steps_present <- if (!is.null(summary_dt) && nrow(summary_dt) > 0) {
        unique(summary_dt$step)
    } else {
        character(0)
    }
    module_steps <- c(
        processing    = "processing",
        prestructure  = "prestructure",
        climate       = "climate",
        traits        = "traits",
        structure     = "structure",
        pregea        = "pregea",
        gea           = "gea",
        gwas          = "gwas",
        gea_x_gwas    = "gea_x_gwas",
        maladaptation = "maladaptation"
    )
    vapply(module_steps, function(s) s %in% steps_present, logical(1))
}

#' Find zoom region tags for a GF suffix
#' @noRd
find_gf_zooms <- function(project, suffix, method = "gradient_forest") {
    zoom_dir <- mod_path(project, MOD_MALAD, "plots", method, suffix, "zoom")
    if (!dir.exists(zoom_dir)) return(character(0))
    files <- list.files(zoom_dir, pattern = "\\.png$")
    gsub("\\.png$", "", files)
}

#' Find swept condition_pcs values that have per-model RDA artifacts
#' (PreGEA/plots/rda/models/pc{n}/). Backs the RDA tab's model selector —
#' one biplot/screeplot/anova set per value pregea_rda_setup.R swept.
#' Returns an integer vector sorted ascending, or empty if RDA hasn't run yet.
#' @noRd
find_rda_models <- function(project) {
    models_dir <- mod_path(project, MOD_PREGEA, "plots", "rda", "models")
    if (!dir.exists(models_dir)) return(integer(0))
    dirs <- list.dirs(models_dir, full.names = FALSE, recursive = FALSE)
    ns <- suppressWarnings(as.integer(sub("^pc", "", dirs)))
    sort(ns[!is.na(ns)])
}

#' Find haplotype regions for a tag
#' @noRd
find_haplotype_regions <- function(project, tag) {
    sel_file <- hap_selected_regions_path(project, tag)
    if (!file.exists(sel_file)) return(character(0))
    tryCatch({
        dt <- data.table::fread(sel_file, select = "region_id")
        as.character(dt$region_id)
    }, error = function(e) character(0))
}

#' Find haplotype regions that have clustree plots (scan completed)
#' @noRd
find_haplotype_scan_regions <- function(project, tag) {
    clustree_dir <- mod_path(project, MOD_INTER, "haplotype", tag, "plots", "clustree")
    if (!dir.exists(clustree_dir)) return(character(0))
    files <- list.files(clustree_dir, pattern = "^Region_(.*)_clustree_MG\\.png$")
    gsub("^Region_(.*)_clustree_MG\\.png$", "\\1", files)
}

#' Find traits available for a haplotype region (from boxplot files)
#' @noRd
find_haplotype_traits <- function(project, tag, region_id) {
    plots_dir <- mod_path(project, MOD_INTER, "haplotype", tag, "plots")
    if (!dir.exists(plots_dir)) return(character(0))
    rid_esc <- gsub("([\\-\\.])", "\\\\\\1", region_id, perl = TRUE)
    pat <- paste0("^Region_", rid_esc, "_boxplot_(.*)\\.png$")
    files <- list.files(plots_dir, pattern = pat)
    sort(gsub(pat, "\\1", files))
}

#' Find traits for a haplotype region with coordinate-overlap fallback.
#' Returns list(traits, matched_rid) or NULL if nothing found.
#' Calls find_haplotype_traits() for exact match first, then scans all
#' Region_*_boxplot_*.png files and checks coordinate overlap.
#' @noRd
.find_haplotype_traits_fuzzy <- function(project, tag, region_id) {
    # Exact match first
    exact <- find_haplotype_traits(project, tag, region_id)
    if (length(exact) > 0L) return(list(traits = exact, matched_rid = region_id))

    plots_dir <- mod_path(project, MOD_INTER, "haplotype", tag, "plots")
    if (!dir.exists(plots_dir)) return(NULL)

    all_files <- list.files(plots_dir, pattern = "^Region_.*_boxplot_.*\\.png$")
    if (length(all_files) == 0L) return(NULL)

    target <- .parse_region_id(region_id)
    if (is.null(target)) return(NULL)

    # Extract unique region ids present on disk
    rids <- unique(sub("^Region_(.+)_boxplot_.+\\.png$", "\\1", all_files))
    for (crid in rids) {
        candidate <- .parse_region_id(crid)
        if (is.null(candidate)) next
        if (.regions_overlap(target, candidate)) {
            traits <- find_haplotype_traits(project, tag, crid)
            if (length(traits) > 0L)
                return(list(traits = traits, matched_rid = crid))
        }
    }
    NULL
}
