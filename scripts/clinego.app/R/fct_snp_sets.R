# fct_snp_sets.R
# Curated SNP sets for Maladaptation: saved from the GEA tab, or declared in config
# (GEA.snp_sets / GWAS.snp_sets, built by mode=gea / mode=gwas). Both live in one store.
#
# Shared contract with the Snakemake pipeline:
#   Set store:  {PROJECT}_results/_intermediate/snp_sets/{name}/selected_snps.tsv
#               header'd TSV, columns: SNPID, chr, pos, min_pvalue; unique by SNPID,
#               ordered chr, pos, SNPID (collapse_snp_set(), scripts/R/lib/snp_set_select.R)
#   Record:     {name}/snp_set.yaml — what reproduces the set: name, origin (shiny |
#               config), source_module, regime, strategy, clumping_distance, pairs,
#               cells, k_best, n_snps, created. A Shiny record also carries the master
#               rule (threshold_type / threshold_value) and config_snippet, the
#               GEA.snp_sets block that rebuilds it. Absent for sets saved before
#               records existed.
#   Config sets are named gea.<name> / gwas.<name>; save_snp_set() refuses both prefixes.
#   Manifest:   {PROJECT}_results/_intermediate/snp_sets/manifest.json
#               array of {name, n_snps, source_module, threshold_type, threshold_value,
#               threshold_overrides, regime, strategy, clumping_distance, traits,
#               methods, pairs, k_best, created}
#               Shiny sets only, written by Shiny ONLY — the pipeline resolves sets by glob.
#   GF suffixes emitted by pipeline: {name}_spatial | {name}_nospatial

# ─── Manifest / record read ───────────────────────────────────────────────────

#' Read the SNP sets manifest (returns list of objects; empty list on missing/corrupt)
#' @noRd
read_snp_sets_manifest <- function(project) {
    if (is.null(project)) return(list())
    path <- snp_sets_manifest_path(project)
    if (!file_ok(path)) return(list())
    tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(e) list()
    )
}

#' Read a set's snp_set.yaml as parsed (no normalisation); NULL if absent or corrupt
#' @noRd
read_snp_set_record <- function(project, name) {
    if (is.null(project)) return(NULL)
    path <- snp_set_record_path(project, name)
    if (!file_ok(path)) return(NULL)
    rec <- tryCatch(yaml::read_yaml(path), error = function(e) NULL)
    if (is.list(rec)) rec else NULL
}

# ─── List for picker ──────────────────────────────────────────────────────────

#' Every set in the store: data.frame(name, n_snps, origin)
#'
#' A set is a store directory holding a selected_snps.tsv — the same rule as the
#' pipeline's glob, so a config set (no manifest entry) is listed and a manifest
#' entry whose directory is gone is not. Order: manifest order, then the rest by
#' name. Metadata from snp_set.yaml when present, else the manifest; origin is
#' "shiny" / "config" from the record, "shiny" for a manifest-only set, "unknown"
#' for a directory with neither (e.g. written by promote_snp_set.R).
#' @noRd
list_snp_sets <- function(project) {
    empty <- data.frame(name = character(0), n_snps = integer(0), origin = character(0),
                        stringsAsFactors = FALSE)
    if (is.null(project)) return(empty)
    store <- snp_sets_dir(project)
    if (!dir.exists(store)) return(empty)
    dirs <- list.dirs(store, full.names = FALSE, recursive = FALSE)
    dirs <- dirs[file.exists(file.path(store, dirs, "selected_snps.tsv"))]
    if (length(dirs) == 0) return(empty)

    man       <- read_snp_sets_manifest(project)
    man_names <- vapply(man, function(x) as.character(x$name %||% ""), character(1))
    ordered   <- c(intersect(man_names, dirs), sort(setdiff(dirs, man_names), method = "radix"))

    rows <- lapply(ordered, function(nm) {
        rec <- read_snp_set_record(project, nm)
        ent <- if (nm %in% man_names) man[[match(nm, man_names)]] else NULL
        n <- rec$n_snps %||% ent$n_snps
        if (is.null(n)) {
            n <- tryCatch(nrow(data.table::fread(snp_set_path(project, nm), select = 1L)),
                          error = function(e) NA_integer_)
        }
        origin <- rec$origin %||% (if (!is.null(ent)) "shiny" else "unknown")
        data.frame(name = nm, n_snps = as.integer(n), origin = as.character(origin),
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' The fewest SNPs a set needs to enter mode=maladaptation (the RDA offset's floor).
#' Mirrors MIN_OFFSET_SNPS in workflow/rules/common.smk; keep the two in step.
#' @noRd
MIN_OFFSET_SNPS <- 3L

#' The sets the Maladaptation picker offers: every listed set with at least
#' MIN_OFFSET_SNPS SNPs
#'
#' A config set that selected nothing is built as a header-only TSV, and a strict
#' one can select one or two SNPs. The pipeline skips such a set under
#' `snp_sets: all` and refuses it when named, so offering it would let one
#' unchecked box write a list that fails mode=maladaptation at parse time.
#' @noRd
pickable_snp_sets <- function(project) {
    sets <- list_snp_sets(project)
    sets[!is.na(sets$n_snps) & sets$n_snps >= MIN_OFFSET_SNPS, , drop = FALSE]
}

#' Picker label: "<name> (<n> SNPs) · app" / "· config"
#' @noRd
snp_set_label <- function(sets) {
    origin <- c(shiny = "app", config = "config")[sets$origin]
    origin[is.na(origin)] <- "file"
    sprintf("%s (%d SNPs) · %s", sets$name, sets$n_snps, origin)
}

#' Return TRUE if a set with this name already exists
#' @noRd
set_exists <- function(project, name) {
    isTRUE(name %in% list_snp_sets(project)$name)
}

#' Why a set name cannot be saved from the app, or NULL when it can
#' @noRd
snp_set_name_problem <- function(name) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name))
        return("Enter a name.")
    if (!grepl("^[A-Za-z0-9_.]+$", name))
        return("Only letters, digits, underscore and dot allowed.")
    if (grepl("^(gea|gwas)\\.", name, ignore.case = TRUE))
        return("Names starting with gea. or gwas. are reserved for sets declared in config.")
    NULL
}

# ─── Save ─────────────────────────────────────────────────────────────────────

#' Save a curated SNP set (writes TSV + record + atomically upserts manifest).
#'
#' @param project  Project name string.
#' @param name     Set name — must match [A-Za-z0-9_.]+ and not start with gea./gwas.
#' @param sigsnps_dt  interactive_sigsnps() result (8-col long data.table).
#' @param params_list List of manifest metadata fields (see contract above).
#' @param record   snp_set_record_for_save() result, or NULL. When given, it is
#'   written as snp_set.yaml (with name, n_snps and created filled in here) and
#'   its pairs and k_best join the manifest entry. When NULL, a record left by an
#'   earlier save under this name is removed rather than left describing other SNPs.
#' @return Invisible integer: number of unique SNPs written.
#' @noRd
save_snp_set <- function(project, name, sigsnps_dt, params_list, record = NULL) {
    # Enforce the documented rule HERE, not only in the UI caller. mod_gea.R
    # gates the modal, but the store is also written by direct calls, and a name
    # carrying a glob metacharacter is a deletion hazard for its siblings
    # (see delete_snp_set below). promote_snp_set.R:49 refuses the same pattern
    # pipeline-side, so an unvalidated name here is one the pipeline rejects.
    # gea./gwas. belong to config sets: the pipeline's snp_set_build rule would
    # overwrite such a directory on the next mode=gea.
    problem <- snp_set_name_problem(name)
    if (!is.null(problem)) stop("Invalid SNP set name '", name, "': ", problem)

    dir <- file.path(snp_sets_dir(project), name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    # Long dt → UNIQUE by SNPID, min p across methods/traits (shared lib, so a
    # set rebuilt from config writes the same bytes)
    uniq <- collapse_snp_set(sigsnps_dt)

    # Atomic TSV write
    tsv     <- snp_set_path(project, name)
    tmp_tsv <- tempfile(tmpdir = dir, fileext = ".tsv")
    data.table::fwrite(uniq, tmp_tsv, sep = "\t", quote = FALSE)
    file.rename(tmp_tsv, tsv)

    n       <- nrow(uniq)
    created <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")

    rec_path <- snp_set_record_path(project, name)
    if (!is.null(record)) {
        record$name    <- name
        record$origin  <- "shiny"
        record$n_snps  <- n
        record$created <- created
        write_snp_set_record(record, rec_path)
    } else if (file.exists(rec_path)) {
        unlink(rec_path, expand = FALSE)
    }

    # Atomic manifest upsert (replace same-name entry or append)
    entry <- c(list(
        name          = name,
        n_snps        = n,
        created       = created
    ), params_list)
    if (!is.null(record)) {
        # I(): write_json(auto_unbox = TRUE) would turn a one-pair list into a scalar
        entry$pairs  <- I(as.character(unlist(record$pairs)))
        entry$k_best <- record$k_best
    }

    man  <- read_snp_sets_manifest(project)
    keep <- Filter(function(x) !identical(x$name, name), man)
    .write_manifest_atomic(project, c(keep, list(entry)))

    invisible(n)
}

# ─── Record + config snippet ──────────────────────────────────────────────────

#' The trait::method pairs a set saved from the GEA tab records
#'
#' The matrix selection as compute_interactive_sigsnps() reads it
#' (resolve_selected_pairs()), restricted to cells that exist in the p-value
#' tables, so a stale selection cannot name a column build_snp_set.R refuses.
#' A pair that selected nothing stays: it contributes no SNPs on either side.
#' @param pvalues_list Named list method -> wide p-value table.
#' @noRd
snp_set_pairs_for_save <- function(tm_selection_json, combo_counts, known_traits,
                                   pvalues_list) {
    pairs <- resolve_selected_pairs(tm_selection_json, combo_counts, known_traits)
    sp <- split_snp_set_pairs(pairs[grepl("^.+::", pairs)])
    if (nrow(sp) == 0) return(character(0))
    exists <- mapply(function(t, m) m %in% names(pvalues_list) && t %in% names(pvalues_list[[m]]),
                     sp$trait, sp$method)
    sp$pair[exists]
}

#' Build the record of a set about to be saved from the GEA tab
#'
#' @param pairs The trait::method pairs the selection used
#'   (resolve_selected_pairs()); stored sorted, so the record does not depend
#'   on the order the matrix sent them in.
#' @param overrides,master_type,master_value,registry_defaults The tab's rule
#'   state; each pair's cell gets effective_rule_for() of it.
#' @return A record list for save_snp_set(record =), config_snippet included.
#' @noRd
snp_set_record_for_save <- function(name, source_module, pairs, overrides,
                                    master_type, master_value, registry_defaults,
                                    strategy, clumping_distance, regime, k_best) {
    pairs <- sort(unique(as.character(pairs)), method = "radix")
    sp <- split_snp_set_pairs(pairs)
    cells <- lapply(seq_len(nrow(sp)), function(i) {
        r <- effective_rule_for(sp$method[i], sp$trait[i], overrides,
                                master_type, master_value, registry_defaults)
        list(type = as.character(r$type), value = as.numeric(r$value))
    })
    names(cells) <- pairs
    rec <- list(
        name              = name,
        origin            = "shiny",
        source_module     = source_module,
        regime            = regime,
        strategy          = .normalize_strategy(strategy),
        clumping_distance = as.integer(clumping_distance),
        pairs             = pairs,
        cells             = cells,
        threshold_type    = as.character(master_type),
        threshold_value   = as.numeric(master_value)
    )
    if (length(k_best) == 1L && !is.na(k_best)) rec$k_best <- as.integer(k_best)
    rec$config_snippet <- snp_set_config_snippet(rec, registry_defaults)$yaml
    rec
}

#' The config block that rebuilds a set: YAML for GEA.snp_sets / GWAS.snp_sets
#'
#' Read by common.smk's resolve_snp_set_declarations(), so it is written in the
#' terms that resolution uses:
#'   - threshold_type / threshold_value always (the master rule; for a config
#'     record, which stores only cells, the most common rule of the univariate
#'     cells). Declaring one also pins RDA to its registry rule, as in the app.
#'   - overrides: one "trait::method" entry per cell whose rule differs from what
#'     config gives it without one (registry rule for a pinned method, else master).
#'   - traits + methods when the pairs are a full traits x methods rectangle,
#'     else the explicit pairs. Both always written out: omitting them means
#'     "every GEA.configs method", which need not be what was on disk here.
#'   - strategy and clumping_distance always.
#' traits / methods / pairs are written as sequences even when they hold one item.
#'
#' @param record A record (snp_set_record_for_save() or a snp_set.yaml).
#' @param registry_defaults method -> list(adjust=, threshold=, family=).
#' @return list(yaml = string or NULL, note = string or NULL).
#' @noRd
snp_set_config_snippet <- function(record, registry_defaults = list()) {
    if (!identical(as.character(record$regime %||% "snp"), "snp")) {
        return(list(yaml = NULL, note = paste(
            "WZA-regime sets cannot be declared in config yet (only per-SNP sets),",
            "so no config block reproduces this one.")))
    }
    group <- if (identical(record$source_module, "GWAS")) "GWAS" else "GEA"
    pairs <- as.character(unlist(record$pairs, use.names = FALSE))
    sp    <- split_snp_set_pairs(pairs)
    rule_of <- function(k) list(type  = as.character(record$cells[[k]]$type),
                                value = as.numeric(record$cells[[k]]$value))
    pinned <- function(m) {
        rd <- registry_defaults[[m]]
        if (is.null(rd) || identical(rd$family %||% "univariate_pvalue", "univariate_pvalue"))
            return(NULL)
        list(type = as.character(rd$adjust), value = as.numeric(rd$threshold))
    }
    same_rule <- function(a, b) identical(a$type, b$type) && isTRUE(a$value == b$value)

    master <- if (!is.null(record$threshold_type) && !is.null(record$threshold_value)) {
        list(type = as.character(record$threshold_type),
             value = as.numeric(record$threshold_value))
    } else {
        free <- pairs[vapply(sp$method, function(m) is.null(pinned(m)), logical(1))]
        keys <- if (length(free)) free else pairs
        sig  <- vapply(keys, function(k) {
            r <- rule_of(k); paste(r$type, format(r$value, digits = 15))
        }, character(1))
        tab <- table(factor(sig, levels = unique(sig)))
        rule_of(keys[match(names(tab)[which.max(tab)], sig)])
    }

    overrides <- list()
    for (i in seq_along(pairs)) {
        r <- rule_of(pairs[i])
        if (!same_rule(r, pinned(sp$method[i]) %||% master))
            overrides[[pairs[i]]] <- list(type = r$type, value = r$value)
    }

    name <- gsub("[^A-Za-z0-9_]", "_", sub("^(gea|gwas)\\.", "", as.character(record$name)))
    entry <- list(name = name, threshold_type = master$type, threshold_value = master$value)
    if (length(overrides)) entry$overrides <- overrides
    traits  <- unique(sp$trait)
    methods <- unique(sp$method)
    if (length(pairs) == length(traits) * length(methods)) {
        entry$methods <- as.list(methods)
        entry$traits  <- as.list(traits)
    } else {
        entry$pairs <- as.list(pairs)
    }
    entry$strategy          <- .normalize_strategy(as.character(record$strategy %||% "Union"))
    entry$clumping_distance <- as.integer(record$clumping_distance)

    block <- stats::setNames(list(list(snp_sets = list(entry))), group)
    list(yaml = yaml::as.yaml(block, precision = 15L), note = NULL)
}

#' Config block for a set saved before records existed, from its manifest entry
#'
#' The manifest keeps only the traits and methods left after combining, not the
#' cells that were selected, so the block selects every cell of those traits x
#' methods — possibly more than went into the set. The note says so.
#' @noRd
snp_set_legacy_snippet <- function(entry) {
    if (!identical(as.character(entry$regime %||% "snp"), "snp")) {
        return(list(yaml = NULL, note = paste(
            "WZA-regime sets cannot be declared in config yet (only per-SNP sets),",
            "so no config block reproduces this one.")))
    }
    traits  <- as.character(unlist(entry$traits))
    methods <- as.character(unlist(entry$methods))
    if (!length(traits) || !length(methods) || is.null(entry$threshold_type)) {
        return(list(yaml = NULL, note = paste(
            "This set was saved before SNP-set records existed and its manifest entry",
            "lacks the traits, methods or threshold needed to rebuild it.")))
    }
    group <- if (identical(entry$source_module, "GWAS")) "GWAS" else "GEA"
    out <- list(name = gsub("[^A-Za-z0-9_]", "_", as.character(entry$name)),
                threshold_type  = as.character(entry$threshold_type),
                threshold_value = as.numeric(entry$threshold_value))
    ov <- entry$threshold_overrides
    if (is.list(ov) && length(ov) && !is.null(names(ov))) {
        keep <- vapply(names(ov), function(k) {
            sp <- strsplit(k, "::", fixed = TRUE)[[1]]
            length(sp) == 2L && sp[1] %in% traits && sp[2] %in% methods
        }, logical(1))
        ov <- lapply(ov[keep], function(r) list(type = as.character(r$type),
                                                value = as.numeric(r$value)))
        if (length(ov)) out$overrides <- ov
    }
    out$methods <- as.list(methods)
    out$traits  <- as.list(traits)
    out$strategy <- .normalize_strategy(as.character(entry$strategy %||% "Union"))
    if (!is.null(entry$clumping_distance))
        out$clumping_distance <- as.integer(entry$clumping_distance)
    block <- stats::setNames(list(list(snp_sets = list(out))), group)
    list(yaml = yaml::as.yaml(block, precision = 15L), note = paste(
        "Saved before SNP-set records existed: the exact trait x method cells are",
        "unknown, so this block selects every cell of these traits x methods,",
        "which can be more than went into the set."))
}

#' Everything the Maladaptation details popover shows for one set
#'
#' @param registry_defaults Used only to build a snippet for a record that does
#'   not carry one (a config record).
#' @return list(meta = record or manifest entry, has_record, origin,
#'   snippet = list(yaml, note)), or NULL when the set is unknown.
#' @noRd
snp_set_details <- function(project, name, registry_defaults = list()) {
    rec <- read_snp_set_record(project, name)
    if (!is.null(rec)) {
        snippet <- if (!is.null(rec$config_snippet)) {
            list(yaml = rec$config_snippet, note = NULL)
        } else {
            tryCatch(snp_set_config_snippet(rec, registry_defaults),
                     error = function(e) list(yaml = NULL, note = paste(
                         "Could not build a config block:", conditionMessage(e))))
        }
        return(list(meta = rec, has_record = TRUE,
                    origin = as.character(rec$origin %||% "shiny"), snippet = snippet))
    }
    man <- Filter(function(x) identical(x$name, name), read_snp_sets_manifest(project))
    if (length(man) == 0) return(NULL)
    list(meta = man[[1]], has_record = FALSE, origin = "shiny",
         snippet = snp_set_legacy_snippet(man[[1]]))
}

# ─── UI ───────────────────────────────────────────────────────────────────────

#' The "Reproduce in config" block: the YAML in a <pre> with a copy button, or the note
#'
#' @param snippet list(yaml, note) from snp_set_config_snippet() and friends.
#' @param id DOM id for the <pre>, unique on the page.
#' @noRd
snp_set_snippet_ui <- function(snippet, id) {
    note <- if (!is.null(snippet$note))
        htmltools::p(class = "small text-muted mb-1", snippet$note)
    if (is.null(snippet$yaml)) return(note)
    # navigator.clipboard exists only in a secure context (https or localhost).
    # Reached by IP on the shared host, the button selects the text instead.
    js <- sprintf(paste0(
        "var el=document.getElementById('%s');",
        "if(navigator.clipboard&&window.isSecureContext){",
        "navigator.clipboard.writeText(el.innerText);this.innerText='Copied';",
        "}else{var r=document.createRange();r.selectNodeContents(el);",
        "var s=window.getSelection();s.removeAllRanges();s.addRange(r);",
        "this.innerText='Selected — press Ctrl+C';}"), id)
    htmltools::tagList(
        note,
        htmltools::tags$pre(
            id = id, class = "small border rounded p-2 mb-1 bg-body-tertiary",
            style = "max-height: 40vh; overflow: auto; white-space: pre;",
            snippet$yaml),
        htmltools::tags$button(
            type = "button", class = "btn btn-outline-secondary btn-sm", onclick = js,
            bsicons::bs_icon("clipboard"), " Copy")
    )
}

#' Modal shown after a set is saved from the GEA tab: what was saved + its config block
#' @noRd
show_snp_set_saved_modal <- function(ns, name, n, snippet, module = MOD_GEA) {
    mode <- if (identical(module, MOD_GWAS)) "mode=gwas" else "mode=gea"
    shiny::showModal(shiny::modalDialog(
        title = sprintf("Saved SNP set ‘%s’ (%d SNPs)", name, n),
        htmltools::p(class = "small",
            "Maladaptation lists it now. To rebuild the same set without the app, ",
            "add this block to the project config and run ", htmltools::code(mode), ":"),
        snp_set_snippet_ui(snippet, ns("snp_set_snippet")),
        footer = shiny::modalButton("Close"),
        easyClose = TRUE, size = "l"
    ))
}

# ─── Delete ───────────────────────────────────────────────────────────────────

#' Delete a curated SNP set and (optionally) its GF pipeline results.
#'
#' GF result directories use suffixes {name}_spatial and {name}_nospatial
#' (both variants from spatial_correction: both). We match EXACTLY these two
#' suffixes rather than a raw glob prefix to avoid over-deleting a set named
#' "foo" when "foo_bar" also exists.
#'
#' Every unlink() here passes expand = FALSE. The default is TRUE, which treats
#' the path as a WILDCARD: a set named "foo*" used to delete foo_bar and foobaz
#' along with itself, while the manifest below — pruned by name equality — kept
#' listing them. A store path is never a glob. save_snp_set() now refuses such
#' names, but this is what protects sets already on disk.
#'
#' A config set (gea.* / gwas.*) is deleted the same way; it comes back on the
#' next mode=gea / mode=gwas for as long as it is declared.
#' @noRd
delete_snp_set <- function(project, name, remove_gf_results = TRUE) {
    unlink(file.path(snp_sets_dir(project), name),
           recursive = TRUE, force = TRUE, expand = FALSE)

    man  <- read_snp_sets_manifest(project)
    man  <- Filter(function(x) !identical(x$name, name), man)
    .write_manifest_atomic(project, man)

    if (isTRUE(remove_gf_results)) {
        # Known suffixes the pipeline emits for this set name
        suffixes <- c(name,
                      paste0(name, "_spatial"),
                      paste0(name, "_nospatial"))
        # Remove results for all registered maladaptation methods. Must stay in
        # step with workflow/methods/maladaptation.py — the output directories
        # are keyed by REGISTRY KEY, not by engine (common.smk:1382-1389), so
        # rda_offset_corrected has its own tree despite sharing rda_offset's
        # script. Asserted by tests/testthat/test-equivalence-app-pipeline.R,
        # which greps for this assignment: keep it on one line.
        all_methods <- c("gradient_forest", "geometric_offset", "rda_offset", "rda_offset_corrected")
        for (method in all_methods) {
            for (base in c(
                mod_path(project, MOD_MALAD, "plots",  method),
                mod_path(project, MOD_MALAD, "tables", method),
                mod_path(project, MOD_INTER, method)
            )) {
                for (suf in suffixes) {
                    target <- file.path(base, suf)
                    if (file.exists(target) || dir.exists(target))
                        unlink(target, recursive = TRUE, force = TRUE,
                               expand = FALSE)
                }
            }
        }
    }

    invisible(TRUE)
}

# ─── Private helpers ──────────────────────────────────────────────────────────

#' Atomically write the manifest JSON (tempfile + rename idiom from fct_region_params.R)
#' @noRd
.write_manifest_atomic <- function(project, man) {
    path <- snp_sets_manifest_path(project)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".json")
    tryCatch({
        jsonlite::write_json(man, tmp, auto_unbox = TRUE, pretty = TRUE, digits = NA)
        file.rename(tmp, path)
    }, error = function(e) {
        if (file.exists(tmp)) unlink(tmp)
        message("WARNING: Could not write snp_sets manifest: ", conditionMessage(e))
    })
}
