# On-the-fly region computation for interactive Shiny exploration.
# Sources the shared lib so Shiny uses the same clustering logic as the pipeline.

# Guarded exactly like fct_data_loading.R's source() of pval_threshold.R, and for
# the same reason: top-level code in R/ runs during R CMD INSTALL's lazy-load
# step, where the /pipeline bind mount does not exist. Unguarded, this line threw
# "cannot open file ... regions.R", R reported "lazy loading failed" and deleted
# the half-installed package, and remotes::install_local() downgraded that to a
# warning — so `docker build` exited 0 while shipping an image with no
# clinego.app in it. dev.R source()s this file at runtime, where /pipeline is
# mounted, which is why dev mode never noticed.
if (file.exists("/pipeline/scripts/R/lib/regions.R")) {
    source("/pipeline/scripts/R/lib/regions.R")
}

#' Compute all regions from sig SNPs using single-linkage clustering
#'
#' Two SNPs merge when gap <= 2*distance, matching the pipeline's GRanges
#' extend-by-distance + reduce() behavior (lib/regions.R: cluster_snps_to_regions).
#'
#' @param sig_snps data.table with columns SNPID, chr, pos, pvalue, method, trait
#' @param distance integer. Maximum gap / 2 between neighboring SNPs (REGION_DISTANCE).
#' @return data.table with region_id, chr, start, end, length, snp_count,
#'   snp_ids, traits, methods, min_pvalue. Empty data.table if no SNPs.
#' @noRd
compute_all_regions <- function(sig_snps, distance = 2000000L) {
    distance <- as.integer(distance)
    if (is.null(sig_snps) || nrow(sig_snps) == 0) return(.empty_regions())

    # Deduplicate to one row per SNP with min pvalue
    unique_snps <- unique(sig_snps[, .(SNPID, chr = as.character(chr),
                                       pos = as.integer(pos),
                                       min_pvalue = pvalue)])
    unique_snps[, min_pvalue := min(min_pvalue, na.rm = TRUE), by = "SNPID"]
    unique_snps <- unique(unique_snps, by = "SNPID")

    # Cluster via lib (gap > 2*distance criterion, pipeline-identical)
    regions <- cluster_snps_to_regions(unique_snps, distance, trait_label = NULL)
    if (is.null(regions) || nrow(regions) == 0) return(.empty_regions())
    regions[, c("trait", "methods") := NULL]  # re-derived below from long-format

    # Collect traits and methods per region from long-format sig_snps
    snp_pos <- unique(sig_snps[, .(SNPID, chr = as.character(chr),
                                    pos = as.integer(pos), method, trait)])
    snp_pos[, s := pos][, e := pos]
    reg_ov <- data.table::copy(regions[, .(region_id, chr, start, end)])
    data.table::setkey(reg_ov, chr, start, end)
    data.table::setkey(snp_pos, chr, s, e)
    ov <- data.table::foverlaps(snp_pos, reg_ov,
        by.x = c("chr", "s", "e"),
        by.y = c("chr", "start", "end"),
        type = "within", mult = "first", nomatch = NA)

    trait_methods <- ov[!is.na(region_id), .(
        traits  = paste(sort(unique(trait)),  collapse = ","),
        methods = paste(sort(unique(method)), collapse = ",")
    ), by = "region_id"]

    regions <- trait_methods[regions, on = "region_id"]

    data.table::setcolorder(regions,
        c("region_id", "chr", "start", "end", "length",
          "snp_count", "snp_ids", "traits", "methods", "min_pvalue"))
    regions
}

#' Find genes overlapping a region (data.table-based, no GenomicRanges)
#'
#' @param gff_genes data.table from load_gff_genes(): gene_id, chr, start, end, ...
#' @param region_row single-row data.table with chr, start, end
#' @param promoter_length integer. Extend region upstream to capture promoters.
#' @return data.table of genes overlapping [start - promoter_length, end].
#'   Empty data.table if gff_genes is NULL or no overlap found.
#' @noRd
find_genes_in_region <- function(gff_genes, region_row, promoter_length = 10000L) {
    if (is.null(gff_genes) || nrow(gff_genes) == 0) return(data.table::data.table())
    if (is.null(region_row) || nrow(region_row) == 0) return(data.table::data.table())

    promoter_length <- as.integer(promoter_length)
    chr_val <- as.character(region_row$chr[1])
    r_start <- pmax(1L, as.integer(region_row$start[1]) - promoter_length)
    r_end   <- as.integer(region_row$end[1])

    # Extend region window
    query <- data.table::data.table(
        chr   = chr_val,
        start = r_start,
        end   = r_end
    )

    genes_chr <- data.table::copy(gff_genes[chr == chr_val])
    if (nrow(genes_chr) == 0) return(data.table::data.table())

    data.table::setkey(genes_chr, chr, start, end)
    data.table::setkey(query, chr, start, end)

    ov <- data.table::foverlaps(query, genes_chr,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type = "any", mult = "all", nomatch = NULL)

    if (is.null(ov) || nrow(ov) == 0) return(data.table::data.table())

    # Remove query columns (i.start, i.end), keep gene columns
    drop_cols <- intersect(c("i.start", "i.end"), names(ov))
    if (length(drop_cols) > 0) ov[, (drop_cols) := NULL]
    ov
}

#' Run GO enrichment for a single region via processx subprocess
#'
#' Calls scripts/run_enrichment.R then scripts/plot_enrichment.R.
#' Results written to a session-specific temp directory.
#'
#' @param genes_dt data.table of genes for this region (from find_genes_in_region())
#' @param region_id character combined region_id (e.g. "1_1000000-2000000")
#' @param trait character trait name (for per-trait enrichment)
#' @param project_data reactive project data bundle (name, config)
#' @return list(table, dotplot_path, emapplot_path) or NULL on failure
#' @noRd
run_enrichment_subprocess <- function(genes_dt, region_id, trait, project_data) {
    pd <- project_data

    go_field   <- config_get(pd$config, "GFF",        "go_field",   default = NULL)
    gff_feat   <- config_get(pd$config, "GFF",        "feature",    default = "mRNA")
    top_terms  <- config_get(pd$config, "Enrichment", "top_terms",  default = 10L)
    plot_w     <- config_get(pd$config, "Enrichment", "plot_width",  default = 10L)
    plot_h     <- config_get(pd$config, "Enrichment", "plot_height", default = 8L)
    cnet_label <- config_get(pd$config, "Enrichment", "cnet_label", default = "gene_id")

    if (is.null(go_field) || go_field == "" || go_field == "NULL") {
        return(list(error = "GO field not configured (GFF.go_field)"))
    }

    gff <- .resolve_gff_path(pd)
    if (!file_ok(gff)) {
        return(list(error = paste0("GFF not found: ", gff)))
    }

    if (is.null(genes_dt) || nrow(genes_dt) == 0) {
        return(list(error = "No genes in region"))
    }

    pipeline_path <- get_pipeline_path()

    # Temp session directory for this region's enrichment output
    tmp_base  <- file.path(tempdir(), "clinego_enrichment",
                           paste0(pd$name, "_", gsub("[^a-z0-9]", "", tolower(region_id)),
                                  "_", trait))
    tmp_genes   <- file.path(tmp_base, "genes.tsv")
    tmp_tables  <- file.path(tmp_base, "tables")
    tmp_inter   <- file.path(tmp_base, "intermediate")
    tmp_plots   <- file.path(tmp_base, "plots")

    dir.create(tmp_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_inter,  recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_plots,  recursive = TRUE, showWarnings = FALSE)

    # Write per-trait genes file (run_enrichment.R expects region_id with trait suffix)
    per_trait_id <- paste0(region_id, "_", trait)
    genes_out <- data.table::copy(genes_dt)
    genes_out[, region_id := per_trait_id]
    data.table::fwrite(genes_out, tmp_genes, sep = "\t")

    # ── Step 1: Run enrichment ─────────────────────────────────────────────────
    res_enrich <- tryCatch(
        processx::run(
            command = "Rscript",
            args    = c(
                file.path(pipeline_path, "scripts", "run_enrichment.R"),
                tmp_genes,
                gff,
                go_field,
                gff_feat,
                tmp_tables,
                tmp_inter
            ),
            wd      = pipeline_path,
            timeout = 120
        ),
        error = function(e) list(status = 1L, stderr = conditionMessage(e))
    )

    if (res_enrich$status != 0L) {
        return(list(error = paste0("Enrichment failed: ",
                                   utils::tail(strsplit(res_enrich$stderr, "\n")[[1]], 3))))
    }

    # Check if TSV was produced
    tsv_path <- file.path(tmp_tables, trait,
                          paste0("Region_", per_trait_id, "_enrichment.tsv"))
    if (!file_ok(tsv_path)) {
        return(list(error = "No enrichment results (no significant GO terms)"))
    }

    enrich_tbl <- tryCatch(
        data.table::fread(tsv_path, sep = "\t", header = TRUE),
        error = function(e) NULL
    )
    if (is.null(enrich_tbl) || nrow(enrich_tbl) == 0) {
        return(list(error = "No enrichment results (no significant GO terms)"))
    }

    # ── Step 2: Generate plots ─────────────────────────────────────────────────
    res_plots <- tryCatch(
        processx::run(
            command = "Rscript",
            args    = c(
                file.path(pipeline_path, "scripts", "plot_enrichment.R"),
                tmp_inter,
                tmp_plots,
                as.character(as.integer(top_terms)),
                as.character(as.numeric(plot_w)),
                as.character(as.numeric(plot_h)),
                cnet_label,
                tmp_genes,
                "0"  # top_plot_regions = 0 (all)
            ),
            wd      = pipeline_path,
            timeout = 120
        ),
        error = function(e) list(status = 1L, stderr = conditionMessage(e))
    )

    dotplot_path  <- file.path(tmp_plots, trait,
                               paste0("Region_", per_trait_id, "_dotplot.png"))
    emapplot_path <- file.path(tmp_plots, trait,
                               paste0("Region_", per_trait_id, "_emapplot.png"))

    list(
        table         = enrich_tbl,
        dotplot_path  = if (file_ok(dotplot_path))  dotplot_path  else NULL,
        emapplot_path = if (file_ok(emapplot_path)) emapplot_path else NULL
    )
}

#' Run regionplot for a single region via processx subprocess
#'
#' Calls scripts/plot_regionplot.R with CUSTOM_REGION set.
#'
#' @param region_row single-row data.table with chr, start, end
#' @param trait character trait name
#' @param project_data project data bundle
#' @param module character pipeline module (MOD_GEA or MOD_GWAS)
#' @return character path to generated PNG, or NULL on failure
#' @noRd
run_regionplot_subprocess <- function(region_row, trait, project_data, module = MOD_GEA) {
    pd <- project_data

    # Convert GFF on-the-fly if the topr annotation doesn't exist yet
    gff_topr_p <- .ensure_gff_topr(pd)
    if (!file_ok(gff_topr_p)) {
        return(list(error = "Could not find or generate topr gene annotation. Check GFF config."))
    }

    # Build ASSOC_TABLES: method:adjust:filepath using the full pvalue TSVs
    method_files <- find_method_pvalue_files(pd$name, module)
    if (length(method_files) == 0) {
        return(list(error = "No method result files found. Run the association mode first."))
    }

    assoc_tables <- .build_assoc_tables_arg(method_files)
    if (is.null(assoc_tables)) {
        return(list(error = "Could not build association table arguments."))
    }

    custom_region <- paste0(region_row$chr[1], ":",
                             region_row$start[1], "-",
                             region_row$end[1])

    genes_highlight <- config_get(pd$config, "regionplot", "genes",
                                  default = "NULL") %||% "NULL"
    if (length(genes_highlight) > 1) genes_highlight <- paste(genes_highlight, collapse = "|")

    pipeline_path <- get_pipeline_path()

    tmp_base <- file.path(tempdir(), "clinego_regionplot",
                          paste0(pd$name, "_",
                                 gsub("[^a-z0-9]", "", tolower(custom_region)),
                                 "_", trait))
    dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)

    res <- tryCatch(
        processx::run(
            command = "Rscript",
            args    = c(
                file.path(pipeline_path, "scripts", "plot_regionplot.R"),
                "NULL",          # REGIONS_FILE (not needed for custom region)
                gff_topr_p,
                assoc_tables,
                "0",             # TOP_REGIONS
                genes_highlight,
                paste0(tmp_base, "/"),
                custom_region,
                trait,           # CUSTOM_TRAITS
                "NULL"           # CUSTOM_METHODS (all)
            ),
            wd      = pipeline_path,
            timeout = 180
        ),
        error = function(e) list(status = 1L, stderr = conditionMessage(e))
    )

    if (res$status != 0L) {
        stderr_tail <- paste(utils::tail(strsplit(res$stderr %||% "", "\n")[[1]], 5), collapse = "\n")
        return(list(error = paste0("Regionplot failed:\n", stderr_tail)))
    }

    # Find generated PNG in tmp_base
    pngs <- Sys.glob(file.path(tmp_base, "*.png"))
    if (length(pngs) == 0) {
        return(list(error = "No plot was generated. Check that the region contains SNPs in the p-value files."))
    }
    list(path = pngs[1])
}

# ── Async subprocess launchers ────────────────────────────────────────────────

#' Launch GO enrichment as a non-blocking background process.
#'
#' Returns a handle list immediately; the caller polls \code{handle$process$is_alive()}.
#' On completion, pass the handle to \code{.collect_enrichment_result()}.
#'
#' @return list(process, log_file, type, tmp_tables, tmp_plots, per_trait_id, trait)
#'   or list(error=) if pre-flight validation fails.
#' @noRd
launch_enrichment_subprocess <- function(genes_dt, region_id, trait, project_data, module = MOD_GEA) {
    pd <- project_data

    go_field   <- config_get(pd$config, "GFF",        "go_field",   default = NULL)
    gff_feat   <- config_get(pd$config, "GFF",        "feature",    default = "mRNA")
    top_terms  <- config_get(pd$config, "Enrichment", "top_terms",  default = 10L)
    plot_w     <- config_get(pd$config, "Enrichment", "plot_width",  default = 10L)
    plot_h     <- config_get(pd$config, "Enrichment", "plot_height", default = 8L)
    cnet_label <- config_get(pd$config, "Enrichment", "cnet_label", default = "gene_id")

    if (is.null(go_field) || go_field == "" || go_field == "NULL")
        return(list(error = "GO field not configured (GFF.go_field)"))

    gff <- .resolve_gff_path(pd)
    if (!file_ok(gff))
        return(list(error = paste0("GFF not found: ", gff)))

    if (is.null(genes_dt) || nrow(genes_dt) == 0)
        return(list(error = "No genes in region"))

    pipeline_path <- get_pipeline_path()
    # Intermediate and temp files stay in tempdir — only final tables/plots persist
    tmp_base   <- file.path(tempdir(), "clinego_enrichment",
                            paste0(pd$name, "_", gsub("[^a-z0-9]", "", tolower(region_id)),
                                   "_", trait))
    tmp_genes  <- file.path(tmp_base, "genes.tsv")
    tmp_inter  <- file.path(tmp_base, "intermediate")
    log_file   <- file.path(tmp_base, "enrichment.log")
    # Persistent output dirs — survive app restarts; match enrichment_table_path / enrichment_plot_path
    tmp_tables <- mod_path(pd$name, module, "tables", "enrichment")
    tmp_plots  <- mod_path(pd$name, module, "plots", "enrichment")

    dir.create(tmp_inter,  recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_plots,  recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(tmp_genes), recursive = TRUE, showWarnings = FALSE)

    per_trait_id <- paste0(region_id, "_", trait)
    genes_out <- data.table::copy(genes_dt)
    genes_out[, region_id := per_trait_id]
    data.table::fwrite(genes_out, tmp_genes, sep = "\t")

    enrich_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "run_enrichment.R"),
        tmp_genes, gff, go_field, gff_feat, tmp_tables, tmp_inter
    ), shQuote), collapse = " ")

    plot_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "plot_enrichment.R"),
        tmp_inter, tmp_plots,
        as.character(as.integer(top_terms)),
        as.character(as.numeric(plot_w)),
        as.character(as.numeric(plot_h)),
        cnet_label, tmp_genes, "0"
    ), shQuote), collapse = " ")

    # Chain both steps; log file captures stdout+stderr of both
    cmd <- paste0(
        "Rscript ", enrich_args,
        " >> ", shQuote(log_file), " 2>&1",
        " && Rscript ", plot_args,
        " >> ", shQuote(log_file), " 2>&1"
    )

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch enrichment process"))

    list(
        process      = proc,
        log_file     = log_file,
        type         = "enrichment",
        tmp_tables   = tmp_tables,
        tmp_plots    = tmp_plots,
        per_trait_id = per_trait_id,
        trait        = trait
    )
}

#' Collect enrichment result after the background process has finished.
#' @noRd
.collect_enrichment_result <- function(handle) {
    tsv_path <- file.path(handle$tmp_tables, handle$trait,
                          paste0("Region_", handle$per_trait_id, "_enrichment.tsv"))
    if (!file_ok(tsv_path))
        return(list(error = "No enrichment results (no significant GO terms)"))

    enrich_tbl <- tryCatch(
        data.table::fread(tsv_path, sep = "\t", header = TRUE),
        error = function(e) NULL
    )
    if (is.null(enrich_tbl) || nrow(enrich_tbl) == 0)
        return(list(error = "No enrichment results (no significant GO terms)"))

    dotplot_path  <- file.path(handle$tmp_plots, handle$trait,
                               paste0("Region_", handle$per_trait_id, "_dotplot.png"))
    emapplot_path <- file.path(handle$tmp_plots, handle$trait,
                               paste0("Region_", handle$per_trait_id, "_emapplot.png"))

    list(
        table         = enrich_tbl,
        dotplot_path  = if (file_ok(dotplot_path))  dotplot_path  else NULL,
        emapplot_path = if (file_ok(emapplot_path)) emapplot_path else NULL
    )
}

#' Build trait:method1|method2,trait2:method3 arg from region sig SNPs.
#' @noRd
.build_trait_method_arg <- function(region_snps) {
    if (is.null(region_snps) || nrow(region_snps) == 0) return(NULL)
    pairs <- unique(region_snps[, .(trait, method)])
    if (nrow(pairs) == 0) return(NULL)
    # Group methods by trait
    by_trait <- split(pairs$method, pairs$trait)
    paste(mapply(function(trait, methods) {
        paste0(trait, ":", paste(sort(unique(methods)), collapse = "|"))
    }, names(by_trait), by_trait), collapse = ",")
}

#' Launch regionplot as a non-blocking background process.
#' @param region_snps data.table with columns trait, method — sig SNPs in the region
#' @return list(process, log_file, type, tmp_base) or list(error=)
#' @noRd
launch_regionplot_subprocess <- function(region_row, region_snps, project_data, module = MOD_GEA,
                                         regime = "snp") {
    pd <- project_data

    gff_topr_p <- .ensure_gff_topr(pd)
    if (!file_ok(gff_topr_p))
        return(list(error = "Could not find or generate topr gene annotation. Check GFF config."))

    # For MOD_GEAXGWAS, pull method pvalue files from both GEA + GWAS sources.
    method_files <- if (identical(module, MOD_GEAXGWAS)) {
        c(find_method_pvalue_files(pd$name, MOD_GEA),
          find_method_pvalue_files(pd$name, MOD_GWAS))
    } else {
        find_method_pvalue_files(pd$name, module)
    }
    if (length(method_files) == 0)
        return(list(error = "No method result files found. Run the association mode first."))
    # Filter to k_best only — glob returns files for all K values, causing each method
    # to appear twice in the regionplot legend (once per K) with different colors
    k_best <- pd$k_best
    if (!is.na(k_best)) {
        k_pattern <- paste0("_K", k_best, ".")
        method_files <- method_files[grepl(k_pattern, basename(method_files), fixed = TRUE)]
    }

    assoc_tables <- .build_assoc_tables_arg(method_files)
    if (is.null(assoc_tables))
        return(list(error = "Could not build association table arguments."))

    # Build structured trait:method pairs from region sig SNPs
    custom_traits <- .build_trait_method_arg(region_snps)
    if (is.null(custom_traits))
        return(list(error = "No significant trait/method pairs found in this region."))

    # Hash of trait/method combination for per-filter caching
    combo_hash <- substr(digest::digest(custom_traits, algo = "md5"), 1, 8)

    custom_region <- paste0(region_row$chr[1], ":", region_row$start[1], "-", region_row$end[1])
    region_safe <- gsub("[:\\-]", "_", custom_region)
    genes_highlight <- config_get(pd$config, "regionplot", "genes", default = "NULL") %||% "NULL"
    if (length(genes_highlight) > 1) genes_highlight <- paste(genes_highlight, collapse = "|")

    pipeline_path <- get_pipeline_path()
    # Persistent output dir — survives app restarts; subdirectory per combo_hash
    out_dir <- ondemand_regionplot_dir(pd$name, module, combo_hash)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    # Log stays in tempdir to avoid polluting results
    log_file <- file.path(tempdir(), paste0("rplot_", gsub("[^a-z0-9]", "", tolower(custom_region)), ".log"))

    rplot_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "plot_regionplot.R"),
        "NULL", gff_topr_p, assoc_tables, "0",
        genes_highlight, paste0(out_dir, "/"),
        custom_region, custom_traits, "NULL",
        regime
    ), shQuote), collapse = " ")

    cmd <- paste0("Rscript ", rplot_args, " >> ", shQuote(log_file), " 2>&1")

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch regionplot process"))

    list(
        process     = proc,
        log_file    = log_file,
        type        = "regionplot",
        out_dir     = out_dir,
        region_safe = region_safe
    )
}

#' Collect regionplot result after the background process has finished.
#' @noRd
.collect_regionplot_result <- function(handle) {
    expected_png <- file.path(handle$out_dir,
                              paste0("regionplot_custom_", handle$region_safe, ".png"))
    if (!file_ok(expected_png))
        return(list(error = "No plot was generated. Check that the region contains SNPs in the p-value files."))
    list(path = expected_png)
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Empty regions skeleton
#' @noRd
.empty_regions <- function() {
    data.table::data.table(
        region_id  = character(),
        chr        = character(),
        start      = integer(),
        end        = integer(),
        length     = integer(),
        snp_count  = integer(),
        snp_ids    = character(),
        traits     = character(),
        methods    = character(),
        min_pvalue = numeric()
    )
}

#' Resolve full GFF path from project config
#' @noRd
.resolve_gff_path <- function(pd) {
    cfg     <- pd$config
    inp_dir <- config_get(cfg, "Input", "dir",  default = "data")
    gff_rel <- config_get(cfg, "Input", "gff",  default = "")
    if (is.null(gff_rel) || gff_rel == "") return("")
    file.path(get_pipeline_path(), inp_dir, gff_rel)
}

#' Resolve topr-converted GFF path (from intermediate directory)
#' @noRd
.resolve_gff_topr_path <- function(pd) {
    # Pipeline places gff2topr output in _intermediate/annotation/
    Sys.glob(mod_path(pd$name, MOD_INTER, "annotation", "*.tsv"))[1] %||% ""
}

#' Ensure topr-format GFF exists, converting from normalized GFF3 if needed.
#' Returns the path to the topr TSV, or "" on failure.
#' @noRd
.ensure_gff_topr <- function(pd) {
    # Return immediately if already exists (pipeline ran gff2topr rule)
    existing <- .resolve_gff_topr_path(pd)
    if (file_ok(existing)) return(existing)

    # Locate the normalized GFF3 written by the processing mode
    norm_gff <- mod_path(pd$name, MOD_INTER, "annotation", "normalized.gff3")
    if (!file_ok(norm_gff)) return("")

    # Read GFF config values
    cfg       <- pd$config
    feature   <- config_get(cfg, "GFF", "feature",   default = "mRNA")
    gene_name <- config_get(cfg, "GFF", "gene_name", default = "description")
    biotype   <- config_get(cfg, "GFF", "biotype",   default = "biotype")

    out_path <- mod_path(pd$name, MOD_INTER, "annotation", "topr_gene_annotation.tsv")
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

    pipeline_path <- get_pipeline_path()
    script <- file.path(pipeline_path, "scripts", "gff2topr.py")
    if (!file.exists(script)) return("")

    res <- tryCatch(
        processx::run(
            command = "python3",
            args    = c(script, norm_gff, feature, gene_name, biotype, out_path),
            wd      = pipeline_path,
            timeout = 30
        ),
        error = function(e) list(status = 1L)
    )

    if (res$status == 0L && file_ok(out_path)) out_path else ""
}

# ── Haplotype inline helpers ──────────────────────────────────────────────────

#' Count SNPs in a VCF within a chr:start-end region (via awk).
#' Returns NA_integer_ if VCF missing or awk fails.
#' @noRd
count_vcf_region_snps <- function(vcf_path, chr, start, end) {
    if (is.null(vcf_path) || !nzchar(vcf_path) || !file_ok(vcf_path))
        return(NA_integer_)
    cmd <- sprintf(
        "awk -v c=%s -v s=%d -v e=%d 'BEGIN{n=0} !/^#/ && $1==c && $2>=s && $2<=e {n++} END{print n}' %s",
        shQuote(as.character(chr)),
        as.integer(start), as.integer(end),
        shQuote(vcf_path)
    )
    out <- suppressWarnings(as.integer(system(cmd, intern = TRUE, ignore.stderr = TRUE)))
    if (length(out) == 0 || is.na(out)) NA_integer_ else out
}

#' Module → haplotype tag source suffix mapping
#' @noRd
.module_to_hap_source <- function(module) {
    switch(module,
        GWAS     = "gwas",
        GEAxGWAS = "gea",  # GEAxGWAS regions use GEA haplotype source
        NULL
    )
}

#' Find the first haplotype tag matching the current module's source.
#' Returns NULL if no matching tag found.
#' @noRd
find_matching_hap_tag <- function(project, module) {
    source_suffix <- .module_to_hap_source(module)
    if (is.null(source_suffix)) return(NULL)
    all_tags <- find_haplotype_tags(project)
    if (length(all_tags) == 0) return(NULL)
    # Tag format: "{meta_type}_{source}" e.g. "site_association", "cluster_K5_association"
    # meta_type may contain underscores (e.g. "cluster_K5"), so match by suffix instead of split.
    suffix_pat <- paste0("_", source_suffix)
    matched <- Filter(function(tag) endsWith(tag, suffix_pat), all_tags)
    if (length(matched) == 0) NULL else matched[1]
}

# ── Region-id parsing helpers (used for fuzzy overlap matching) ────────────────

#' Parse a region_id string "{chr}_{start}-{end}" into its components.
#' Returns list(chr, start, end) or NULL if format doesn't match.
#' @noRd
.parse_region_id <- function(rid) {
    i <- regexpr("_", rid)
    if (i < 1L) return(NULL)
    chr <- substr(rid, 1L, i - 1L)
    se  <- strsplit(substr(rid, i + 1L, nchar(rid)), "-", fixed = TRUE)[[1]]
    if (length(se) != 2L) return(NULL)
    s <- suppressWarnings(as.integer(se[[1L]]))
    e <- suppressWarnings(as.integer(se[[2L]]))
    if (is.na(s) || is.na(e)) return(NULL)
    list(chr = chr, start = s, end = e)
}

#' Check whether two parsed regions (from .parse_region_id) overlap.
#' @noRd
.regions_overlap <- function(r1, r2) {
    r1$chr == r2$chr && r1$start <= r2$end && r2$start <= r1$end
}

#' Check if haplotype scan results exist for a region (clustree PNG present).
#' Returns list(found = TRUE/FALSE, matched_rid = region_id_used_on_disk).
#' When pipeline region differs from Shiny region, matched_rid holds the pipeline's id.
#' @noRd
check_hap_scan_results <- function(project, tag, region_id) {
    if (is.null(tag) || is.null(region_id))
        return(list(found = FALSE, matched_rid = NULL))

    # Exact match: clustree PNG (full success) OR HapObject (partial: 1 epsilon)
    inter_dir <- hap_intermediate_dir(project, tag)
    if (file_ok(hap_clustree_path(project, tag, region_id, "MG")) ||
        file_ok(file.path(inter_dir, paste0("Region_", region_id, "_HapObject.qs"))))
        return(list(found = TRUE, matched_rid = region_id))

    # Fuzzy: scan all HapObjects for coordinate overlap with the requested region
    if (!dir.exists(inter_dir)) return(list(found = FALSE, matched_rid = NULL))
    hap_files <- list.files(inter_dir, pattern = "^Region_.*_HapObject\\.qs$")
    if (length(hap_files) == 0L) return(list(found = FALSE, matched_rid = NULL))
    target <- .parse_region_id(region_id)
    if (is.null(target)) return(list(found = FALSE, matched_rid = NULL))
    for (f in hap_files) {
        crid <- sub("^Region_(.+)_HapObject\\.qs$", "\\1", f)
        candidate <- .parse_region_id(crid)
        if (!is.null(candidate) && .regions_overlap(target, candidate))
            return(list(found = TRUE, matched_rid = crid))
    }
    list(found = FALSE, matched_rid = NULL)
}

#' Check if haplotype viz results exist for a region.
#' Returns list(traits = character(), matched_rid = region_id_used_on_disk).
#' @noRd
check_hap_viz_results <- function(project, tag, region_id) {
    if (is.null(tag) || is.null(region_id))
        return(list(traits = character(0), matched_rid = region_id))
    res <- .find_haplotype_traits_fuzzy(project, tag, region_id)
    if (is.null(res)) return(list(traits = character(0), matched_rid = region_id))
    res
}

#' Launch haplotype scan as a non-blocking background process.
#'
#' Creates a single-row selected_regions.tsv in a temp dir and runs
#' run_haplotype_scan.R for exactly this region.
#'
#' @param region_row   single-row data.table with region_id, chr, start, end, snp_count
#' @param project_data project data bundle (name, config, k_best)
#' @param tag          haplotype tag string (e.g. "site_association")
#' @param params       list with: epsilon_range (chr), mgmin, minhap, min_snps, meta_type
#' @return list(process, log_file, type, plots_dir, inter_dir, region_id, tag)
#'   or list(error=) on pre-flight failure.
#' @noRd
launch_hap_scan_subprocess <- function(region_row, project_data, tag, params) {
    pd <- project_data

    filtered_vcf <- hap_filtered_vcf_path(pd$name, pd$config)
    if (!file_ok(filtered_vcf))
        return(list(error = paste0("Filtered VCF not found. Run mode=processing first.")))

    ld_vcf <- hap_imputed_vcf_path(pd$name, pd$k_best, pd$config)
    if (!file_ok(ld_vcf)) ld_vcf <- filtered_vcf  # fallback: use filtered VCF for LD

    metadata <- hap_metadata_path(pd$name)
    if (!file_ok(metadata))
        return(list(error = "Aligned metadata not found. Run mode=processing first."))

    meta_type    <- params$meta_type %||% "site"
    clusters_file <- if (grepl("^cluster", meta_type)) {
        cp <- hap_clusters_path(pd$name, pd$k_best)
        if (file_ok(cp)) cp else "NULL"
    } else "NULL"

    inter_dir  <- hap_intermediate_dir(pd$name, tag)
    plots_dir  <- hap_scan_plots_dir(pd$name, tag)
    dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

    # Single-row selected_regions.tsv in temp dir
    tmp_base <- file.path(tempdir(), "clinego_hapscan",
                          paste0(pd$name, "_", tag, "_",
                                 gsub("[^a-z0-9]", "", tolower(region_row$region_id[1]))))
    dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)
    regions_file <- file.path(tmp_base, "selected_regions.tsv")
    data.table::fwrite(region_row, regions_file, sep = "\t")

    log_file <- file.path(tmp_base, "hap_scan.log")
    pipeline_path <- get_pipeline_path()

    scan_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "run_haplotype_scan.R"),
        regions_file,
        filtered_vcf,
        ld_vcf,
        metadata,
        meta_type,
        clusters_file,
        as.character(as.integer(params$mgmin %||% 50L)),
        as.character(as.integer(params$minhap %||% 15L)),
        as.character(params$epsilon_range %||% "0.3,0.5,0.7,0.9"),
        as.character(as.integer(params$min_snps %||% 3L)),
        paste0(inter_dir, "/"),
        paste0(plots_dir, "/")
    ), shQuote), collapse = " ")

    cmd <- paste0("Rscript ", scan_args, " >> ", shQuote(log_file), " 2>&1")

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch haplotype scan process"))

    list(
        process   = proc,
        log_file  = log_file,
        type      = "hap_scan",
        plots_dir = plots_dir,
        inter_dir = inter_dir,
        region_id = region_row$region_id[1],
        tag       = tag,
        params    = params  # stored for persistence on completion
    )
}

#' Collect haplotype scan result after process finishes.
#' @noRd
.collect_hap_scan_result <- function(handle) {
    rid      <- handle$region_id
    pdir     <- handle$plots_dir
    mg_path  <- file.path(pdir, paste0("Region_", rid, "_clustree_MG.png"))
    hap_path <- file.path(pdir, paste0("Region_", rid, "_clustree_hap.png"))
    if (!file_ok(mg_path)) {
        # Check if HapObject exists — scan succeeded with only 1 epsilon (clustree skipped)
        hap_obj <- file.path(handle$inter_dir, paste0("Region_", rid, "_HapObject.qs"))
        if (file_ok(hap_obj)) {
            # Try to read epsilon from scan_status.tsv
            eps_hint <- tryCatch({
                status_file <- file.path(handle$inter_dir, "scan_status.tsv")
                if (file_ok(status_file)) {
                    st  <- data.table::fread(status_file, sep = "\t", header = TRUE)
                    row <- st[st$region_id == rid, ]
                    if (nrow(row) > 0 && !is.na(row$epsilon_succeeded[1]) && nzchar(row$epsilon_succeeded[1])) {
                        as.numeric(gsub(".*_E", "", row$epsilon_succeeded[1]))
                    } else NULL
                } else NULL
            }, error = function(e) NULL)
            return(list(skip_clustree = TRUE, epsilon_hint = eps_hint, region_id = rid,
                        epsilons = if (!is.null(eps_hint)) eps_hint else numeric(0)))
        }
        # No HapObject — scan failed; try to surface reason from scan_status.tsv
        detail_msg <- tryCatch({
            status_file <- file.path(handle$inter_dir, "scan_status.tsv")
            if (file_ok(status_file)) {
                st  <- data.table::fread(status_file, sep = "\t", header = TRUE)
                row <- st[st$region_id == rid, ]
                if (nrow(row) > 0) {
                    n   <- row$snp_count[1]
                    msg <- row$failure_reason[1]
                    paste0("Region has ", n, " SNPs in VCF. ",
                           if (!is.na(msg) && nzchar(msg)) msg else
                           "Check epsilon_range and min_snps settings.")
                } else NULL
            } else NULL
        }, error = function(e) NULL)
        return(list(error = detail_msg %||% paste0(
            "No clustree plots generated. Check that the region has enough SNPs ",
            "(\u2265 min_snps) and try adjusting epsilon_range.")))
    }

    # Discover per-trait clustree variants (Region_{rid}_clustree_MG_{trait}.png)
    trait_files <- Sys.glob(file.path(pdir,
                                       paste0("Region_", rid, "_clustree_MG_*.png")))
    traits <- if (length(trait_files) > 0) {
        fnames <- basename(trait_files)
        sub("\\.png$", "", sub(paste0("^Region_", rid, "_clustree_MG_"), "", fnames))
    } else character(0)

    # Parse successful epsilons from scan_status.tsv
    epsilons <- tryCatch({
        status_file <- file.path(handle$inter_dir, "scan_status.tsv")
        if (file_ok(status_file)) {
            st  <- data.table::fread(status_file, sep = "\t", header = TRUE)
            row <- st[st$region_id == rid, ]
            if (nrow(row) > 0 && !is.na(row$epsilon_succeeded[1]) && nzchar(row$epsilon_succeeded[1])) {
                vals <- unlist(strsplit(row$epsilon_succeeded[1], ","))
                sort(as.numeric(gsub(".*_E", "", vals)))
            } else numeric(0)
        } else numeric(0)
    }, error = function(e) numeric(0))

    list(
        mg_path   = mg_path,
        hap_path  = if (file_ok(hap_path)) hap_path else NULL,
        traits    = traits,
        plots_dir = pdir,
        region_id = rid,
        epsilons  = epsilons
    )
}

#' Launch haplotype visualization as a non-blocking background process.
#'
#' @param region_row   single-row data.table with region_id, chr, start, end
#' @param project_data project data bundle (name, config, k_best)
#' @param tag          haplotype tag
#' @param params       list with: epsilon_selected (numeric), meta_type
#' @return list(process, log_file, type, region_id, tag) or list(error=)
#' @noRd
launch_hap_viz_subprocess <- function(region_row, project_data, tag, params) {
    pd <- project_data

    epsilon <- as.numeric(params$epsilon_selected)
    if (is.null(epsilon) || is.na(epsilon))
        return(list(error = "Please enter an epsilon value (e.g. 0.6)."))

    inter_dir  <- hap_intermediate_dir(pd$name, tag)
    hapobj <- Sys.glob(file.path(inter_dir,
                                  paste0("Region_", region_row$region_id[1], "_HapObject.qs")))
    if (length(hapobj) == 0)
        return(list(error = "HapObject not found. Run the Scan step first."))

    metadata <- hap_metadata_path(pd$name)
    if (!file_ok(metadata))
        return(list(error = "Aligned metadata not found. Run mode=processing first."))

    meta_type    <- params$meta_type %||% "site"
    clusters_file <- if (grepl("^cluster", meta_type)) {
        cp <- hap_clusters_path(pd$name, pd$k_best)
        if (file_ok(cp)) cp else "NULL"
    } else "NULL"

    plots_dir  <- hap_viz_plots_dir(pd$name, tag)
    tables_dir <- hap_viz_tables_dir(pd$name, tag)
    dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

    # Single-row selected_regions.tsv
    tmp_base <- file.path(tempdir(), "clinego_hapviz",
                          paste0(pd$name, "_", tag, "_",
                                 gsub("[^a-z0-9]", "", tolower(region_row$region_id[1]))))
    dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)
    regions_file <- file.path(tmp_base, "selected_regions.tsv")
    data.table::fwrite(region_row, regions_file, sep = "\t")

    # Raster path (use first found climate raster or "NULL")
    raster_path  <- Sys.glob(file.path(project_base(pd$name),
                                        "climate", "rasters", "present", "*.tif"))[1] %||% "NULL"
    if (is.null(raster_path) || !file_ok(raster_path)) raster_path <- "NULL"
    raster_layer_raw <- config_get(pd$config, "climate", "predictors", default = "bio_1")
    raster_layer <- strsplit(as.character(raster_layer_raw[[1]]), ",")[[1]][1]

    pie_alpha      <- config_get(pd$config, "piemap", "alpha",        default = 0.8)
    show_labels    <- config_get(pd$config, "piemap", "show_labels",  default = "F")
    label_size     <- config_get(pd$config, "piemap", "label_size",   default = 3)
    pie_scale      <- config_get(pd$config, "piemap", "pie_scale",    default = 1)
    use_points     <- config_get(pd$config, "piemap", "use_points",   default = "F")
    regionmap_ext  <- config_get(pd$config, "map", "zoom_extent",     default = "NULL") %||% "NULL"
    piemap_args    <- paste(pie_alpha, show_labels, label_size, pie_scale, use_points, sep = ",")

    log_file <- file.path(tmp_base, "hap_viz.log")
    pipeline_path <- get_pipeline_path()

    viz_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "run_haplotype_viz.R"),
        regions_file,
        metadata,
        as.character(epsilon),
        meta_type,
        tag,
        paste0(inter_dir, "/"),
        paste0(plots_dir, "/"),
        paste0(tables_dir, "/"),
        raster_path,
        raster_layer,
        piemap_args,
        regionmap_ext,
        clusters_file
    ), shQuote), collapse = " ")

    cmd <- paste0("Rscript ", viz_args, " >> ", shQuote(log_file), " 2>&1")

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch haplotype visualization process"))

    list(
        process   = proc,
        log_file  = log_file,
        type      = "hap_viz",
        plots_dir = plots_dir,
        region_id = region_row$region_id[1],
        tag       = tag,
        params    = list(epsilon_selected = epsilon, meta_type = meta_type)  # for persistence
    )
}

#' Collect haplotype viz result after process finishes.
#' @noRd
.collect_hap_viz_result <- function(handle) {
    rid    <- handle$region_id
    pdir   <- handle$plots_dir
    # Glob for any trait's crosshap_viz PNG for this region
    vizs   <- Sys.glob(file.path(pdir, paste0("Region_", rid, "_crosshap_viz_*.png")))
    if (length(vizs) == 0)
        return(list(error = paste0(
            "No visualization plots generated for region ", rid,
            ". Try a different epsilon or check the log for errors.")))
    # Extract available traits from filenames
    traits <- gsub(paste0("^Region_", rid, "_crosshap_viz_(.*)\\.png$"),
                   "\\1", basename(vizs))
    list(traits = sort(traits))
}

# ── Pairwise trait overlap (on-demand async) ──────────────────────────────────

#' Flatten a named per-method sig-SNP list into a long TSV in tempdir.
#' Input list: names = method names, values = data.tables with columns
#' SNPID, chr, pos, pvalue, method, trait (as returned by compute_method_sigsnps_cached).
#' Returns the path to the written file, or NULL if no SNPs.
#' @noRd
.write_sigsnps_tmp <- function(sig_list, source_label) {
    if (length(sig_list) == 0) return(NULL)
    rows <- lapply(names(sig_list), function(m) {
        dt <- sig_list[[m]]
        if (is.null(dt) || nrow(dt) == 0) return(NULL)
        # Ensure the method column is set (compute_method_sigsnps_cached adds it)
        dt <- data.table::copy(dt)
        if (!"method" %in% names(dt)) dt[, method := m]
        dt[, .(SNPID, chr, pos, pvalue, method, trait)]
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    flat <- data.table::rbindlist(rows, use.names = TRUE)
    if (nrow(flat) == 0) return(NULL)
    tmp <- tempfile(pattern = paste0("pw_sigsnps_", source_label, "_"), fileext = ".tsv")
    data.table::fwrite(flat, tmp, sep = "\t", quote = FALSE)
    tmp
}

#' Launch the pairwise overlap computation as a non-blocking background process.
#' @param gea_sig_list Named list of per-method sig-SNP data.tables (GEA side).
#' @param gwas_sig_list Named list of per-method sig-SNP data.tables (GWAS side).
#' @param gea_window Clumping distance from the GEA filter bar (bp).
#' @param gwas_window Clumping distance from the GWAS filter bar (bp).
#' @param min_snps Minimum SNPs per side per cluster for a window overlap.
#' @param hash 8-char MD5 hash of the current parameter set (used as output dir key).
#' @param project_data Current project_data() bundle.
#' @return list(process, log_file, type="pairwise", out_dir, gea_window, gwas_window)
#'   or list(error=) on pre-flight failure.
#' @noRd
launch_pairwise_subprocess <- function(gea_sig_list, gwas_sig_list,
                                        gea_window, gwas_window, min_snps,
                                        hash, project_data) {
    pd <- project_data

    n_gea  <- sum(vapply(gea_sig_list,  function(d) if (!is.null(d)) nrow(d) else 0L, integer(1)))
    n_gwas <- sum(vapply(gwas_sig_list, function(d) if (!is.null(d)) nrow(d) else 0L, integer(1)))
    if (n_gea + n_gwas == 0)
        return(list(error = "No significant SNPs at the current thresholds. Adjust the filter bars."))

    gea_tmp  <- .write_sigsnps_tmp(gea_sig_list,  "climate")
    gwas_tmp <- .write_sigsnps_tmp(gwas_sig_list, "phenotype")

    # At least one side must have sig SNPs
    gea_arg  <- if (!is.null(gea_tmp))  gea_tmp  else "NULL"
    gwas_arg <- if (!is.null(gwas_tmp)) gwas_tmp else "NULL"

    pipeline_path <- get_pipeline_path()

    out_dir <- ondemand_pairwise_dir(pd$name, hash)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    out_collapsed <- file.path(out_dir, "pairwise_collapsed_snps.tsv")
    out_table     <- file.path(out_dir, "pairwise_overlap_table.tsv")
    log_file      <- file.path(tempdir(), paste0("pairwise_", hash, ".log"))

    script_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "compute_pairwise_ondemand.R"),
        gea_arg, gwas_arg,
        as.character(gea_window), as.character(gwas_window), as.character(min_snps),
        out_collapsed, out_table
    ), shQuote), collapse = " ")

    cmd <- paste0("Rscript ", script_args, " >> ", shQuote(log_file), " 2>&1")

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch pairwise computation process."))

    list(
        process    = proc,
        log_file   = log_file,
        type       = "pairwise",
        out_dir    = out_dir,
        gea_window = gea_window,
        gwas_window = gwas_window
    )
}

#' Collect pairwise result after the background process has finished.
#' @return list(table=dt, collapsed=dt, windows=list(gea=, gwas=)) or list(error=)
#' @noRd
.collect_pairwise_result <- function(handle) {
    tbl_path  <- file.path(handle$out_dir, "pairwise_overlap_table.tsv")
    coll_path <- file.path(handle$out_dir, "pairwise_collapsed_snps.tsv")
    if (!file_ok(tbl_path) || !file_ok(coll_path))
        return(list(error = "Pairwise output files not found. Check the computation log."))
    tbl  <- tryCatch(
        data.table::fread(tbl_path,  sep = "\t",
                          colClasses = c(trait_a="character", trait_b="character",
                                         source_a="character", source_b="character",
                                         comparison_type="character")),
        error = function(e) NULL
    )
    coll <- tryCatch(
        data.table::fread(coll_path, sep = "\t",
                          colClasses = c(SNPID="character", chr="character",
                                         trait="character", source="character")),
        error = function(e) NULL
    )
    if (is.null(tbl) || is.null(coll))
        return(list(error = "Could not read pairwise output files."))
    list(
        table     = tbl,
        collapsed = coll,
        windows   = list(gea = handle$gea_window, gwas = handle$gwas_window)
    )
}

#' Build colon-separated ASSOC_TABLES arg for plot_regionplot.R
#' Format: "method:adjust_str:pvalues_K{k}.tsv,..."
#' Accepts plain pvalue TSVs ({method}_pvalues_K{k}.tsv).
#' @noRd
.build_assoc_tables_arg <- function(method_files) {
    if (length(method_files) == 0) return(NULL)
    entries <- lapply(method_files, function(f) {
        m <- regmatches(basename(f),
            regexec("^(.+?)_pvalues_(K[0-9.]+)\\.tsv$", basename(f)))[[1]]
        if (length(m) < 3) return(NULL)
        method <- m[2]
        # Use a fixed adjust token; the regionplot background is threshold-independent
        paste(method, "bonf_0.05", f, sep = ":")
    })
    entries <- Filter(Negate(is.null), entries)
    if (length(entries) == 0) return(NULL)
    paste(entries, collapse = ",")
}
