# Source shared pipeline utilities for interactive threshold computation.
# compute_pval_threshold() is used directly; lib/sig_snps.R is not needed
# here (we implement a simpler single-threaded version inline).
if (file.exists("/pipeline/scripts/R/utils/pval_threshold.R")) {
    source("/pipeline/scripts/R/utils/pval_threshold.R")
}

# ── Processing QC loaders ────────────────────────────────────────────────────

#' Load filtering summary table
#' @noRd
load_filtering_summary <- function(project) {
    p <- qc_table_path(project, "filtering_summary.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE),
        error = function(e) data.table::data.table()
    )
}

#' Load PreGEA recommendations as a method -> param -> list(value=, evidence=) lookup.
#' Feeds the GEA config sidebar's method-param editor (pre-fill/Apply badges) —
#' recommendation keys (method, param) match the gea.py registry's own params
#' spec names 1:1 (LFMM.K, EMMAX.n_pcs, RDA.condition_pcs/axes/predictor_set),
#' since pregea_recommend.R was written against that same registry vocabulary.
#' @noRd
pregea_recommendations_lookup <- function(project) {
    p <- pregea_table_path(project, "", "pregea_recommendations")
    if (!file_ok(p)) return(list())
    dt <- tryCatch(data.table::fread(p, sep = "\t", header = TRUE),
                   error = function(e) data.table::data.table())
    if (nrow(dt) == 0 || !all(c("method", "param", "recommended_value") %in% names(dt))) return(list())
    out <- list()
    for (i in seq_len(nrow(dt))) {
        m <- dt$method[i]; pr <- dt$param[i]
        if (is.null(out[[m]])) out[[m]] <- list()
        out[[m]][[pr]] <- list(
            value    = dt$recommended_value[i],
            evidence = if ("evidence_metric" %in% names(dt))
                paste0(dt$evidence_metric[i], "=", dt$evidence_value[i]) else NULL
        )
    }
    out
}

#' Load LD decay half-distance summary table (half-decay bp, r2=0.2 bp, per group/scope)
#' @noRd
load_ld_decay_table <- function(project) {
    p <- ld_decay_table_path(project)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(group = "character", scope = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load sample heterozygosity table
#' @noRd
load_sample_heterozygosity <- function(project) {
    p <- qc_table_path(project, "sample_heterozygosity.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(sample = "character", site = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load relatedness pairs table (plink IBS allele-sharing, always produced by Processing)
#' @noRd
load_relatedness_pairs <- function(project) {
    p <- qc_table_path(project, "relatedness_pairs.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(IID1 = "character", IID2 = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load related-samples removal preview table. Always produced once Filter.relatedness
#' is set (even in 'keep' mode, as a preview) — see relatedness_filter rule.
#' @noRd
load_relatedness_removed <- function(project) {
    p <- qc_table_path(project, "relatedness_removed.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(IID = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load sNMF Q-matrix (clusters_K{k}.tsv: sample, site, C1..CK — soft ancestry
#' proportions, no discrete cluster column). @noRd
load_clusters <- function(project, k) {
    p <- hap_clusters_path(project, k)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(sample = "character", site = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Summarize a Q-matrix into per-cluster sample counts via argmax over the C1..CK
#' columns (there is no discrete cluster column on disk — see load_clusters()).
#' Returns NULL when dt is empty. @noRd
cluster_pop_summary <- function(dt, k) {
    if (nrow(dt) == 0) return(NULL)
    q_cols <- paste0("C", seq_len(k))
    q_cols <- intersect(q_cols, names(dt))
    if (length(q_cols) == 0) return(NULL)
    assign <- q_cols[max.col(as.matrix(dt[, ..q_cols]), ties.method = "first")]
    counts <- table(assign)
    list(
        n_samples = nrow(dt),
        n_pops    = length(counts),
        min_n     = as.integer(min(counts)),
        max_n     = as.integer(max(counts))
    )
}

#' Load pipeline summary TSV
#' @noRd
load_pipeline_summary <- function(project) {
    p <- pipeline_summary_path(project)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE),
        error = function(e) data.table::data.table()
    )
}

#' Load regions per trait
#' @noRd
load_regions_per_trait <- function(project, module = MOD_GEA) {
    p <- regions_per_trait_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", region_id = "character"))
        dt
    }, error = function(e) data.table::data.table())
}

#' Load combined regions
#' @noRd
load_regions_combined <- function(project, module = MOD_GEA) {
    p <- regions_combined_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character"))
    }, error = function(e) data.table::data.table())
}

#' Assign region_id to a SNP data.table via coordinate-based foverlaps lookup.
#'
#' Uses regions_combined.tsv to map each SNP position to its enclosing region.
#' Region IDs are in combined format (chr_start-end, no trait suffix).
#' @param snps_dt data.table with columns SNPID, chr, pos
#' @param project project name
#' @param module pipeline module (MOD_GEA, MOD_GWAS, MOD_GEAXGWAS)
#' @return snps_dt with region_id column added (NA for SNPs outside any region)
#' @noRd
assign_region_ids <- function(snps_dt, project, module) {
    reg_combined_path <- regions_combined_path(project, module)
    snps_dt[, region_id := NA_character_]
    if (!file.exists(reg_combined_path)) return(snps_dt)
    # Cache the regions file read so repeated calls (e.g., threshold changes re-triggering
    # interactive_sigsnps) don't re-read from disk. Fingerprint = mtime so a pipeline
    # re-run that writes new regions automatically invalidates the cache entry.
    fp   <- tryCatch(as.character(file.info(reg_combined_path)$mtime), error = function(e) "unknown")
    regs <- load_cached(
        paste0("regions_combined_", project, "_", module),
        function() tryCatch(
            data.table::fread(reg_combined_path, sep = "\t", header = TRUE,
                              colClasses = c(chr = "character", region_id = "character")),
            error = function(e) data.table::data.table()
        ),
        fingerprint = fp
    )
    if (nrow(regs) == 0 || !all(c("chr", "start", "end", "region_id") %in% names(regs)))
        return(snps_dt)

    regs_ov <- data.table::data.table(
        chr       = as.character(regs$chr),
        start     = as.integer(regs$start),
        end       = as.integer(regs$end),
        region_id = as.character(regs$region_id)
    )
    snp_pos <- unique(snps_dt[, .(SNPID, chr = as.character(chr), pos = as.integer(pos))])
    snp_ov  <- data.table::data.table(
        SNPID = snp_pos$SNPID,
        chr   = snp_pos$chr,
        start = snp_pos$pos,
        end   = snp_pos$pos
    )
    data.table::setkey(regs_ov, chr, start, end)
    data.table::setkey(snp_ov,  chr, start, end)
    ov <- data.table::foverlaps(snp_ov, regs_ov,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type = "within", mult = "first", nomatch = NA)
    rid_map <- ov[, .(SNPID, region_id)]
    snps_dt[rid_map, region_id := i.region_id, on = "SNPID"]
    snps_dt[, region_id := as.character(region_id)]
    snps_dt
}

#' Load selected SNPs table (all sig SNPs) and pivot wide→long
#'
#' Pipeline outputs wide format: SNPID, chr, pos, {METHOD1, METHOD2, ...}, min_pvalue
#' where method columns contain the trait name (or "" if not significant for that method).
#' Returns long format: SNPID, chr, pos, pvalue, method, trait, region_id
#' @noRd
load_selected_snps <- function(project, module = MOD_GEA, k_best = NA) {
    p <- selected_snps_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", SNPID = "character"))

        # Identify method columns (everything except the fixed columns)
        fixed_cols  <- c("SNPID", "chr", "pos", "min_pvalue")
        method_cols <- setdiff(names(dt), fixed_cols)

        if (length(method_cols) == 0) return(data.table::data.table())

        # Pivot wide→long: method column name → method, value → trait
        long <- data.table::melt(dt,
            id.vars       = c("SNPID", "chr", "pos", "min_pvalue"),
            measure.vars  = method_cols,
            variable.name = "method",
            value.name    = "trait",
            variable.factor = FALSE
        )

        # Keep only rows where this method flagged the SNP (trait != "")
        long <- long[trait != ""]

        # Some methods (e.g. LFMM) store comma-separated traits when a SNP is
        # significant for multiple traits. Split those into one row per trait.
        if (any(grepl(",", long$trait, fixed = TRUE))) {
            long <- long[, {
                trs <- unlist(strsplit(trait, ",", fixed = TRUE))
                .(SNPID = SNPID, chr = chr, pos = pos, min_pvalue = min_pvalue,
                  method = method, trait = trimws(trs))
            }, by = seq_len(nrow(long))][, seq_len := NULL]
        }

        # Join actual per-method p-values from sig_snps files.
        # min_pvalue is kept for combined Manhattan y-range but NOT shown in tables.
        method_dts <- load_all_method_sigsnps(project, module, k = k_best)
        if (length(method_dts) > 0) {
            method_pvals <- data.table::rbindlist(method_dts, use.names = TRUE)
            method_pvals <- method_pvals[, .(SNPID, method, trait, pvalue)]
            long <- merge(long, method_pvals,
                          by = c("SNPID", "method", "trait"), all.x = TRUE)
            # Rows without a match keep min_pvalue as fallback
            long[is.na(pvalue), pvalue := min_pvalue]
        } else {
            long[, pvalue := min_pvalue]
        }
        # Keep min_pvalue for combined Manhattan overlay y-placement
        # (combined background PNG y-axis is calibrated to min across all methods)

        long <- assign_region_ids(long, project, module)
        long
    }, error = function(e) data.table::data.table())
}

#' Load per-method sig SNPs (with actual per-method p-values)
#'
#' Reads the per-method sig_snps TSV produced by find_sig_snps.R, which contains
#' the actual p-value for that specific method (not min across all methods).
#' This is required for correct y-axis placement in per-method Manhattan overlays,
#' because each method's background PNG has a y-range calibrated to its own p-values.
#' Returns long format: SNPID, chr, pos, pvalue, method, trait, region_id
#' @noRd
load_method_sigsnps <- function(project, module = MOD_GEA, method, k, adjust) {
    p <- method_sigsnps_direct_path(project, module, method, k, adjust)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", SNPID = "character",
                                               method = "character", trait = "character"))
        if (nrow(dt) == 0) return(data.table::data.table())

        # Keep only the columns needed for the overlay
        keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(dt))
        dt <- dt[, ..keep]

        dt <- assign_region_ids(dt, project, module)
        dt
    }, error = function(e) data.table::data.table())
}

#' Load RDA diagnostics table (long key/value) as a named list.
#' Returns list() on a missing/old project (no RDA method configured yet).
#' @noRd
load_rda_diagnostics <- function(project, module = MOD_GEA) {
    hits <- Sys.glob(rda_diagnostics_path(project, module))
    if (length(hits) == 0) return(list())
    tryCatch({
        dt <- data.table::fread(hits[[1]], sep = "\t", header = TRUE, colClasses = "character")
        stats::setNames(as.list(dt$value), dt$key)
    }, error = function(e) list())
}

#' Load geometric-offset conditioning diagnostics (long quantity/value) as a named list.
#' Returns list() when geometric_offset has not been run for this suffix.
#' @noRd
load_geometric_offset_diagnostics <- function(project, suffix) {
    p <- geometric_offset_diagnostics_path(project, suffix)
    if (!file_ok(p)) return(list())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE, colClasses = "character")
        stats::setNames(as.list(dt$value), dt$quantity)
    }, error = function(e) list())
}

#' Load the Gradient-Forest imputation-sensitivity diagnostics.
#'
#' Returns list() when the check has not been run for this suffix (older projects, or
#' Maladaptation.methods.gradient_forest.sensitivity_check: false).
#' @noRd
load_gf_sensitivity <- function(project, suffix, method = "gradient_forest") {
    p <- gf_sensitivity_path(project, suffix, method)
    if (!file_ok(p)) return(list())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE, colClasses = "character")
        stats::setNames(as.list(dt$value), dt$quantity)
    }, error = function(e) list())
}

#' Load RDA candidates side table (Mahalanobis/p/q/loadings/assigned predictor).
#' @noRd
load_rda_candidates <- function(project, module = MOD_GEA) {
    hits <- Sys.glob(rda_candidates_path(project, module))
    if (length(hits) == 0) return(data.table::data.table())
    tryCatch(
        data.table::fread(hits[[1]], sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", SNPID = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load RDA anova table (full/by-axis/by-margin, long format with a status column).
#' @noRd
load_rda_anova <- function(project, module = MOD_GEA) {
    hits <- Sys.glob(rda_anova_path(project, module))
    if (length(hits) == 0) return(data.table::data.table())
    tryCatch(
        data.table::fread(hits[[1]], sep = "\t", header = TRUE),
        error = function(e) data.table::data.table()
    )
}

#' Load genes per region (collapsed)
#' @noRd
load_genes <- function(project, module = MOD_GEA) {
    p <- genes_per_region_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", region_id = "character"))
        # genes_per_region_collapsed uses per-trait region_ids ({chr}_{start}-{end}_{trait}).
        # Normalize to combined format ({chr}_{start}-{end}) to match region dropdown.
        if ("region_id" %in% names(dt) && nrow(dt) > 0) {
            dt[, region_id := sub("_[A-Za-z].*$", "", region_id)]
            dt <- unique(dt, by = intersect(c("gene_id", "region_id"), names(dt)))
        }
        dt
    }, error = function(e) data.table::data.table())
}

#' Load GO enrichment table for a specific region and trait
#' @noRd
load_enrichment_table <- function(project, module = MOD_GEA, trait, region_id) {
    p <- enrichment_table_path(project, module, trait, region_id)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE)
    }, error = function(e) data.table::data.table())
}

#' Load haplotype scan status table
#' @noRd
load_haplotype_status <- function(project, tag) {
    p <- hap_scan_status_path(project, tag)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(region_id = "character", chr = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load GF site offset table
#' @noRd
load_gf_site_table <- function(project, suffix) {
    p <- gf_site_table_path(project, suffix)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(site = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load pairwise trait overlap table (pipeline-computed)
#' @noRd
load_pairwise_table <- function(project) {
    p <- pairwise_table_path(project)
    if (!file_ok(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
            colClasses = c(trait_a = "character", source_a = "character",
                           trait_b = "character", source_b = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load pairwise collapsed sig SNPs (long format: one row per SNP per trait)
#' @noRd
load_pairwise_collapsed <- function(project) {
    p <- pairwise_collapsed_path(project)
    if (!file_ok(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
            colClasses = c(chr = "character", SNPID = "character",
                           trait = "character", source = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Parse GFF9 attribute string into a data.table of key-value columns
#'
#' Each row's attributes are parsed independently, so variable key order and
#' missing keys (filled with NA) are handled correctly.
#'
#' @param attr_strings character vector of raw GFF attribute strings
#' @return data.table with one column per attribute key
#' @noRd
parse_gff_attributes <- function(attr_strings) {
    parsed <- lapply(attr_strings, function(s) {
        if (is.na(s) || !nzchar(s)) return(list())
        pairs <- strsplit(s, ";", fixed = TRUE)[[1]]
        pairs <- pairs[nzchar(pairs)]
        keys  <- sub("=.*",     "", pairs)
        vals  <- sub("^[^=]+=", "", pairs)
        stats::setNames(as.list(vals), keys)
    })
    data.table::rbindlist(parsed, fill = TRUE)
}

#' Load and parse the project GFF as a gene coordinate table
#'
#' Reads the GFF file, filters to the configured feature type (default "mRNA"),
#' extracts gene_id from the Parent= or ID= attribute, and parses all attribute
#' fields into columns. Cached per project to avoid re-loading the full GFF on
#' every reactive trigger.
#'
#' @param project project name
#' @param config parsed YAML config (from project_data()$config)
#' @return data.table with gene_id, chr, start, end, + attribute columns.
#'   Empty data.table if GFF not found or not readable.
#' @noRd
load_gff_genes <- function(project, config) {
    cache_key <- paste0("gff_genes_", project)
    load_cached(cache_key, function() {
        inp_dir  <- config_get(config, "Input", "dir",     default = "data")
        gff_rel  <- config_get(config, "Input", "gff",     default = "")
        feat     <- config_get(config, "GFF",   "feature", default = "mRNA")

        if (is.null(gff_rel) || gff_rel == "") return(data.table::data.table())

        gff_abs <- file.path(get_pipeline_path(), inp_dir, gff_rel)
        if (!file.exists(gff_abs)) return(data.table::data.table())

        tryCatch({
            # Read GFF, skip comment lines, filter to target feature
            dt <- data.table::fread(
                cmd    = paste0("grep -v '^#' ", shQuote(gff_abs)),
                header = FALSE,
                col.names = c("chr", "source", "feature", "start", "end",
                              "score", "strand", "phase", "attributes"),
                colClasses = c("character", "character", "character",
                               "integer", "integer",
                               "character", "character", "character", "character")
            )
            dt <- dt[feature == feat]
            if (nrow(dt) == 0) return(data.table::data.table())

            dt[, chr   := sub("^chr", "", as.character(chr))]
            dt[, start := as.integer(start)]
            dt[, end   := as.integer(end)]

            # Extract gene_id: prefer Parent=, fall back to ID=
            dt[, gene_id := {
                id <- regmatches(attributes,
                    regexpr("(?<=Parent=)[^;]+", attributes, perl = TRUE))
                na_idx <- !nzchar(id)
                if (any(na_idx)) {
                    id2 <- regmatches(attributes,
                        regexpr("(?<=ID=)[^;]+", attributes, perl = TRUE))
                    id[na_idx] <- id2[na_idx]
                }
                # Strip isoform suffixes (.1, .2, _exon_1 etc.)
                sub("\\.[0-9]+(_exon_[0-9]+)?$", "", id)
            }]

            # Parse attributes into individual columns
            attr_dt <- parse_gff_attributes(dt$attributes)
            # Drop ID/Parent — redundant with gene_id
            drop_cols <- intersect(c("ID", "Parent"), names(attr_dt))
            if (length(drop_cols) > 0) attr_dt[, (drop_cols) := NULL]
            cbind(dt[, .(gene_id, chr, start, end)], attr_dt)
        }, error = function(e) data.table::data.table())
    })
}

#' Load all per-method sig SNPs for a module (for interactive Shiny combine)
#'
#' Reads all {method}_pvalues_K*_sig_snps_*.tsv files found under
#' {module}/tables/methods/, returning a named list (method → data.table).
#' Each data.table has columns: SNPID, chr, pos, pvalue, method, trait.
#' Returns an empty list when no files exist (old pipeline run / fallback).
#' @noRd
load_all_method_sigsnps <- function(project, module = MOD_GEA, k = NA) {
    files <- find_method_sigsnps_files(project, module)
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, "_")
        files <- files[grepl(k_pattern, basename(files), fixed = TRUE)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        # Extract method from path: .../tables/methods/{method}/...
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr    = "character",
                                             SNPID  = "character",
                                             method = "character",
                                             trait  = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(dt))
        if (length(keep) < 5L) next
        result[[method_name]] <- dt[, ..keep]
    }
    result
}

#' Load all per-method WZA sig-windows (WZA regime counterpart of load_all_method_sigsnps)
#' @noRd
load_all_method_wza_sigwindows <- function(project, module = MOD_GEA, k = NA) {
    glob <- mod_path(project, module, "tables", "methods", "*",
                     "*_wza_K*_sig_windows_*.tsv")
    files <- Sys.glob(glob)
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, "_")
        files <- files[grepl(k_pattern, basename(files), fixed = TRUE)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr    = "character",
                                             SNPID  = "character",
                                             method = "character",
                                             trait  = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(dt))
        if (length(keep) < 5L) next
        result[[method_name]] <- dt[, ..keep]
    }
    result
}

#' Load phenotype missing summary (GWAS mode)
#' @noRd
load_pheno_missing_summary <- function(project) {
    p <- pheno_missing_summary_path(project)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(trait = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load trait names from any pvalue TSV (column names after chr/pos/SNPID)
#' Used to detect traits that have pvalue output but zero significant SNPs.
#' Returns character(0) when no pvalue files exist.
#' @noRd
load_all_trait_names <- function(project, module = MOD_GEA) {
    f <- find_pvalue_tsv(project, module)
    if (!nzchar(f) || !file.exists(f)) return(character(0))
    tryCatch({
        hdr <- data.table::fread(f, nrows = 0L)
        fixed <- c("chr", "pos", "SNPID")
        setdiff(names(hdr), fixed)
    }, error = function(e) character(0))
}

#' Load all per-method sig SNPs for the Overlapping tab
#'
#' Combines per-method sig SNPs from both MOD_GEA and MOD_GWAS.
#' GEA traits (bio_*) and GWAS traits (phenotype columns) are naturally distinct,
#' so same method names across modules are merged without prefixing.
#'
#' @param project Character: project name
#' @param k       Optional K value to filter by (NA = all)
#' @return Named list of data.tables (method -> dt), same format as load_all_method_sigsnps()
#' @noRd
load_all_overlap_method_sigsnps <- function(project, k = NA) {
    gea  <- load_all_method_sigsnps(project, MOD_GEA,  k)
    gwas <- load_all_method_sigsnps(project, MOD_GWAS, k)

    all_methods <- union(names(gea), names(gwas))
    result <- list()
    for (m in all_methods) {
        parts <- Filter(
            function(x) !is.null(x) && nrow(x) > 0,
            list(gea[[m]], gwas[[m]])
        )
        if (length(parts) > 0)
            result[[m]] <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
    }
    result
}

#' Load full p-value TSVs for all methods (wide format: one column per trait).
#'
#' Returns a named list (method → data.table with SNPID, chr, pos, trait1..N).
#' Excludes _sig_snps_, _wza_, and _block_ files.
#' @noRd
load_all_method_pvalues <- function(project, module = MOD_GEA, k = NA) {
    glob <- mod_path(project, module, "tables", "methods", "*", "*_pvalues_K*.tsv")
    files <- Sys.glob(glob)
    # Exclude downstream outputs (sig_snps, wza, block) and GWAS per-trait files
    files <- grep("_sig_snps_|_wza_|_block_", files, value = TRUE, invert = TRUE)
    files <- keep_merged_method_files(files, regime = "snp")
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, ".tsv")
        files <- files[endsWith(basename(files), k_pattern)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr = "character", SNPID = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        result[[method_name]] <- dt
    }
    result
}

#' qs2-backed disk cache for load_all_method_pvalues.
#'
#' On first call (or after a pipeline re-run that changes the file fingerprint),
#' parses the TSVs and serialises the result list with qs2::qd_save() to
#' {project}_results/_intermediate/qs_cache/.  On subsequent calls in the same
#' or a new R session, qs2::qd_read() returns the pre-parsed list in ~milliseconds
#' regardless of the 200 MB cachem in-memory cap.
#'
#' @noRd
load_all_method_pvalues_cached <- function(project, module = MOD_GEA, k = NA) {
    fp         <- pvalues_file_fingerprint(project, module, k, regime = "snp")
    k_str      <- if (is.na(k)) "NA" else as.character(k)
    cache_dir  <- file.path(project_base(project), MOD_INTER, "qs_cache")
    cache_file <- file.path(cache_dir,
                            paste0("pvalues_", module, "_K", k_str, "_", fp, ".qs2"))

    if (file.exists(cache_file)) {
        result <- tryCatch(qs2::qd_read(cache_file),
                           error = function(e) NULL)
        if (!is.null(result)) return(result)
    }

    result <- load_all_method_pvalues(project, module, k)

    tryCatch({
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        qs2::qd_save(result, cache_file)
    }, error = function(e) {
        warning("qs2 cache write failed (pvalues): ", conditionMessage(e))
    })

    result
}

#' Load full WZA p-value TSVs for all methods (wide format: one column per trait).
#' @noRd
load_all_method_wza_pvalues <- function(project, module = MOD_GEA, k = NA) {
    glob <- mod_path(project, module, "tables", "methods", "*", "*_wza_K*.tsv")
    files <- Sys.glob(glob)
    files <- grep("_sig_windows_", files, value = TRUE, invert = TRUE)
    files <- keep_merged_method_files(files, regime = "wza")
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, ".tsv")
        files <- files[endsWith(basename(files), k_pattern)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr = "character", SNPID = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        result[[method_name]] <- dt
    }
    result
}

#' qs2-backed disk cache for load_all_method_wza_pvalues.
#' Same caching strategy as load_all_method_pvalues_cached — see that function.
#' @noRd
load_all_method_wza_pvalues_cached <- function(project, module = MOD_GEA, k = NA) {
    fp         <- pvalues_file_fingerprint(project, module, k, regime = "wza")
    k_str      <- if (is.na(k)) "NA" else as.character(k)
    cache_dir  <- file.path(project_base(project), MOD_INTER, "qs_cache")
    cache_file <- file.path(cache_dir,
                            paste0("wza_", module, "_K", k_str, "_", fp, ".qs2"))

    if (file.exists(cache_file)) {
        result <- tryCatch(qs2::qd_read(cache_file),
                           error = function(e) NULL)
        if (!is.null(result)) return(result)
    }

    result <- load_all_method_wza_pvalues(project, module, k)

    tryCatch({
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        qs2::qd_save(result, cache_file)
    }, error = function(e) {
        warning("qs2 cache write failed (wza): ", conditionMessage(e))
    })

    result
}

#' Summarise the SNP -> WZA-window collapse for a module.
#'
#' Pure helper that reads an already-loaded WZA p-value list (from
#' load_all_method_wza_pvalues()) and returns the collapse statistics.
#' Window positions are identical across methods (positional fixed tiles),
#' so we only inspect the first method's table.
#'
#' @param wza_list named list (method -> wide WZA data.table)
#' @return list(n_windows, n_snps, window_bp) or NULL when no WZA data
#' @noRd
wza_collapse_stats <- function(wza_list) {
    if (length(wza_list) == 0) return(NULL)
    dt <- wza_list[[1]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    n_snps <- if ("n_snps" %in% names(dt)) sum(dt$n_snps, na.rm = TRUE) else NA_integer_
    win_bp <- NA_integer_
    m <- regmatches(dt$SNPID[1], regexec("^.*:(\\d+)-(\\d+)$", dt$SNPID[1]))[[1]]
    if (length(m) == 3) win_bp <- as.integer(m[3]) - as.integer(m[2]) + 1L
    list(n_windows = nrow(dt), n_snps = n_snps, window_bp = win_bp)
}

#' Compute a fingerprint string for a set of files (mtime + size).
#'
#' Used to detect when source p-value files have been regenerated by the pipeline,
#' so disk caches keyed on this fingerprint automatically miss after a Run.
#' @noRd
pvalues_file_fingerprint <- function(project, module, k, regime = "snp") {
    suffix_pat <- if (regime == "wza") "*_wza_K*.tsv" else "*_pvalues_K*.tsv"
    glob_pat   <- mod_path(project, module, "tables", "methods", "*", suffix_pat)
    files      <- Sys.glob(glob_pat)
    files      <- grep("_sig_snps_|_block_|_sig_windows_", files,
                       value = TRUE, invert = TRUE)
    files      <- keep_merged_method_files(files, regime = regime)
    if (!is.na(k)) {
        k_pat  <- paste0("_K", k, ".tsv")
        files  <- files[endsWith(basename(files), k_pat)]
    }
    if (length(files) == 0) return("empty")
    fi <- file.info(files[order(files)])  # stable sort so fingerprint is deterministic
    digest::digest(list(mtime = fi$mtime, size = fi$size), algo = "md5")
}

#' Combined file fingerprint for both GEA and GWAS pvalue files.
#' Used as part of the pairwise_params hash so a pipeline re-run invalidates cached results.
#' @noRd
pairwise_file_fingerprint <- function(project) {
    fp_gea  <- tryCatch(pvalues_file_fingerprint(project, MOD_GEA,  NA), error = function(e) "err")
    fp_gwas <- tryCatch(pvalues_file_fingerprint(project, MOD_GWAS, NA), error = function(e) "err")
    paste(fp_gea, fp_gwas, sep = "__")
}

#' Compute significant SNPs interactively using a user-chosen threshold.
#'
#' Wraps compute_pval_threshold() (sourced from pipeline utils) per method × trait.
#' Disk-caches results keyed by (module, k, regime, type, value, file_fingerprint)
#' so qval (slow) is paid once per unique combination; bonf/top/custom are sub-second.
#' Including the file fingerprint ensures the cache is automatically bypassed when the
#' pipeline regenerates p-value files with a different trait set.
#'
#' @param pvalues_list Named list of wide data.tables (method → dt) from
#'   load_all_method_pvalues() or load_all_method_wza_pvalues().
#' @param type   Character: "bonf" | "qval" | "top" | "custom"
#' @param value  Numeric threshold value (alpha for bonf/qval, N for top, raw p for custom)
#' @param k      K-best value (for cache key)
#' @param regime Character: "snp" or "wza"
#' @param project Project name (for cache dir)
#' @param module  Pipeline module (for cache dir)
#' @param overrides Named list: "trait::method" -> list(type=, value=) for
#'   cells that deviate from the fallback (registry-or-master). A cell absent
#'   from this map follows registry_defaults (if the method's family is not
#'   univariate) or master. Default list() is byte-identical to the
#'   pre-rework single-global behaviour (see R/fct_threshold_rules.R).
#' @param registry_defaults Named list: method -> list(adjust=,threshold=,family=),
#'   from gea_method_significance_defaults(). Pins non-univariate methods
#'   (RDA) to their registry rule unless a cell override exists.
#' @return Named list (method → data.table with SNPID/chr/pos/pvalue/method/trait),
#'   same shape as load_all_method_sigsnps().
#' @noRd
compute_method_sigsnps_cached <- function(pvalues_list, type, value,
                                           k, regime, project, module,
                                           cutoffs = NULL, overrides = list(),
                                           registry_defaults = list()) {
    if (length(pvalues_list) == 0 || is.null(pvalues_list)) return(list())

    # Validate that compute_pval_threshold is available (sourced above)
    if (!exists("compute_pval_threshold", mode = "function")) {
        warning("compute_pval_threshold not available — returning empty sig set")
        return(list())
    }

    # Fingerprint of the source pvalue files (md5 of mtime+size) — ensures a
    # pipeline re-run (new trait set) bypasses stale cache entries.
    fp        <- pvalues_file_fingerprint(project, module, k, regime)
    cache_dir <- interactive_sigsnps_dir(project, module)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    # In-memory front for the disk cache — with per-cell (method, trait) disk
    # files, a naive read-every-time would mean up to n_methods * n_traits
    # fread() calls on every recompute of this reactive; get_app_cache() keeps
    # the hot set in memory across renders within the session.
    mem_cache <- get_app_cache()

    all_rows <- list()
    for (m in names(pvalues_list)) {
        pv_dt <- pvalues_list[[m]]
        if (nrow(pv_dt) == 0) next

        fixed_cols <- c("SNPID", "chr", "pos")
        # WZA tables also have n_snps, mean_maf — skip non-trait cols
        trait_cols <- setdiff(names(pv_dt), c(fixed_cols, "n_snps", "mean_maf"))
        if (length(trait_cols) == 0) next

        for (tr in trait_cols) {
            # Resolved rule for THIS (method, trait) cell — a cell override
            # wins, else the method's registry rule (RDA), else master.
            rule <- effective_rule_for(method = m, trait = tr, overrides = overrides,
                                       master_type = type, master_value = value,
                                       registry_defaults = registry_defaults)

            # Disk-cache key is scoped to exactly this cell's resolved rule —
            # editing one cell's threshold recomputes exactly this one
            # (method, trait), not the whole project's sig-SNP set.
            # `rule_version` invalidates every entry written under an older
            # selection rule. Bump it whenever the mask SEMANTICS change: none of
            # the other key components move when only the comparison operator
            # does, so without it a warm project keeps serving a set computed
            # under the old rule — from disk and from mem_cache alike, since
            # mem_key is derived from this same hash.
            #   v2: exclusive '<' -> inclusive '<=' (boundary SNP + ties dropped)
            cache_key  <- list(module = module, k = k, regime = regime,
                               method = m, trait = tr,
                               rule_type = rule$type, rule_value = rule$value, fp = fp,
                               rule_version = 2L)
            cache_hash <- digest::digest(cache_key, algo = "md5")
            mem_key    <- paste0("cellsig_", cache_hash)

            cell_dt <- mem_cache$get(mem_key)
            if (cachem::is.key_missing(cell_dt)) {
                cache_file <- file.path(cache_dir, paste0(cache_hash, ".tsv"))
                cell_dt <- if (file.exists(cache_file)) {
                    tryCatch(
                        data.table::fread(cache_file, sep = "\t", header = TRUE,
                                          colClasses = c(chr = "character", SNPID = "character",
                                                         method = "character", trait = "character")),
                        error = function(e) NULL
                    )
                } else NULL

                if (is.null(cell_dt)) {
                    pvec <- pv_dt[[tr]]

                    # Use precomputed cutoff from combo_thresholds if available
                    # (avoids a redundant full-vector scan on cold cache; same
                    # value, different path) — only valid when it was computed
                    # under the SAME resolved rule as this cell's.
                    precomp_key <- paste0(tr, "::", m)
                    precomp_thr <- if (!is.null(cutoffs)) cutoffs[[precomp_key]] else NULL
                    thresh_result <- if (!is.null(precomp_thr) && !is.na(precomp_thr) && precomp_thr > 0) {
                        list(status = "ok", threshold = precomp_thr)
                    } else {
                        tryCatch(
                            compute_pval_threshold(pvec, rule$type, rule$value),
                            error = function(e) list(status = "error", threshold = NA_real_)
                        )
                    }

                    cell_dt <- if (thresh_result$status != "ok" || is.na(thresh_result$threshold)) {
                        data.table::data.table(SNPID = character(), chr = character(), pos = integer(),
                                               pvalue = numeric(), method = character(), trait = character())
                    } else {
                        # Inclusive '<=', same rule as the pipeline's sig_snps.R mask —
                        # compute_pval_threshold() returns a cutoff that is a member
                        # of the call set under 'qval'/'top'. A strict '<' made the
                        # app's interactive set one SNP (plus ties) short of the
                        # pipeline TSVs.
                        mask <- pvec <= thresh_result$threshold & !is.na(pvec)
                        if (!any(mask)) {
                            data.table::data.table(SNPID = character(), chr = character(), pos = integer(),
                                                   pvalue = numeric(), method = character(), trait = character())
                        } else {
                            pv_dt[mask, .(
                                SNPID  = SNPID, chr = chr, pos = pos,
                                pvalue = .SD[[tr]], method = m, trait = tr
                            ), .SDcols = tr]
                        }
                    }
                    data.table::fwrite(cell_dt, cache_file, sep = "\t")
                }
                mem_cache$set(mem_key, cell_dt)
            }

            if (nrow(cell_dt) > 0) all_rows[[length(all_rows) + 1L]] <- cell_dt
        }
    }

    if (length(all_rows) == 0) return(list())
    flat <- data.table::rbindlist(all_rows)
    if (nrow(flat) == 0) return(list())
    split(flat, by = "method", keep.by = TRUE)
}

#' Compute significance threshold per trait × method (no IO, no cache)
#'
#' Returns a named numeric vector keyed "trait::method" -> threshold p-value.
#' NA when too few tests or compute_pval_threshold status != "ok".
#' Works for all modes:
#'   bonf -> 0.05 / n_tested
#'   qval -> max p still passing FDR
#'   top  -> Nth smallest p
#' For WZA regime, pass the WZA window pvalue list — denominator is n_windows.
#'
#' @param pvalues_list Named list: method -> data.table with SNPID/chr/pos + trait cols
#' @param type         Adjustment type: "bonf", "qval", "top"
#' @param value        Significance level / FDR / top-N (numeric or character)
#' @param overrides    Named list: "trait::method" -> list(type=, value=) for
#'   cells that deviate from the fallback. Default list() = identical output
#'   to before per-cell rules existed — protects any caller that doesn't pass
#'   the new argument.
#' @param registry_defaults Named list: method -> list(adjust=,threshold=,family=).
#'   Pins non-univariate methods (RDA) to their registry rule unless a cell
#'   override exists. Default list() = no method is pinned.
#' @return Named numeric vector "trait::method" -> threshold (NA if unavailable)
#' @noRd
compute_method_thresholds <- function(pvalues_list, type, value, overrides = list(),
                                      registry_defaults = list()) {
    if (length(pvalues_list) == 0 || is.null(pvalues_list)) return(list())
    fixed_cols <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
    out <- list()
    for (m in names(pvalues_list)) {
        pv_dt <- pvalues_list[[m]]
        if (is.null(pv_dt) || nrow(pv_dt) == 0) next
        trait_cols <- setdiff(names(pv_dt), fixed_cols)
        if (length(trait_cols) == 0) next
        for (tr in trait_cols) {
            rule <- effective_rule_for(method = m, trait = tr, overrides = overrides,
                                       master_type = type, master_value = value,
                                       registry_defaults = registry_defaults)
            pvec <- pv_dt[[tr]]
            result <- tryCatch(
                compute_pval_threshold(pvec, rule$type, rule$value),
                error = function(e) list(status = "error", threshold = NA_real_)
            )
            key <- paste0(tr, "::", m)
            out[[key]] <- if (identical(result$status, "ok")) result$threshold else NA_real_
        }
    }
    out
}

#' Cached data loader using cachem
#' Creates a session-level cache automatically on first call.
#' @noRd
.cache_env <- new.env(parent = emptyenv())

get_app_cache <- function() {
    if (is.null(.cache_env$cache)) {
        .cache_env$cache <- cachem::cache_mem(max_size = 200 * 1024^2)  # 200 MB
    }
    .cache_env$cache
}

#' Load with session cache (avoids re-reading the same file across tab switches)
#'
#' @param cache_key  Base cache key string.
#' @param loader_fn  Zero-argument function that loads the data.
#' @param fingerprint Optional string appended to the key. Pass a file's mtime or
#'   a digest of mtime+size so the cache is bypassed when the file is regenerated
#'   by the pipeline without changing its path.
#' @noRd
load_cached <- function(cache_key, loader_fn, fingerprint = NULL) {
    cache <- get_app_cache()
    full_key <- if (!is.null(fingerprint)) paste0(cache_key, "_", fingerprint) else cache_key
    safe_key <- gsub("[^a-z0-9]", "", tolower(full_key))
    cached <- cache$get(safe_key)
    if (!cachem::is.key_missing(cached)) return(cached)
    val <- loader_fn()
    cache$set(safe_key, val)
    val
}
