#!/usr/bin/env Rscript
# Write pipeline summary table for the current mode
# Appends/updates a running Pipeline_summary.tsv with metrics from the current step

library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
################
MODE = args[1]           # Pipeline mode: processing, prestructure, climate, traits, pregea, structure, gea, gwas, gea_x_gwas, maladaptation
OUTPUT = args[2]         # Pipeline_summary.tsv path
# Remaining args are mode-specific input files
################

message(paste0('INFO: Writing summary for mode: ', MODE))

######################## Functions

# Read existing summary or create empty one
read_or_create_summary <- function(output_path) {
    if (file.exists(output_path) && file.size(output_path) > 0) {
        existing <- fread(output_path)
        message(paste0('INFO: Loaded existing summary with ', nrow(existing), ' rows'))
        return(existing)
    }
    return(data.table(step = character(), metric = character(), value = character()))
}

# Add rows to summary, replacing any existing rows for the same step
update_summary <- function(existing, new_rows) {
    # Remove old rows for this step
    step_name <- unique(new_rows$step)
    existing <- existing[!step %in% step_name]
    # Append new rows
    rbind(existing, new_rows)
}

# Create a row
row <- function(step, metric, value) {
    data.table(step = step, metric = metric, value = as.character(value))
}

#################### Main

# Read existing summary
summary_dt <- read_or_create_summary(OUTPUT)

if (MODE == 'processing') {
    # args: MODE OUTPUT vcf_filt vcf_ld samples_list samples_filtered samples_removed
    #       qc_raw_summary filtering_summary het_outlier_sd has_dp depth_summary
    #       pihat_thresh relatedness_removed relatedness_action
    VCF_FILT          = args[3]
    VCF_LD            = args[4]
    SAMPLES_LIST      = args[5]
    SAMPLES_FILTERED  = args[6]
    SAMPLES_REMOVED   = args[7]
    QC_RAW_SUMMARY    = args[8]
    FILTERING_SUMMARY = args[9]
    HET_OUTLIER_SD    = args[10]
    HAS_DP            = args[11]
    DEPTH_SUMMARY     = if (length(args) >= 12) args[12] else 'NULL'
    PI_HAT_THRESH     = if (length(args) >= 13) args[13] else 'NULL'
    RELATEDNESS_REMOVED = if (length(args) >= 14) args[14] else 'NULL'
    RELATEDNESS_ACTION  = if (length(args) >= 15) args[15] else 'keep'
    COORD_MISSING_SUMMARY = if (length(args) >= 16) args[16] else 'NULL'

    # Count samples
    n_samples_total    <- length(readLines(SAMPLES_LIST))
    n_samples_filtered <- length(readLines(SAMPLES_FILTERED))
    n_samples_removed  <- length(readLines(SAMPLES_REMOVED))

    # Count SNPs from VCF files (non-header lines)
    count_vcf_snps <- function(vcf_path) {
        con <- file(vcf_path, 'r')
        n <- 0L
        while (length(line <- readLines(con, n = 10000)) > 0)
            n <- n + sum(!startsWith(line, '#'))
        close(con)
        return(n)
    }
    n_snps_filtered <- count_vcf_snps(VCF_FILT)
    n_snps_ld       <- count_vcf_snps(VCF_LD)

    # Parse raw bcftools stats output (SN/TSTV lines)
    parse_bcftools_stats <- function(path) {
        if (!file.exists(path) || file.size(path) == 0) return(list(snps=NA, titv=NA))
        lines <- readLines(path)
        snps <- NA; titv <- NA
        for (ln in lines) {
            if (startsWith(ln, 'SN')) {
                parts <- strsplit(ln, '\t')[[1]]
                if (length(parts) >= 4 && grepl('number of SNPs', parts[3]))
                    snps <- trimws(parts[4])
            } else if (startsWith(ln, 'TSTV')) {
                parts <- strsplit(ln, '\t')[[1]]
                if (length(parts) >= 5) titv <- trimws(parts[5])
            }
        }
        list(snps=snps, titv=titv)
    }
    bcf <- parse_bcftools_stats(QC_RAW_SUMMARY)
    n_raw_snps <- bcf$snps
    titv       <- bcf$titv

    # Read filtering summary for detailed attrition
    filt_sum <- fread(FILTERING_SUMMARY)

    # Het outlier count from filtering summary
    n_het_removed <- 0L
    if (HET_OUTLIER_SD != 'NULL' && 'After het outlier removal' %in% filt_sum$stage) {
        n_before <- filt_sum[stage == 'After sample missingness filter', n_samples]
        n_after  <- filt_sum[stage == 'After het outlier removal', n_samples]
        if (length(n_before) > 0 && length(n_after) > 0)
            n_het_removed <- as.integer(n_before) - as.integer(n_after)
    }

    # Depth info — three states:
    #   not_provided:   VCF has no FORMAT/DP field
    #   declared_empty: FORMAT/DP declared but genotype DP values were empty/non-finite
    #                   (typical for GBS) — plot_qc_processing.R skips depth_summary.tsv in this case
    #   available:      usable depth values found, depth_summary.tsv was written
    depth_summary_exists <- HAS_DP == 'TRUE' && DEPTH_SUMMARY != 'NULL' && file.exists(DEPTH_SUMMARY)
    depth_qc_status <- if (HAS_DP != 'TRUE') {
        'not_provided'
    } else if (depth_summary_exists) {
        'available'
    } else {
        'declared_empty'
    }

    # Relatedness (IBS allele-sharing) removal count — 0 unless Filter.relatedness_action is 'remove'.
    # relatedness_removed.tsv is always computed once a threshold is set (even in 'keep' mode,
    # for the Shiny hover-note preview), so gate on the action here to avoid reporting samples
    # as "removed" when they were only previewed.
    n_relatedness_removed <- 0L
    if (tolower(RELATEDNESS_ACTION) == 'remove' &&
        PI_HAT_THRESH != 'NULL' && RELATEDNESS_REMOVED != 'NULL' && file.exists(RELATEDNESS_REMOVED)) {
        rel_removed <- fread(RELATEDNESS_REMOVED)
        n_relatedness_removed <- nrow(rel_removed)
    }

    new_rows <- rbind(
        row('processing', 'samples_total',              n_samples_total),
        row('processing', 'samples_after_filtering',    n_samples_filtered),
        row('processing', 'samples_removed',            n_samples_removed),
        row('processing', 'samples_het_outliers_removed', n_het_removed),
        row('processing', 'samples_removed_relatedness', n_relatedness_removed),
        row('processing', 'snps_raw',                   n_raw_snps),
        row('processing', 'snps_after_filtering',       n_snps_filtered),
        row('processing', 'snps_after_ld_pruning',      n_snps_ld),
        row('processing', 'titv_ratio',                 titv),
        row('processing', 'depth_qc',                   depth_qc_status)
    )

    # Append depth stats if available
    if (depth_summary_exists) {
        depth_dt <- fread(DEPTH_SUMMARY)
        for (i in seq_len(nrow(depth_dt))) {
            new_rows <- rbind(new_rows, row('processing', depth_dt$metric[i], depth_dt$value[i]))
        }
    }

    # Coordinate availability (climate/GEA auto-DROP of samples with missing lat/lon)
    if (COORD_MISSING_SUMMARY != 'NULL' && file.exists(COORD_MISSING_SUMMARY)) {
        coord_sum <- fread(COORD_MISSING_SUMMARY)
        new_rows <- rbind(new_rows,
            row('processing', 'samples_with_coordinates',        coord_sum$n_available[1]),
            row('processing', 'samples_dropped_missing_coordinates', coord_sum$n_missing[1])
        )
    }

    # The sample accounting must close: total == after + removed + het + related. It cannot
    # when plink --keep silently dropped IDs it could not find (audit 2026-09-13 B41) —
    # samples_total is counted from the keep-list, the other three from the .imiss. Only this
    # one identity is promoted from the opt-in checker; the rest of
    # check_summary_accounting() stays there (climate_predictor_count_disagrees is a filed,
    # known-red defect of this script's structure branch).
    source('/pipeline/scripts/R/lib/invariants.R')
    acct <- check_summary_accounting(new_rows)
    acct <- acct[check == 'sample_accounting_does_not_close']
    if (nrow(acct) > 0) {
        stop(paste0('Sample accounting does not close: ', acct$detail[1],
                    '. A keep-list ID missing from the VCF header is the usual cause — ',
                    'see calculate_sample_missing.log.'))
    }

} else if (MODE == 'prestructure') {
    # args: MODE OUTPUT cross_entropy_plot K_START K_END
    K_START = args[3] %>% as.numeric
    K_END = args[4] %>% as.numeric

    new_rows <- rbind(
        row('prestructure',  'K_range', paste0(K_START, '-', K_END))
    )

} else if (MODE == 'pregea') {
    # args: MODE OUTPUT recs_path lfmm_path emmax_path rda_path collin_path guard_path
    # varpart/dbMEM metrics moved to the 'climate' branch below (module split —
    # PreGEA no longer produces those artifacts). guard_path may be 'NULL'
    # when PreGEA.TransferGuard.enabled is off.
    RECS_PATH    = args[3]
    LFMM_PATH    = args[4]
    EMMAX_PATH   = args[5]
    RDA_PATH     = args[6]
    COLLIN_PATH  = args[7]
    GUARD_PATH   = args[8]
    # DEFLATION_FLOOR mirrors PreGEA's fixed constant (common.smk PREGEA_DEFLATION_FLOOR,
    # not user-configurable) — a rung counts as "deflated" only below the floor
    # the recommender itself uses, not against a raw lambda_gc < 1 cutoff (that
    # flagged even the recommended EMMAX rung, e.g. lambda_gc=0.93 > 0.90 floor).
    DEFLATION_FLOOR = 0.90

    read_opt <- function(path) if (identical(path, 'NULL') || !file.exists(path)) NULL else fread(path, sep = '\t', header = TRUE)

    new_rows <- list()
    add <- function(metric, value) new_rows[[length(new_rows) + 1]] <<- row('pregea', metric, value)

    recs <- read_opt(RECS_PATH)
    if (!is.null(recs) && nrow(recs) > 0) {
        for (m in c('LFMM', 'EMMAX', 'RDA')) {
            mrows <- recs[method == m]
            for (i in seq_len(nrow(mrows))) {
                add(sprintf('%s_%s_recommended', m, mrows$param[i]), mrows$recommended_value[i])
                add(sprintf('%s_%s_evidence', m, mrows$param[i]),
                   sprintf('%s=%s', mrows$evidence_metric[i], mrows$evidence_value[i]))
            }
        }
    }

    lfmm <- read_opt(LFMM_PATH)
    if (!is.null(lfmm)) {
        pooled <- lfmm[trait == '__pooled__']
        add('LFMM_K_range', paste(sort(unique(pooled$rung_value_num)), collapse = ','))
        if (nrow(pooled) > 0) add('LFMM_hist_flatness_best', min(pooled$hist_flatness_ks, na.rm = TRUE))
    }

    emmax <- read_opt(EMMAX_PATH)
    if (!is.null(emmax)) {
        pooled <- emmax[trait == '__pooled__']
        add('EMMAX_n_pcs_range', paste(sort(unique(pooled$rung_value_num)), collapse = ','))
        z <- pooled[rung_value_num == 0]
        if (nrow(z) > 0) add('EMMAX_lambda_at_zero_pcs', z$lambda_gc[1])
        add('EMMAX_n_rungs_deflated', sum(pooled$lambda_gc < DEFLATION_FLOOR, na.rm = TRUE))
    }

    rda <- read_opt(RDA_PATH)
    if (!is.null(rda)) {
        ok <- rda[status == 'ok']
        if (nrow(ok) > 0) add('RDA_condition_pcs_range', paste(sort(unique(ok$condition_pcs)), collapse = ','))
    }
    collin <- read_opt(COLLIN_PATH)
    if (!is.null(collin)) {
        add('RDA_predictors_input', nrow(collin))
        add('RDA_predictors_retained', sum(collin$action == 'kept'))
        add('RDA_predictors_dropped_collinear', sum(collin$action == 'dropped'))
    }

    guard <- read_opt(GUARD_PATH)
    if (!is.null(guard) && nrow(guard) > 0) {
        add('transfer_guard_status', 'ok')
    } else {
        add('transfer_guard_status', 'disabled')
    }

    new_rows <- if (length(new_rows) > 0) rbindlist(new_rows) else data.table(step = character(), metric = character(), value = character())

} else if (MODE == 'climate') {
    # args: MODE OUTPUT invariant_path dbmem_diag_path variance_partition_path confounding_path px_path
    INVARIANT_PATH = args[3]
    DBMEM_PATH     = args[4]
    VARPART_PATH   = args[5]
    CONFOUND_PATH  = args[6]
    PX_PATH        = args[7]

    read_opt <- function(path) if (identical(path, 'NULL') || !file.exists(path)) NULL else fread(path, sep = '\t', header = TRUE)

    # Human-readable component names (e.g. "Climate", "Climate ∩ Structure")
    # -> a clean summary-key suffix (climate, climate_structure). Same
    # transform for every branch's rows now that scripts/pregea_varpart.R
    # writes one uniform schema regardless of which model (3-way/2-way/
    # climate-only) actually ran — no more hardcoded fraction-name list.
    comp_key <- function(comp) {
        k <- gsub("[^A-Za-z0-9]+", "_", comp)
        k <- gsub("^_+|_+$", "", k)
        tolower(k)
    }

    new_rows <- list()
    add <- function(metric, value) new_rows[[length(new_rows) + 1]] <<- row('climate', metric, value)

    invariant <- read_opt(INVARIANT_PATH)
    if (!is.null(invariant)) add('n_invariant_predictors', nrow(invariant))

    dbmem <- read_opt(DBMEM_PATH)
    if (!is.null(dbmem)) {
        dv <- setNames(dbmem$value, dbmem$key)
        for (k in c('spatial_level', 'n_sites', 'n_unique_coords', 'mst_threshold_km', 'n_mem_positive', 'status'))
            if (k %in% names(dv)) add(paste0('dbMEM_', k), dv[[k]])
    }

    varpart <- read_opt(VARPART_PATH)
    if (!is.null(varpart) && nrow(varpart) > 0) {
        add('varpart_status', unique(varpart$status)[1])
        if ('model' %in% names(varpart)) add('varpart_model', unique(varpart$model)[1])
        for (i in seq_len(nrow(varpart))) {
            add(paste0('varpart_', comp_key(varpart$component[i]), '_R2adj'), varpart$variance_pct[i])
        }
    }

    confound <- read_opt(CONFOUND_PATH)
    if (!is.null(confound) && nrow(confound) > 0) {
        add('varpart_confounding_flag', as.logical(confound$confounded[1]))
    }

    px <- read_opt(PX_PATH)
    if (!is.null(px) && nrow(px) > 0) {
        top <- px[order(-Px)][1]
        add('Px_top_variable', top$variable)
        add('Px_top_value', top$Px)
    }

    new_rows <- if (length(new_rows) > 0) rbindlist(new_rows) else data.table(step = character(), metric = character(), value = character())

} else if (MODE == 'traits') {
    # args: MODE OUTPUT trait_summary_path trait_invariant_path pairs_png_path corr_climate_path
    SUMMARY_PATH   = args[3]
    INVARIANT_PATH = args[4]
    PAIRS_PATH     = args[5]
    CORR_CLIM_PATH = if (length(args) >= 6) args[6] else "NULL"

    read_opt <- function(path) if (identical(path, 'NULL') || !file.exists(path)) NULL else fread(path, sep = '\t', header = TRUE)

    new_rows <- list()
    add <- function(metric, value) new_rows[[length(new_rows) + 1]] <<- row('traits', metric, value)

    trait_summary <- read_opt(SUMMARY_PATH)
    if (!is.null(trait_summary)) {
        add('n_traits', nrow(trait_summary))
        if (nrow(trait_summary) > 0) {
            add('traits', paste(trait_summary$trait, collapse = ','))
            if ('n_missing' %in% names(trait_summary))
                add('n_traits_with_missing', sum(trait_summary$n_missing > 0))
        }
    }

    invariant <- read_opt(INVARIANT_PATH)
    if (!is.null(invariant)) {
        add('n_invariant_traits', nrow(invariant))
        if (nrow(invariant) > 0) add('invariant_traits', paste(invariant$predictor, collapse = ','))
    }

    # Whether the pairs grid was rendered or replaced by the >max_factors
    # placeholder is not recorded in any table — the plot script logs it. What
    # IS reportable here is simply that both figures exist.
    add('pairs_plot_written', file.exists(PAIRS_PATH))
    add('climate_correlogram_written', !identical(CORR_CLIM_PATH, 'NULL') && file.exists(CORR_CLIM_PATH))

    new_rows <- if (length(new_rows) > 0) rbindlist(new_rows) else data.table(step = character(), metric = character(), value = character())

} else if (MODE == 'structure') {
    # args: MODE OUTPUT K_BEST climate_site predictors ld_decay_path ld_decay_group_by ld_decay_scope
    #       climate_na_excluded
    K_BEST = args[3]
    CLIMATE_SITE = args[4]
    PREDICTORS = args[5]
    LD_DECAY_PATH = if (length(args) >= 6) args[6] else "NULL"
    LD_DECAY_GROUP_BY = if (length(args) >= 7) args[7] else "NULL"
    LD_DECAY_SCOPE_VAL = if (length(args) >= 8) args[8] else "NULL"
    CLIMATE_NA_EXCLUDED = if (length(args) >= 9) args[9] else "NULL"

    new_rows <- rbind(row('structure', 'K_best', K_BEST))

    if (CLIMATE_SITE != 'NULL') {
        predictors_list <- str_split(PREDICTORS, ',')[[1]]
        climate <- fread(CLIMATE_SITE)
        n_climate_vars <- ncol(climate) - 4  # Subtract site, sample, lat, lon columns

        new_rows <- rbind(new_rows,
            row('structure', 'climate_predictors', PREDICTORS),
            row('structure', 'n_climate_variables', n_climate_vars)
        )
    }

    # Climate-NA exclusions (samples whose raster extraction returned NA, e.g. ocean/NoData
    # pixel) -- mirrors the processing-mode coord-missing summary. Zero when no sample hits an
    # ocean/NoData pixel, the common case.
    if (CLIMATE_NA_EXCLUDED != 'NULL' && file.exists(CLIMATE_NA_EXCLUDED)) {
        na_excl <- fread(CLIMATE_NA_EXCLUDED)
        new_rows <- rbind(new_rows, row('structure', 'samples_excluded_climate_na', nrow(na_excl)))
        if (nrow(na_excl) > 0) {
            new_rows <- rbind(new_rows,
                row('structure', 'climate_na_excluded_samples', paste(na_excl$sample, collapse = ',')))
        }
    }

    # LD decay summary
    if (LD_DECAY_PATH != 'NULL' && file.exists(LD_DECAY_PATH)) {
        ld_table <- fread(LD_DECAY_PATH)
        gw_all <- ld_table[group == 'All' & scope == 'genome_wide']
        if (nrow(gw_all) > 0) {
            new_rows <- rbind(new_rows,
                row('structure', 'ld_decay_half_decay_genome_wide', gw_all$half_decay_bp[1]),
                row('structure', 'ld_decay_r2_02_genome_wide', gw_all$r2_02_bp[1])
            )
        }
        new_rows <- rbind(new_rows,
            row('structure', 'ld_decay_group_by', LD_DECAY_GROUP_BY),
            row('structure', 'ld_decay_groups_analyzed', nrow(ld_table[scope == 'genome_wide'])),
            row('structure', 'ld_decay_scope', LD_DECAY_SCOPE_VAL)
        )
    }

} else if (MODE %in% c('gea', 'gwas')) {
    # args: MODE OUTPUT selected_snps regions_per_trait regions_combined genes_per_region
    #        [region_dist_spec r2_threshold ld_decay_group ld_decay_path]
    is_gwas <- MODE == 'gwas'
    arg_offset <- if (is_gwas) 1L else 0L  # gwas has an extra missing_summary arg before selected_snps
    MISSING_SUMMARY   <- if (is_gwas) args[3] else NULL
    SELECTED_SNPS     <- args[3 + arg_offset]
    REGIONS_PER_TRAIT <- args[4 + arg_offset]
    REGIONS_COMBINED  <- args[5 + arg_offset]
    GENES_PER_REGION  <- args[6 + arg_offset]
    REGION_DIST_SPEC  <- if (length(args) >= 7 + arg_offset) args[7 + arg_offset] else 'unknown'
    R2_THRESHOLD_VAL  <- if (length(args) >= 8 + arg_offset) as.numeric(args[8 + arg_offset]) else 0.2
    LD_DECAY_GROUP_VAL <- if (length(args) >= 9 + arg_offset) args[9 + arg_offset] else 'All'
    LD_DECAY_PATH_VAL <- if (length(args) >= 10 + arg_offset) args[10 + arg_offset] else 'NULL'

    # Count selected SNPs
    snps <- fread(SELECTED_SNPS)
    n_selected_snps <- nrow(snps)

    # Count per-method SNPs
    method_cols <- setdiff(colnames(snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
    method_rows <- list()
    for (m in method_cols) {
        n_method <- sum(!is.na(snps[[m]]) & snps[[m]] != '' & snps[[m]] != '""')
        method_rows <- c(method_rows, list(row(MODE, paste0('sig_snps_', m), n_method)))
    }

    # Count regions
    regions_trait    <- fread(REGIONS_PER_TRAIT)
    regions_combined <- fread(REGIONS_COMBINED)

    # Count genes
    genes   <- fread(GENES_PER_REGION)
    n_genes <- if (nrow(genes) > 0 && 'gene_id' %in% colnames(genes)) n_distinct(genes$gene_id) else 0

    new_rows <- rbind(
        row(MODE, 'selected_snps_total', n_selected_snps),
        do.call(rbind, method_rows),
        row(MODE, 'regions_per_trait', nrow(regions_trait)),
        row(MODE, 'regions_combined',  nrow(regions_combined)),
        row(MODE, 'genes_found',       n_genes),
        row(MODE, 'region_distance_spec',    REGION_DIST_SPEC),
        row(MODE, 'region_r2_threshold',     R2_THRESHOLD_VAL),
        row(MODE, 'region_ld_decay_group',   LD_DECAY_GROUP_VAL)
    )

    # Per-chromosome distances from LD decay table (when auto mode)
    if (LD_DECAY_PATH_VAL != 'NULL' && file.exists(LD_DECAY_PATH_VAL) &&
        REGION_DIST_SPEC %in% c('auto_per_chromosome', 'auto_genome_wide', 'auto')) {
        ld_table <- fread(LD_DECAY_PATH_VAL)
        if (REGION_DIST_SPEC == 'auto_per_chromosome') {
            chr_rows <- ld_table[group == LD_DECAY_GROUP_VAL & !scope %in% c('genome_wide')]
            if (nrow(chr_rows) > 0) {
                new_rows <- rbind(new_rows,
                    row(MODE, 'region_distance_chr_min_bp',
                        round(min(chr_rows$r2_02_bp, na.rm = TRUE))),
                    row(MODE, 'region_distance_chr_median_bp',
                        round(median(chr_rows$r2_02_bp, na.rm = TRUE))),
                    row(MODE, 'region_distance_chr_max_bp',
                        round(max(chr_rows$r2_02_bp, na.rm = TRUE)))
                )
            }
        } else {
            gw_row <- ld_table[group == LD_DECAY_GROUP_VAL & scope == 'genome_wide']
            if (nrow(gw_row) > 0) {
                new_rows <- rbind(new_rows,
                    row(MODE, 'region_distance_genome_wide_bp', round(gw_row$r2_02_bp[1])))
            }
        }
    }

    # Add per-trait region counts
    if (nrow(regions_trait) > 0 && 'trait' %in% colnames(regions_trait)) {
        trait_counts <- regions_trait %>%
            dplyr::group_by(trait) %>%
            dplyr::summarise(n = dplyr::n(), .groups = 'drop')
        for (i in seq_len(nrow(trait_counts))) {
            new_rows <- rbind(new_rows,
                row(MODE, paste0('regions_', trait_counts$trait[i]), trait_counts$n[i]))
        }
    }

    # GWAS-specific: missing value summary
    if (is_gwas && !is.null(MISSING_SUMMARY) && file.exists(MISSING_SUMMARY)) {
        missing <- fread(MISSING_SUMMARY)
        new_rows <- rbind(new_rows,
            row('gwas', 'n_phenotype_traits', nrow(missing)),
            row('gwas', 'missing_strategy',   missing$strategy[1])
        )
        for (i in seq_len(nrow(missing))) {
            t <- missing$trait[i]
            new_rows <- rbind(new_rows,
                row('gwas', paste0(t, '_n_available'), missing$n_available[i]),
                row('gwas', paste0(t, '_n_missing'),   missing$n_missing[i]))
            if (missing$n_missing[i] > 0 && nchar(missing$missing_samples[i]) > 0) {
                new_rows <- rbind(new_rows,
                    row('gwas', paste0(t, '_dropped_samples'), missing$missing_samples[i]))
            }
        }
    }

} else if (MODE == 'maladaptation') {
    # args: MODE OUTPUT offset_site_values
    #
    # offset_site_values is a MANIFEST FILE (one path per line, *.txt) written by
    # the write_summary rule. It used to be a space-joined string, but once offsets
    # gained the scenario dimension the list grew to methods x sets x scenarios --
    # 1176 paths on a 42-scenario sweep -- and overflowed the command line. The
    # space-joined form is still accepted so a hand-run call keeps working.
    OFFSET_SITE = args[3]

    offset_paths <- if (grepl('\\.txt$', OFFSET_SITE) && file.exists(OFFSET_SITE)) {
        readLines(OFFSET_SITE)
    } else {
        strsplit(trimws(OFFSET_SITE), "\\s+")[[1]]
    }
    offset_paths <- trimws(offset_paths)
    offset_paths <- offset_paths[nchar(offset_paths) > 0]

    # Emit per-tag offset stats (one row per method x SNP set x spatial tag x scenario)
    # Path: .../tables/{method}/{run_label}_{spatial_tag}/{scenario}/genetic_offset_site.tsv
    # tag = "{method}_{run_label}_{spatial_tag}_{scenario}"
    tag_rows <- lapply(offset_paths, function(p) {
        scenario_dir <- basename(dirname(p))                    # e.g. "deme_004"
        spatial_dir  <- basename(dirname(dirname(p)))           # e.g. "EMMAX_bonf005_nospatial"
        method_dir   <- basename(dirname(dirname(dirname(p))))  # e.g. "geometric_offset"
        tag <- paste0(method_dir, '_', spatial_dir, '_', scenario_dir)
        if (!file.exists(p)) {
            return(list(
                row('maladaptation', sprintf('offset_min_%s', tag), NA),
                row('maladaptation', sprintf('offset_max_%s', tag), NA),
                row('maladaptation', sprintf('offset_mean_%s', tag), NA)
            ))
        }
        offset <- fread(p, colClasses = c(site = 'character', sample = 'character'))
        offset_col <- setdiff(colnames(offset), c('site', 'sample', 'latitude', 'longitude'))
        if (length(offset_col) > 0) {
            vals <- offset[[offset_col[1]]]
            list(
                row('maladaptation', sprintf('offset_min_%s', tag),  round(min(vals,  na.rm = TRUE), 4)),
                row('maladaptation', sprintf('offset_max_%s', tag),  round(max(vals,  na.rm = TRUE), 4)),
                row('maladaptation', sprintf('offset_mean_%s', tag), round(mean(vals, na.rm = TRUE), 4))
            )
        } else {
            list(
                row('maladaptation', sprintf('offset_min_%s', tag),  NA),
                row('maladaptation', sprintf('offset_max_%s', tag),  NA),
                row('maladaptation', sprintf('offset_mean_%s', tag), NA)
            )
        }
    })

    new_rows <- rbind(
        row('maladaptation', 'maladaptation_methods', 'completed'),
        rbindlist(lapply(unlist(tag_rows, recursive = FALSE), as.data.table))
    )

} else if (MODE == 'gea_x_gwas') {
    # args: MODE OUTPUT overlap_summary selected_snps regions_per_trait regions_combined genes_per_region [enrichment_done]
    OVERLAP_SUMMARY = args[3]
    SELECTED_SNPS = args[4]
    REGIONS_PER_TRAIT = args[5]
    REGIONS_COMBINED = args[6]
    GENES_PER_REGION = args[7]
    ENRICHMENT_DONE = if (length(args) >= 8) args[8] else "NULL"

    # Read overlap summary (last rows are SUMMARY rows)
    overlap <- fread(OVERLAP_SUMMARY)
    summary_rows <- overlap[gea_region_id == 'SUMMARY']
    overlap_pairs <- overlap[gea_region_id != 'SUMMARY']

    n_gea <- as.integer(summary_rows[gwas_region_id == 'gea_total', overlap_length])
    n_gwas <- as.integer(summary_rows[gwas_region_id == 'gwas_total', overlap_length])
    n_gea_overlap <- as.integer(summary_rows[gwas_region_id == 'gea_overlapping', overlap_length])
    n_gwas_overlap <- as.integer(summary_rows[gwas_region_id == 'gwas_overlapping', overlap_length])
    n_pairs <- as.integer(summary_rows[gwas_region_id == 'overlap_pairs', overlap_length])

    gea_pct <- summary_rows[gwas_region_id == 'gea_total', overlap_pct]
    gwas_pct <- summary_rows[gwas_region_id == 'gwas_total', overlap_pct]

    # Count merged SNPs and regions
    snps <- fread(SELECTED_SNPS)
    regions_trait <- fread(REGIONS_PER_TRAIT)
    regions_combined <- fread(REGIONS_COMBINED)
    genes <- fread(GENES_PER_REGION)
    n_genes <- if (nrow(genes) > 0 && 'gene_id' %in% colnames(genes)) n_distinct(genes$gene_id) else 0

    enrichment_summary <- ''
    if (ENRICHMENT_DONE != "NULL" && file.exists(ENRICHMENT_DONE)) {
        enrichment_summary <- 'completed'
    }

    new_rows <- rbind(
        row('gea_x_gwas', 'gea_regions', n_gea),
        row('gea_x_gwas', 'gwas_regions', n_gwas),
        row('gea_x_gwas', 'overlap_pairs', n_pairs),
        row('gea_x_gwas', 'gea_overlapping', n_gea_overlap),
        row('gea_x_gwas', 'gwas_overlapping', n_gwas_overlap),
        row('gea_x_gwas', 'gea_overlap_pct', gea_pct),
        row('gea_x_gwas', 'gwas_overlap_pct', gwas_pct),
        row('gea_x_gwas', 'merged_snps_total', nrow(snps)),
        row('gea_x_gwas', 'new_regions_per_trait', nrow(regions_trait)),
        row('gea_x_gwas', 'new_regions_combined', nrow(regions_combined)),
        row('gea_x_gwas', 'genes_found', n_genes),
        row('gea_x_gwas', 'enrichment_status', enrichment_summary)
    )

} else {
    message(paste0('WARNING: Unknown mode "', MODE, '", skipping summary'))
    quit(save = "no", status = 0)
}

# Update summary
summary_dt <- update_summary(summary_dt, new_rows)

# Save
fwrite(summary_dt, OUTPUT, sep = '\t')
message(paste0('INFO: Saved summary with ', nrow(summary_dt), ' rows to ', OUTPUT))
message('INFO: Complete')
