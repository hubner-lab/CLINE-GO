#!/usr/bin/env Rscript
# VCF Processing QC Plots
# Generates all QC plots and derived tables for the processing module.
#
# Args:
#   1  raw_summary        - vcf_raw_summary.tsv (bcftools stats output)
#   2  sample_miss        - samples_missing_stats.tsv (plink .imiss)
#   3  snp_miss_raw       - snp_missingness_raw.lmiss (plink .lmiss, pre-filter)
#   4  maf_raw            - maf_raw.frq (plink --freq, pre-filter)
#   5  maf_filtered       - maf_filtered.frq (plink --freq, post-filter)
#   6  het_raw            - het_raw.het (plink --het, raw VCF)
#   7  het_filtered       - het_filtered.het (plink --het, filtered VCF)
#   8  snp_density_raw    - snp_density_raw.snpden (vcftools --SNPdensity, raw VCF)
#   9  snp_density_filt   - snp_density_filtered.snpden (vcftools --SNPdensity, filtered VCF)
#   10 metadata           - metadata.tsv (aligned, site column)
#   11 vcf_filt           - filtered VCF (for SNP count)
#   12 vcf_ld             - LD-pruned VCF (for SNP count)
#   13 sample_miss_thresh - sample_miss threshold (numeric)
#   14 snp_miss_thresh    - snp_miss threshold (numeric)
#   15 maf_thresh         - MAF threshold (numeric)
#   16 het_outlier_sd     - het outlier SD threshold or 'NULL'
#   17 has_dp             - 'TRUE' or 'FALSE'
#   18 depth_sample       - depth_per_sample.idepth or 'NULL'
#   19 depth_site         - depth_per_site.ldepth.mean or 'NULL'
#   20 plots_dir          - output directory for plots
#   21 tables_dir         - output directory for tables
#   22 relatedness_pairs  - relatedness_pairs.tsv (plink --genome, IID1/IID2/IBS; always produced)
#   23 pihat_thresh       - Filter.relatedness IBS threshold or 'NULL' (colors/counts only; removal is a separate Filter.relatedness_action gate)
#   24 relatedness_removed - relatedness_removed.tsv (preview of samples that WOULD be dropped, IID/IBS/F_MISS) or 'NULL' when no threshold is set
#   25 relatedness_action  - Filter.relatedness_action ('keep' or 'remove') -- whether the preview above was actually applied

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
})

source("/pipeline/scripts/R/utils/theme_clinego.R")

args <- commandArgs(trailingOnly = TRUE)

RAW_SUMMARY       <- args[1]
SAMPLE_MISS_FILE  <- args[2]
SNP_MISS_RAW_FILE <- args[3]
MAF_RAW_FILE      <- args[4]
MAF_FILT_FILE     <- args[5]
HET_RAW_FILE      <- args[6]
HET_FILT_FILE     <- args[7]
SNP_DENS_RAW      <- args[8]
SNP_DENS_FILT     <- args[9]
METADATA_FILE     <- args[10]
VCF_FILT          <- args[11]
VCF_LD            <- args[12]
SAMPLE_MISS_THRESH <- as.numeric(args[13])
SNP_MISS_THRESH    <- as.numeric(args[14])
MAF_THRESH         <- as.numeric(args[15])
HET_OUTLIER_SD     <- if (args[16] == 'NULL') NULL else as.numeric(args[16])
HAS_DP             <- args[17] == 'TRUE'
DEPTH_SAMPLE_FILE  <- args[18]
DEPTH_SITE_FILE    <- args[19]
PLOTS_DIR          <- args[20]
TABLES_DIR         <- args[21]
RELATEDNESS_PAIRS  <- args[22]
PI_HAT_THRESH      <- if (args[23] == 'NULL') NULL else as.numeric(args[23])
RELATEDNESS_REMOVED_FILE <- if (args[24] == 'NULL') NULL else args[24]
RELATEDNESS_ACTION       <- tolower(args[25])

#=============================================================================
# HELPERS
#=============================================================================

save_plot <- function(p, filename, width = 8, height = 5) {
    ggsave(file.path(PLOTS_DIR, filename), p, width = width, height = height, dpi = 150)
    message("INFO: Saved ", filename)
}

count_vcf_snps <- function(vcf_path) {
    con <- file(vcf_path, 'r')
    n <- 0L
    while (length(line <- readLines(con, n = 10000)) > 0)
        n <- n + sum(!startsWith(line, '#'))
    close(con)
    n
}

#=============================================================================
# 1. LOAD RAW SUMMARY (baseline counts + Ti/Tv — raw bcftools stats format)
#=============================================================================
parse_bcftools_stats <- function(path) {
    out <- list(snps = NA_integer_, samples = NA_integer_, titv = NA_real_)
    if (!file.exists(path) || file.size(path) == 0) return(out)
    lines <- readLines(path)
    for (ln in lines) {
        if (startsWith(ln, "SN")) {
            parts <- strsplit(ln, "\t")[[1]]
            if (length(parts) >= 4) {
                if (grepl("number of SNPs", parts[3]))    out$snps    <- as.integer(trimws(parts[4]))
                if (grepl("number of samples", parts[3])) out$samples <- as.integer(trimws(parts[4]))
            }
        } else if (startsWith(ln, "TSTV")) {
            parts <- strsplit(ln, "\t")[[1]]
            if (length(parts) >= 5) out$titv <- as.numeric(trimws(parts[5]))
        }
    }
    out
}
bcf_stats    <- parse_bcftools_stats(RAW_SUMMARY)
n_raw_snps    <- bcf_stats$snps
n_raw_samples <- bcf_stats$samples
titv_ratio    <- bcf_stats$titv

#=============================================================================
# 2. SAMPLE MISSINGNESS DISTRIBUTION
#=============================================================================
smiss <- fread(SAMPLE_MISS_FILE, colClasses = list(character = "IID"))
p_smiss <- ggplot(smiss, aes(x = F_MISS)) +
    geom_histogram(
        aes(fill = F_MISS > SAMPLE_MISS_THRESH),
        binwidth = 0.02, color = "white", linewidth = 0.2
    ) +
    scale_fill_manual(
        values = c("FALSE" = CLINEGO_RETAINED, "TRUE" = CLINEGO_REMOVED),
        labels = c("FALSE" = "Retained", "TRUE" = "Removed"),
        name = NULL
    ) +
    geom_vline(xintercept = SAMPLE_MISS_THRESH, linetype = "dashed", color = CLINEGO_THRESHOLD) +
    labs(
        title = "Per-sample missingness distribution",
        x = "Fraction missing genotypes (F_MISS)",
        y = "Number of samples"
    ) +
    theme_clinego()
save_plot(p_smiss, "sample_missingness_distribution.png")

#=============================================================================
# 3. PER-SNP MISSINGNESS DISTRIBUTION
#=============================================================================
lmiss <- fread(SNP_MISS_RAW_FILE)
p_lmiss <- ggplot(lmiss, aes(x = F_MISS)) +
    geom_histogram(
        aes(fill = F_MISS > SNP_MISS_THRESH),
        binwidth = 0.01, color = "white", linewidth = 0.1
    ) +
    scale_fill_manual(
        values = c("FALSE" = CLINEGO_RETAINED, "TRUE" = CLINEGO_REMOVED),
        labels = c("FALSE" = "Retained", "TRUE" = "Removed"),
        name = NULL
    ) +
    geom_vline(xintercept = SNP_MISS_THRESH, linetype = "dashed", color = CLINEGO_THRESHOLD) +
    labs(
        title = "Per-SNP missingness distribution (pre-filter)",
        x = "Fraction missing genotypes (F_MISS)",
        y = "Number of SNPs"
    ) +
    theme_clinego()
save_plot(p_lmiss, "snp_missingness_distribution.png")

#=============================================================================
# 4. MAF DISTRIBUTION (raw only, with threshold annotation)
#=============================================================================
maf_raw  <- fread(MAF_RAW_FILE)
maf_filt <- fread(MAF_FILT_FILE)  # kept for attrition computation below

n_below_maf <- nrow(maf_raw[MAF < MAF_THRESH])

p_maf <- ggplot(maf_raw, aes(x = MAF)) +
    geom_histogram(
        aes(fill = MAF < MAF_THRESH),
        binwidth = 0.01, color = NA
    ) +
    scale_fill_manual(
        values = c("FALSE" = CLINEGO_RETAINED, "TRUE" = CLINEGO_REMOVED),
        labels = c("FALSE" = "Retained", "TRUE" = "Removed"),
        name = NULL
    ) +
    geom_vline(xintercept = MAF_THRESH, linetype = "dashed", color = CLINEGO_THRESHOLD) +
    labs(
        title = "Minor allele frequency distribution (pre-filter)",
        x = "Minor allele frequency (MAF)",
        y = "Number of SNPs"
    ) +
    theme_clinego()
save_plot(p_maf, "maf_distribution.png")

#=============================================================================
# 5. HET VS MISSINGNESS SCATTER
#=============================================================================
het  <- fread(HET_FILT_FILE, colClasses = list(character = "IID"))
meta <- fread(METADATA_FILE, colClasses = list(character = "sample"))

# Compute observed heterozygosity
het[, obs_het := (N_NM - O_HOM) / N_NM]
het[, sample := IID]

# Merge with metadata to get site (kept for the table only \u2014 NOT used to color
# the plot: population studies typically have far more sampling sites than a
# legend can usefully show).
merged <- merge(het, meta[, .(sample, site)], by = "sample", all.x = TRUE)
merged[is.na(site), site := "unknown"]

# Also merge missingness
smiss2 <- copy(smiss)
smiss2[, sample := IID]
merged <- merge(merged, smiss2[, .(sample, F_MISS)], by = "sample", all.x = TRUE)

# Compute outlier bounds
het_mean <- mean(merged$obs_het, na.rm = TRUE)
het_sd   <- sd(merged$obs_het, na.rm = TRUE)
sd_thresh <- if (!is.null(HET_OUTLIER_SD)) HET_OUTLIER_SD else 3
het_lo <- het_mean - sd_thresh * het_sd
het_hi <- het_mean + sd_thresh * het_sd
merged[, outlier := obs_het < het_lo | obs_het > het_hi]

# Build plot \u2014 no site coloring (see comment above); outliers marked red.
p_het <- ggplot(merged, aes(x = F_MISS, y = obs_het)) +
    geom_point(aes(color = outlier, shape = outlier), size = 2, alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = CLINEGO_NEUTRAL, "TRUE" = CLINEGO_REMOVED), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4), guide = "none") +
    geom_hline(yintercept = c(het_lo, het_hi), linetype = "dashed", color = CLINEGO_THRESHOLD, linewidth = 0.5) +
    geom_vline(xintercept = SAMPLE_MISS_THRESH, linetype = "dashed", color = CLINEGO_THRESHOLD, linewidth = 0.5) +
    labs(
        title = "Heterozygosity vs missingness per sample",
        x = "Fraction missing (F_MISS)",
        y = "Observed heterozygosity"
    ) +
    theme_clinego()

# Label outliers
outliers <- merged[outlier == TRUE]
if (nrow(outliers) > 0) {
    p_het <- p_het +
        ggrepel::geom_text_repel(
            data = outliers,
            aes(label = sample),
            size = 2.5, color = CLINEGO_REMOVED, max.overlaps = 20
        )
}
save_plot(p_het, "het_vs_missingness.png", width = 8, height = 6)

# Write sample heterozygosity table
het_table <- merged[, .(sample, site, F_MISS, obs_het, F_inbreeding = F, outlier_flag = outlier)]
fwrite(het_table, file.path(TABLES_DIR, "sample_heterozygosity.tsv"), sep = "\t", quote = FALSE)
message("INFO: Saved sample_heterozygosity.tsv")

#=============================================================================
# 6. SNP DENSITY HEATMAP BARS
#=============================================================================
dens_raw  <- fread(SNP_DENS_RAW)
dens_filt <- fread(SNP_DENS_FILT)

# vcftools --SNPdensity columns: CHROM, BIN_START, SNP_COUNT, VARIANTS/KB
setnames(dens_raw,  c("CHROM", "BIN_START", "SNP_COUNT", "VARIANTS_PER_KB"))
setnames(dens_filt, c("CHROM", "BIN_START", "SNP_COUNT", "VARIANTS_PER_KB"))

dens_raw[,  stage := "Raw"]
dens_filt[, stage := "Filtered"]

# Normalize chr names: raw VCF may have "chr1" prefix while filtered VCF has "1" (sed strips it)
dens_raw[,  CHROM := sub("^chr", "", CHROM)]
dens_filt[, CHROM := sub("^chr", "", CHROM)]

dens_both <- rbind(dens_raw, dens_filt)
# Raw on top, Filtered below: factor levels reversed so Raw is higher on y-axis
dens_both[, stage := factor(stage, levels = c("Filtered", "Raw"))]

# Compute bin midpoint in Mb
dens_both[, pos_mb := (BIN_START + 500000) / 1e6]

# Natural chromosome order
chr_order <- unique(dens_both$CHROM)
chr_order <- chr_order[order(suppressWarnings(as.integer(sub("^(chr)?", "", chr_order))), chr_order)]
dens_both[, CHROM := factor(CHROM, levels = chr_order)]

# facet_grid by chromosome: Raw (top) / Filtered (bottom) rows within each chr facet
p_density <- ggplot(dens_both, aes(x = pos_mb, y = stage, fill = SNP_COUNT)) +
    geom_tile(height = 0.9) +
    scale_fill_gradientn(
        colors = c("#F0F2F5", CLINEGO_RETAINED, "#0B4A42"),
        name = "SNPs/1Mb",
        na.value = "#F0F2F5"
    ) +
    facet_grid(CHROM ~ ., switch = "y") +
    labs(
        title = "SNP density (1Mb windows)",
        x = "Position (Mb)",
        y = NULL
    ) +
    theme_clinego() +
    theme(
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8, hjust = 1)
    )
save_plot(p_density, "snp_density_by_chr.png", width = 10, height = max(4, length(chr_order) * 1.2 + 1))

# Write density tables
fwrite(dens_raw[, .(CHROM, BIN_START, SNP_COUNT, VARIANTS_PER_KB)],
       file.path(TABLES_DIR, "snp_density_raw.tsv"), sep = "\t", quote = FALSE)
fwrite(dens_filt[, .(CHROM, BIN_START, SNP_COUNT, VARIANTS_PER_KB)],
       file.path(TABLES_DIR, "snp_density_filtered.tsv"), sep = "\t", quote = FALSE)
message("INFO: Saved snp_density_raw.tsv and snp_density_filtered.tsv")

#=============================================================================
# 7. FILTERING ATTRITION TABLE + WATERFALL PLOT
#=============================================================================
n_filt_snps <- count_vcf_snps(VCF_FILT)
n_ld_snps   <- count_vcf_snps(VCF_LD)

n_raw_samples_in_meta <- nrow(smiss)
n_after_sample_filter <- nrow(smiss[F_MISS <= SAMPLE_MISS_THRESH])
n_after_het_filter    <- if (!is.null(HET_OUTLIER_SD)) sum(!merged$outlier, na.rm = TRUE) else n_after_sample_filter
n_final_samples       <- nrow(meta)  # post all filtering (from aligned metadata)

# Relatedness removal count (read once here, reused by the Section 10 MDS plot below)
relatedness_removed_applied <- RELATEDNESS_ACTION == "remove" &&
    !is.null(RELATEDNESS_REMOVED_FILE) && file.exists(RELATEDNESS_REMOVED_FILE)
removed_ids <- character(0)
if (!is.null(RELATEDNESS_REMOVED_FILE) && file.exists(RELATEDNESS_REMOVED_FILE)) {
    removed_tbl <- fread(RELATEDNESS_REMOVED_FILE)
    if (nrow(removed_tbl) > 0) removed_ids <- removed_tbl$IID
}

# SNP attrition: separate MAF and SNP missingness stages
# n_below_maf is already computed above from maf_raw
n_after_maf <- n_raw_snps - n_below_maf

attrition <- data.table(
    stage = c(
        "Raw VCF",
        "After sample missingness filter",
        if (!is.null(HET_OUTLIER_SD)) "After het outlier removal" else NULL,
        if (relatedness_removed_applied) "After relatedness filter" else NULL,
        "After MAF filter",
        "After MAF + SNP missingness filter",
        "After LD pruning"
    ),
    n_samples = c(
        n_raw_samples_in_meta,
        n_after_sample_filter,
        if (!is.null(HET_OUTLIER_SD)) n_after_het_filter else NULL,
        if (relatedness_removed_applied) n_final_samples else NULL,
        # filter_vcf applies --keep and --maf/--geno in a single plink call, so
        # sample selection is already final by this point regardless of which
        # optional stages above ran — always n_final_samples, never a stale
        # pre-relatedness/pre-het count (that produced a nonsensical up-jump).
        n_final_samples,
        n_final_samples,
        n_final_samples
    ),
    n_snps = c(
        n_raw_snps,
        n_raw_snps,
        if (!is.null(HET_OUTLIER_SD)) n_raw_snps else NULL,
        if (relatedness_removed_applied) n_raw_snps else NULL,
        n_after_maf,
        n_filt_snps,
        n_ld_snps
    )
)
attrition[, pct_samples := round(100 * n_samples / n_raw_samples_in_meta, 1)]
attrition[, pct_snps    := round(100 * n_snps / n_raw_snps, 1)]
attrition[, stage := factor(stage, levels = stage)]

fwrite(attrition, file.path(TABLES_DIR, "filtering_summary.tsv"), sep = "\t", quote = FALSE)
message("INFO: Saved filtering_summary.tsv")

# SNP-only waterfall plot — exclude sample-only stages (sample miss + het outlier)
snp_stages <- c("Raw VCF", "After MAF filter", "After MAF + SNP missingness filter", "After LD pruning")
attrition_snps <- attrition[stage %in% snp_stages]
attrition_snps[, stage := factor(stage, levels = snp_stages)]

p_attrition <- ggplot(attrition_snps, aes(x = stage, y = n_snps)) +
    geom_col(fill = CLINEGO_RETAINED, width = 0.6) +
    geom_text(aes(label = paste0(n_snps, "\n(", pct_snps, "%)")),
              vjust = -0.3, size = 3) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "SNP retention through filtering pipeline",
        x = NULL,
        y = "Number of SNPs"
    ) +
    theme_clinego() +
    theme(axis.text.x = element_text(size = 9))
save_plot(p_attrition, "filtering_attrition.png", width = 8, height = 5)

#=============================================================================
# 8. DEPTH DISTRIBUTIONS (conditional)
#=============================================================================
if (HAS_DP && DEPTH_SAMPLE_FILE != 'NULL' && file.exists(DEPTH_SAMPLE_FILE)) {
    depth_s <- fread(DEPTH_SAMPLE_FILE)  # columns: INDV, N_SITES, MEAN_DEPTH
    depth_l <- fread(DEPTH_SITE_FILE)    # columns: CHROM, POS, MEAN_DEPTH

    # GBS VCFs often declare FORMAT/DP in the header but leave per-genotype DP empty ('.'),
    # which vcftools reports as -nan/0. Skip plot+table entirely rather than write a blank
    # PNG — Shiny falls back to an explanatory placeholder when these files don't exist.
    depth_ok <- any(is.finite(depth_s$MEAN_DEPTH) & depth_s$MEAN_DEPTH > 0) &&
                any(is.finite(depth_l$MEAN_DEPTH) & depth_l$MEAN_DEPTH > 0)

if (depth_ok) {
    # Bin per-site depth for histogram (avoid plotting millions of points)
    depth_l_hist <- depth_l[, .(count = .N), by = .(depth_bin = round(MEAN_DEPTH))]

    p_depth_s <- ggplot(depth_s, aes(x = MEAN_DEPTH)) +
        geom_histogram(bins = 40, fill = CLINEGO_RETAINED, color = "white", linewidth = 0.2) +
        geom_vline(xintercept = mean(depth_s$MEAN_DEPTH, na.rm = TRUE),
                   linetype = "dashed", color = CLINEGO_THRESHOLD) +
        labs(
            title = "Per-sample mean depth distribution",
            x = "Mean depth per sample",
            y = "Number of samples"
        ) +
        theme_clinego()

    p_depth_l <- ggplot(depth_l_hist, aes(x = depth_bin, y = count)) +
        geom_col(fill = CLINEGO_RETAINED, width = 0.8) +
        labs(
            title = "Per-SNP mean depth distribution",
            x = "Mean depth per SNP",
            y = "Number of SNPs"
        ) +
        theme_clinego()

    # Combine into one figure (2 panels)
    p_depth <- cowplot::plot_grid(p_depth_s, p_depth_l, ncol = 2, labels = c("A", "B"))
    cowplot::save_plot(file.path(PLOTS_DIR, "depth_distribution.png"), p_depth,
                       base_width = 10, base_height = 4)

    # Write depth summary table
    depth_summary <- data.table(
        metric = c("mean_depth_per_sample", "median_depth_per_sample",
                   "mean_depth_per_site",   "median_depth_per_site"),
        value  = c(
            round(mean(depth_s$MEAN_DEPTH, na.rm = TRUE), 2),
            round(median(depth_s$MEAN_DEPTH, na.rm = TRUE), 2),
            round(mean(depth_l$MEAN_DEPTH, na.rm = TRUE), 2),
            round(median(depth_l$MEAN_DEPTH, na.rm = TRUE), 2)
        )
    )
    fwrite(depth_summary, file.path(TABLES_DIR, "depth_summary.tsv"), sep = "\t", quote = FALSE)
    message("INFO: Saved depth_distribution.png and depth_summary.tsv")
} else {
    # Snakemake declares depth_distribution.png a required output whenever HAS_DP is TRUE
    # (a header-only check) — it has no way to know at parse time that genotype DP is
    # actually unpopulated. Leave a 0-byte marker so the rule's output contract is satisfied;
    # Shiny's file_ok() treats file.size==0 as unavailable and still falls back to the
    # explanatory placeholder (see mod_processing.R depth_status).
    file.create(file.path(PLOTS_DIR, "depth_distribution.png"))
    message("WARNING: FORMAT/DP is declared in the VCF header but no finite, positive depth values were found (typical for GBS with unpopulated DP) — skipping depth_distribution.png and depth_summary.tsv.")
}
}

#=============================================================================
# 9. RELATEDNESS (IBS ALLELE-SHARING) DISTRIBUTION — always drawn
#=============================================================================
if (file.exists(RELATEDNESS_PAIRS)) {
    rel_pairs <- fread(RELATEDNESS_PAIRS)
} else {
    rel_pairs <- data.table(IID1 = character(0), IID2 = character(0), IBS = numeric(0))
}

if (nrow(rel_pairs) > 0) {
    if (!is.null(PI_HAT_THRESH)) {
        # Pair count above threshold is reported in the Shiny hover note, not on-plot
        # (rule: plots carry only data-derived annotations, not instructional text).
        p_rel <- ggplot(rel_pairs, aes(x = IBS)) +
            geom_histogram(
                aes(fill = IBS > PI_HAT_THRESH),
                binwidth = 0.02, color = "white", linewidth = 0.1
            ) +
            scale_fill_manual(
                values = c("FALSE" = CLINEGO_RETAINED, "TRUE" = CLINEGO_REMOVED),
                labels = c("FALSE" = "Below threshold", "TRUE" = "Above threshold"),
                name = NULL
            ) +
            geom_vline(xintercept = PI_HAT_THRESH, linetype = "dashed", color = CLINEGO_THRESHOLD)
    } else {
        p_rel <- ggplot(rel_pairs, aes(x = IBS)) +
            geom_histogram(binwidth = 0.02, fill = CLINEGO_RETAINED, color = "white", linewidth = 0.1)
    }

    p_rel <- p_rel +
        labs(
            title = "Pairwise relatedness distribution (IBS allele-sharing)",
            x = "IBS (proportion shared alleles)",
            y = "Number of sample pairs"
        ) +
        theme_clinego()
    save_plot(p_rel, "relatedness_distribution.png")
} else {
    message("WARNING: Relatedness pairs file is empty or missing — skipping relatedness_distribution.png (need >=2 samples for pairwise IBS).")
}

#=============================================================================
# 10. RELATEDNESS MDS (IBS allele-sharing as distance) — always drawn
#=============================================================================
# Reuses rel_pairs from Section 9 and removed_ids from Section 7 (compute
# once). Discarded samples (preview from relatedness_filter, or actual
# removals when Filter.relatedness_action is 'remove') are marked red so
# users see WHICH samples cluster as duplicates/clones, not just that a
# spike exists in the histogram.
all_ids <- unique(c(rel_pairs$IID1, rel_pairs$IID2))

mds_coords <- if (nrow(rel_pairs) > 0 && length(all_ids) >= 3) {
    tryCatch({
        n <- length(all_ids)
        dist_mat <- matrix(0, n, n, dimnames = list(all_ids, all_ids))
        idx1 <- match(rel_pairs$IID1, all_ids)
        idx2 <- match(rel_pairs$IID2, all_ids)
        dist_mat[cbind(idx1, idx2)] <- 1 - rel_pairs$IBS
        dist_mat[cbind(idx2, idx1)] <- 1 - rel_pairs$IBS

        fit <- cmdscale(as.dist(dist_mat), k = 2)
        data.table(sample = rownames(fit), MDS1 = fit[, 1], MDS2 = fit[, 2])
    }, error = function(e) {
        message("WARNING: relatedness MDS failed (", conditionMessage(e), ") — skipping relatedness_mds.png content.")
        NULL
    })
} else {
    message("WARNING: Fewer than 3 samples with pairwise IBS — skipping relatedness_mds.png content.")
    NULL
}

if (!is.null(mds_coords)) {
    mds_coords[, discarded := sample %in% removed_ids]

    p_mds <- ggplot(mds_coords, aes(x = MDS1, y = MDS2, color = discarded)) +
        geom_point(size = 2.5, alpha = 0.85) +
        scale_color_manual(
            values = c("FALSE" = CLINEGO_NEUTRAL, "TRUE" = CLINEGO_REMOVED),
            labels = c("FALSE" = "Retained", "TRUE" = "Discarded (relatedness)"),
            name = NULL
        ) +
        labs(
            title = "Sample relatedness MDS (IBS allele-sharing)",
            x = "MDS1", y = "MDS2"
        ) +
        theme_clinego()

    discarded_pts <- mds_coords[discarded == TRUE]
    if (nrow(discarded_pts) > 0) {
        discarded_pts[, label := sub("_DUP$", "", sample)]
        p_mds <- p_mds +
            ggrepel::geom_text_repel(
                data = discarded_pts, aes(label = label),
                size = 2.5, color = CLINEGO_REMOVED, max.overlaps = Inf,
                show.legend = FALSE
            )
    }
} else {
    # Empty placeholder so the Snakemake output contract is always satisfied.
    p_mds <- ggplot() +
        annotate("text", x = 0, y = 0, label = "Not enough samples for relatedness MDS", color = CLINEGO_NEUTRAL) +
        theme_void()
}
save_plot(p_mds, "relatedness_mds.png")

message("INFO: All QC plots and tables complete.")
