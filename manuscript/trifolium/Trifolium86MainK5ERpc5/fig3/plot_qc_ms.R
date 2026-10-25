#!/usr/bin/env Rscript
# DERIVED COPY of five sections of scripts/plot_qc_processing.R (pipeline, HEAD 971c38a) for the Trifolium supplement
# figure s_qc (plan ~/.claude/plans/let-s-follow-the-plan-lively-hellman.md, Phase 3, item B1). Lines marked
# `# DERIVED` differ; each is an entry in manuscript_fig3/pipeline_deviations.md (status pending):
#   D23 theme_ms() + minou colours (retained MINOU teal, removed MINOU red, threshold lines MS_INK; density tracks
#       MS_SEQUENTIAL) instead of theme_clinego() and the CLINEGO_* constants; ms_save() at 180 mm instead of ggsave()
#   D24 no plot titles
#   D25 panel c (relatedness MDS) without text labels: the pipeline labels discarded plants by sampling id, and
#       manuscript figures print no raw ids; discarded plants drawn last so no retained point covers them
#   D26 panel d (SNP density): raw and filtered tracks each with its own fill scale and legend (ggnewscale), so the
#       filtered track is not lost in the raw scale; the pipeline draws both with one scale
#   D29 axis titles "Number of plants" / "Number of plant pairs" instead of "Number of samples" / "Number of sample
#       pairs" (manuscript vocabulary: the units are plants)
#   D27 one composite (patchwork, tags a-e, bold 8 pt as Figure 5) instead of one PNG per section; count axes with
#       thousands separators and negative ticks with a true minus sign (manuscript labels rule)
# Sections kept (pipeline line numbers at 971c38a): a = 2 sample missingness (:114-134), b = 9 relatedness IBS
# (:464-502), c = 10 relatedness MDS (:505-566), d = 6 SNP density (:247-300), e = 4 MAF (:160-184). The depth panel
# (section 8) is left out (spec B1); the other sections are not drawn.
#
# Args: the 25 positional arguments workflow/rules/processing.smk (rule plot_qc_processing) passes to the pipeline
# script, unchanged, then
#   26 out_stem  - output stem, written as <out_stem>.{svg,png}                                   # DERIVED (D27)
#   27 counts    - TSV of the asserted counts (sources/ of the supplement)                         # DERIVED (D27)
# Args 1, 3, 5-7, 10-12, 14, 16-21 and 25 are read by sections not drawn here and are accepted only so the call keeps
# the rule's argument list.
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(patchwork)   # DERIVED (D27)
    library(ggnewscale)  # DERIVED (D26)
})

source("/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R")  # DERIVED (D23)
source("/pipeline/scripts/R/lib/gff_parsing.R")   # normalize_chr()

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 27)                                          # DERIVED (D27)

SAMPLE_MISS_FILE  <- args[2]
MAF_RAW_FILE      <- args[4]
SNP_DENS_RAW      <- args[8]
SNP_DENS_FILT     <- args[9]
SAMPLE_MISS_THRESH <- as.numeric(args[13])
MAF_THRESH         <- as.numeric(args[15])
RELATEDNESS_PAIRS  <- args[22]
PI_HAT_THRESH      <- if (args[23] == 'NULL') NULL else as.numeric(args[23])
RELATEDNESS_REMOVED_FILE <- if (args[24] == 'NULL') NULL else args[24]
OUT_STEM   <- args[26]                                                 # DERIVED (D27)
COUNTS_TSV <- args[27]                                                 # DERIVED (D27)
stopifnot(!is.null(PI_HAT_THRESH), !is.null(RELATEDNESS_REMOVED_FILE))  # DERIVED (D27): this run sets both

RETAINED  <- MINOU[["teal"]]                                           # DERIVED (D23)
REMOVED   <- MINOU[["red"]]                                            # DERIVED (D23)
THRESHOLD <- MS_INK                                                    # DERIVED (D23)
lab_count <- scales::label_comma()                                     # DERIVED (D27)
lab_signed <- scales::label_number(style_negative = "minus")           # DERIVED (D27)

#=============================================================================
# a. SAMPLE MISSINGNESS DISTRIBUTION (section 2)
#=============================================================================
smiss <- fread(SAMPLE_MISS_FILE, colClasses = list(character = "IID"))
p_smiss <- ggplot(smiss, aes(x = F_MISS)) +
    geom_histogram(
        aes(fill = F_MISS > SAMPLE_MISS_THRESH),
        binwidth = 0.02, color = "white", linewidth = 0.2
    ) +
    scale_fill_manual(
        values = c("FALSE" = RETAINED, "TRUE" = REMOVED),              # DERIVED (D23)
        labels = c("FALSE" = "Retained", "TRUE" = "Removed"),
        name = NULL
    ) +
    geom_vline(xintercept = SAMPLE_MISS_THRESH, linetype = "dashed", color = THRESHOLD, linewidth = 0.3) +  # DERIVED (D23)
    labs(
        x = "Fraction missing genotypes (F_MISS)",                     # DERIVED (D24)
        y = "Number of plants"                                         # DERIVED (D29)
    ) +
    theme_ms()                                                         # DERIVED (D23)

#=============================================================================
# b. RELATEDNESS (IBS ALLELE-SHARING) DISTRIBUTION (section 9)
#=============================================================================
rel_pairs <- fread(RELATEDNESS_PAIRS)
stopifnot(nrow(rel_pairs) > 0)                                         # DERIVED (D27): the empty branch never runs here
p_rel <- ggplot(rel_pairs, aes(x = IBS)) +
    geom_histogram(
        aes(fill = IBS > PI_HAT_THRESH),
        binwidth = 0.02, color = "white", linewidth = 0.1
    ) +
    scale_fill_manual(
        values = c("FALSE" = RETAINED, "TRUE" = REMOVED),              # DERIVED (D23)
        labels = c("FALSE" = "Below threshold", "TRUE" = "Above threshold"),
        name = NULL
    ) +
    geom_vline(xintercept = PI_HAT_THRESH, linetype = "dashed", color = THRESHOLD, linewidth = 0.3) +  # DERIVED (D23)
    scale_y_continuous(labels = lab_count) +                           # DERIVED (D27)
    labs(
        x = "IBS (proportion shared alleles)",                         # DERIVED (D24)
        y = "Number of plant pairs"                                    # DERIVED (D29)
    ) +
    theme_ms()                                                         # DERIVED (D23)

#=============================================================================
# c. RELATEDNESS MDS (section 10; removed_ids as section 7 reads them)
#=============================================================================
removed_ids <- character(0)
if (!is.null(RELATEDNESS_REMOVED_FILE) && file.exists(RELATEDNESS_REMOVED_FILE)) {
    removed_tbl <- fread(RELATEDNESS_REMOVED_FILE)
    if (nrow(removed_tbl) > 0) removed_ids <- removed_tbl$IID
}
all_ids <- unique(c(rel_pairs$IID1, rel_pairs$IID2))
n <- length(all_ids)
dist_mat <- matrix(0, n, n, dimnames = list(all_ids, all_ids))
idx1 <- match(rel_pairs$IID1, all_ids)
idx2 <- match(rel_pairs$IID2, all_ids)
dist_mat[cbind(idx1, idx2)] <- 1 - rel_pairs$IBS
dist_mat[cbind(idx2, idx1)] <- 1 - rel_pairs$IBS
fit <- cmdscale(as.dist(dist_mat), k = 2)
mds_coords <- data.table(sample = rownames(fit), MDS1 = fit[, 1], MDS2 = fit[, 2])
mds_coords[, discarded := sample %in% removed_ids]
setorder(mds_coords, discarded)                                        # DERIVED (D25): discarded drawn last, on top
p_mds <- ggplot(mds_coords, aes(x = MDS1, y = MDS2, color = discarded)) +
    geom_point(size = 1.2, alpha = 0.85) +                             # DERIVED (D23): 2.5 at 8 x 5 in
    scale_color_manual(
        values = c("FALSE" = RETAINED, "TRUE" = REMOVED),              # DERIVED (D23)
        labels = c("FALSE" = "Retained", "TRUE" = "Discarded (relatedness)"),
        name = NULL
    ) +
    scale_x_continuous(labels = lab_signed) +                          # DERIVED (D27)
    scale_y_continuous(labels = lab_signed) +                          # DERIVED (D27)
    labs(x = "MDS1", y = "MDS2") +                                     # DERIVED (D24)
    theme_ms()                                                         # DERIVED (D23)
# discarded-sample text labels not drawn                               # DERIVED (D25)

#=============================================================================
# d. SNP DENSITY HEATMAP BARS (section 6)
#=============================================================================
dens_raw  <- fread(SNP_DENS_RAW)
dens_filt <- fread(SNP_DENS_FILT)
setnames(dens_raw,  c("CHROM", "BIN_START", "SNP_COUNT", "VARIANTS_PER_KB"))
setnames(dens_filt, c("CHROM", "BIN_START", "SNP_COUNT", "VARIANTS_PER_KB"))
dens_raw[,  stage := "Raw"]
dens_filt[, stage := "Filtered"]
dens_raw[,  CHROM := normalize_chr(CHROM)]
dens_filt[, CHROM := normalize_chr(CHROM)]
dens_both <- rbind(dens_raw, dens_filt)
dens_both[, stage := factor(stage, levels = c("Filtered", "Raw"))]
dens_both[, pos_mb := (BIN_START + 500000) / 1e6]
chr_order <- unique(dens_both$CHROM)
chr_order <- chr_order[order(suppressWarnings(as.integer(sub("^(chr)?", "", chr_order))), chr_order)]
dens_both[, CHROM := factor(CHROM, levels = chr_order)]
dens_scale <- function(name, order) scale_fill_gradientn(              # DERIVED (D26)
    colors = MS_SEQUENTIAL, name = name, labels = lab_count, na.value = "#FFFFFF",  # DERIVED (D23)
    guide = guide_colourbar(order = order, theme = theme(legend.key.width = grid::unit(32, "mm"), legend.key.height = grid::unit(2, "mm"),
                                          legend.title.position = "top")))
p_density <- ggplot(mapping = aes(x = pos_mb, y = stage)) +          # DERIVED (D26): fill mapped per layer
    geom_tile(data = dens_both[stage == "Raw"], aes(fill = SNP_COUNT), height = 0.9) +  # DERIVED (D26)
    dens_scale("Raw, SNPs/1Mb", 1) +                                   # DERIVED (D26)
    new_scale_fill() +                                                 # DERIVED (D26)
    geom_tile(data = dens_both[stage == "Filtered"], aes(fill = SNP_COUNT), height = 0.9) +  # DERIVED (D26)
    dens_scale("Filtered, SNPs/1Mb", 2) +                              # DERIVED (D26)
    scale_y_discrete(limits = levels(dens_both$stage)) +               # DERIVED (D26): Raw on top, as the pipeline
    facet_grid(CHROM ~ ., switch = "y") +
    labs(
        x = "Position (Mb)",                                           # DERIVED (D24)
        y = NULL
    ) +
    theme_ms() +                                                       # DERIVED (D23)
    theme(
        legend.position = "bottom",                                    # DERIVED (D23): was right
        legend.spacing.x = grid::unit(8, "mm"),                        # DERIVED (D26): keeps the two keys' labels apart
        panel.grid = element_blank(),
        panel.spacing = grid::unit(0.6, "mm"),                         # DERIVED (D23)
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = MS_TICK, hjust = 1)  # DERIVED (D23)
    )

#=============================================================================
# e. MAF DISTRIBUTION (section 4, raw only)
#=============================================================================
maf_raw <- fread(MAF_RAW_FILE)
n_below_maf <- nrow(maf_raw[MAF < MAF_THRESH])
p_maf <- ggplot(maf_raw, aes(x = MAF)) +
    geom_histogram(
        aes(fill = MAF < MAF_THRESH),
        binwidth = 0.01, color = NA
    ) +
    scale_fill_manual(
        values = c("FALSE" = RETAINED, "TRUE" = REMOVED),              # DERIVED (D23)
        labels = c("FALSE" = "Retained", "TRUE" = "Removed"),
        name = NULL
    ) +
    geom_vline(xintercept = MAF_THRESH, linetype = "dashed", color = THRESHOLD, linewidth = 0.3) +  # DERIVED (D23)
    scale_y_continuous(labels = lab_count) +                           # DERIVED (D27)
    labs(
        x = "Minor allele frequency (MAF)",                            # DERIVED (D24)
        y = "Number of SNPs"
    ) +
    theme_ms()                                                         # DERIVED (D23)

#=============================================================================
# asserts + counts (spec B1)                                           # DERIVED (D27): this whole block
#=============================================================================
pass_miss <- smiss[F_MISS <= SAMPLE_MISS_THRESH, IID]
n_pairs_above <- rel_pairs[IBS > PI_HAT_THRESH, .N]
dsum <- dens_both[, .(total = sum(SNP_COUNT), min = min(SNP_COUNT), max = max(SNP_COUNT), windows = .N), by = stage]
counts <- data.table(
    panel = c("a", "a", "a", "b", "b", "c", "c", "c", "c", "e", "e", rep("d", 8)),
    quantity = c("plants", "plants above the sample missingness threshold (drawn red)", "sample missingness threshold",
                 "plant pairs", "plant pairs above the IBS threshold (drawn red)",
                 "plants in the MDS", "plants in the MDS that passed the missingness filter",
                 "plants drawn red (relatedness_removed.tsv)", "IBS threshold",
                 "raw SNPs", "raw SNPs below the MAF threshold (drawn red)",
                 paste(rep(c("Raw", "Filtered"), each = 4), c("SNPs (sum over windows)", "1 Mb windows",
                                                               "SNPs per window, min", "SNPs per window, max"))),
    value = c(nrow(smiss), smiss[F_MISS > SAMPLE_MISS_THRESH, .N], SAMPLE_MISS_THRESH,
              nrow(rel_pairs), n_pairs_above,
              nrow(mds_coords), sum(mds_coords$sample %in% pass_miss), mds_coords[discarded == TRUE, .N], PI_HAT_THRESH,
              nrow(maf_raw), n_below_maf,
              unlist(lapply(c("Raw", "Filtered"), function(g) dsum[stage == g, c(total, windows, min, max)]))))
print(counts)
chk <- c(a_plants = nrow(smiss) == 86L, a_red = smiss[F_MISS > SAMPLE_MISS_THRESH, .N] == 22L,
         c_red = mds_coords[discarded == TRUE, .N] == 7L,
         c_red_ids = setequal(mds_coords[discarded == TRUE, sample], removed_ids),
         c_red_ibs = all(removed_tbl$IBS > PI_HAT_THRESH),
         d_raw_total = dsum[stage == "Raw", total] == 7163793L,        # = raw SNPs (filtering_summary.tsv)
         d_filt_total = dsum[stage == "Filtered", total] == 209466L)   # = after MAF + SNP missingness
message(sprintf("MDS holds %d plants; %d of them passed the missingness filter (%d passed it); identical sets: %s",
                nrow(mds_coords), sum(mds_coords$sample %in% pass_miss), length(pass_miss),
                setequal(mds_coords$sample, pass_miss)))
if (!all(chk)) stop("s_qc assertion failed: ", paste(names(chk)[!chk], collapse = ", "))
message("s_qc asserts passed: ", paste(names(chk), collapse = ", "))
fwrite(counts, COUNTS_TSV, sep = "\t")

#=============================================================================
# composite                                                            # DERIVED (D27): this whole block
#=============================================================================
W_MM <- 180; H_MM <- 125
bottom <- (p_density | p_maf) + plot_layout(widths = c(2, 1))
fig <- ((p_smiss | p_rel | p_mds) / bottom) + plot_layout(heights = c(1, 1.45)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 8, colour = MS_INK))
ms_save(OUT_STEM, fig, W_MM, H_MM)
message(sprintf("s_qc: %d x %d mm -> %s.{svg,png}", W_MM, H_MM, OUT_STEM))
