#!/usr/bin/env Rscript
# Find significant SNPs from association analysis

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(qvalue)
    library(stringr)
    library(parallel)
})

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/io_selected_snps.R")
source("/pipeline/scripts/R/lib/regions.R")          # distance_for_chr()
source("/pipeline/scripts/R/lib/region_distance.R")  # resolve_clumping_distance_strict()
source("/pipeline/scripts/R/lib/sig_snps.R")

args = commandArgs(trailingOnly = TRUE)
################
ASSOC_TABLE       = args[1]  # p-values table
ADJUST            = args[2]  # e.g., "bonf_0.05"
CLUMPING_SPEC     = args[3]  # bp, or auto_genome_wide | auto_per_chromosome
METHOD = args[4]  # EMMAX or LFMM
CPU    = as.numeric(args[5])
OUTPUT = args[6]
REGIME         = if (length(args) >= 7) args[7] else "--snp"   # --snp | --wza
LD_DECAY_PATH  = if (length(args) >= 8) args[8] else "NULL"    # resolves an auto_* spec
R2_THRESHOLD   = if (length(args) >= 9) as.numeric(args[9]) else 0.2
LD_DECAY_GROUP = if (length(args) >= 10) args[10] else "All"
################
if (!REGIME %in% c("--snp", "--wza"))
    stop("find_sig_snps.R: argument 7 must be --snp or --wza, got '", REGIME, "'")
IS_WZA = REGIME == "--wza"

# The cross-trait overlap annotation reads the shared clumping distance; an auto_*
# spec is resolved from the LD-decay table exactly as create_regions.R resolves it.
CLUMPING_DISTANCE <- resolve_clumping_distance_strict(CLUMPING_SPEC, LD_DECAY_PATH,
                                                      R2_THRESHOLD, LD_DECAY_GROUP)

message(paste0('INFO: Finding significant SNPs for ', METHOD))
message(paste0('INFO: Adjustment: ', ADJUST))
message(paste0('INFO: SNP clumping distance: ', format_clumping_distance(CLUMPING_DISTANCE)))

adjustment <- str_split(ADJUST, '_')[[1]][1]
pval_value <- as.numeric(str_split(ADJUST, '_')[[1]][2])

wza_meta_cols <- if (IS_WZA) c("n_snps", "mean_maf") else character(0)

pvals <- fread(ASSOC_TABLE, colClasses = c(chr = "character"))
message(str(pvals))

result <- find_significant_snps_per_trait(pvals, adjustment, pval_value, CPU,
                                           exclude_cols = wza_meta_cols,
                                           is_wza = IS_WZA)

# An empty WZA sig-windows table has two very different causes: every window was
# tested and none passed (no signal), or the WZA statistic's own floor sits above the
# threshold so no window COULD pass (finding 6ab65e). compute_wza.R logs why the floor
# is where it is; this says, per trait, which of the two happened here. Per-SNP tables
# are untouched — the same comparison there is the ordinary "nothing was significant".
if (IS_WZA) {
    for (tr in setdiff(colnames(pvals), c("SNPID", "chr", "pos", wza_meta_cols))) {
        thr  <- compute_pval_threshold(pvals[[tr]], adjustment, pval_value)
        pv   <- pvals[[tr]][!is.na(pvals[[tr]])]
        if (length(pv) == 0L || thr$status != "ok") next
        n_called <- sum(pv <= thr$threshold)
        if (n_called == 0L) {
            message(paste0("WARNING: WZA trait=", tr,
                           " — 0 of ", length(pv), " windows called. Best window p=",
                           signif(min(pv), 4), " vs threshold ", signif(thr$threshold, 4),
                           " (ratio ", signif(min(pv) / thr$threshold, 3),
                           "x). Check the WZA floor reported by compute_wza.R before ",
                           "reading this as absence of signal."))
        }
    }
}

sig <- annotate_cross_trait_overlaps(result$sig_snps, CLUMPING_DISTANCE)
sig[, method := METHOD]

final <- sig[, .(SNPID, chr, pos, pvalue, pval_threshold, method,
                  trait, overlap_traits, overlap_snps, overlap_distance)]
write_selected_snps(final, OUTPUT)

message(paste0('INFO: Saved ', nrow(final), ' significant SNPs to ', OUTPUT))
