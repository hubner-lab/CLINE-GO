#!/usr/bin/env Rscript
# LD-decay pairs per group: PLINK 1.9 --r2 (squared allele-count correlation, missing
# calls excluded pairwise) between SNPs of one chromosome within max_distance, summed
# per bp distance into PopLDdecay-layout .stat.gz files that ld_decay_analyze.R reads:
# stat_gw/<stem>.stat.gz and, when the scope includes per_chromosome,
# stat_chr/<stem>_<chr>.stat.gz. <stem> is the group's sample-list name in
# manifest.tsv. When all pairs within the window fit the pair budget the run is
# exhaustive; otherwise PLINK pairs a seeded sample of anchor SNPs with every SNP of
# their window (--ld-snp-list), which thins all distances evenly.
#
# Chromosomes are read with --chr-set 95 no-x no-y no-xy no-mt: in PLINK's default
# human mode a numeric chromosome 23-26 is X/Y/XY/MT (Y calls of non-males are
# dropped) and a code above 26 is an error.
library(data.table)
source("/pipeline/scripts/R/lib/ld_decay.R")

args <- commandArgs(trailingOnly = TRUE)
################
VCF_PATH         <- args[1]  # filtered VCF (observed calls)
SAMPLE_LISTS_DIR <- args[2]  # ld_decay_prepare.R output: manifest.tsv + <stem>.txt
MAX_DIST_KB      <- as.numeric(args[3])
MAX_PAIRS        <- as.numeric(args[4])
SEED             <- as.integer(args[5])
SCOPE            <- args[6]  # genome_wide | per_chromosome | both
STAT_GW_DIR      <- args[7]
STAT_CHR_DIR     <- args[8]
OUT_SUMMARY      <- args[9]  # one row per group: SNPs, anchors, pairs
################
stopifnot(length(args) == 9, MAX_DIST_KB > 0, MAX_PAIRS > 0,
          SCOPE %in% c("genome_wide", "per_chromosome", "both"))
PLINK_CHR <- c("--chr-set", "95", "no-x", "no-y", "no-xy", "no-mt", "--allow-extra-chr")
max_dist <- MAX_DIST_KB * 1000

manifest <- fread(file.path(SAMPLE_LISTS_DIR, "manifest.tsv"), colClasses = "character")
if (!"list_file" %in% names(manifest)) manifest[, list_file := group]
dir.create(STAT_GW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_CHR_DIR, recursive = TRUE, showWarnings = FALSE)
# ld_decay_analyze.R reads whatever file a group or chromosome has, so every one must
# come from this run: a chromosome without pairs now would keep an earlier run's.
unlink(list.files(c(STAT_GW_DIR, STAT_CHR_DIR), pattern = "\\.stat\\.gz$", full.names = TRUE))
# PLINK's pair table is gigabytes at the default budget: it stays beside the output.
tmp <- ld_scratch_dir(OUT_SUMMARY, "pairs_tmp")

plink <- function(...) {
    a <- c(...)
    out <- suppressWarnings(system2("plink", a, stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) {
        stop("plink failed (status ", status, "): plink ", paste(a, collapse = " "),
             "\n", paste(tail(out, 15), collapse = "\n"))
    }
    invisible(out)
}

all_set <- file.path(tmp, "all")
plink("--vcf", VCF_PATH, "--const-fid", PLINK_CHR, "--vcf-half-call", "missing",
      "--make-bed", "--threads", "1", "--out", all_set)
bim <- fread(paste0(all_set, ".bim"), header = FALSE,
             colClasses = c(V1 = "character", V2 = "character"))
setnames(bim, c("chr", "id", "cm", "pos", "a1", "a2"))
# PLINK rewrites chromosome names ("chr1" -> "1"); the per-chromosome files must carry
# the VCF's names, which ld_decay_prepare.R lists in chromosomes.txt.
vcf_chr <- system2("bcftools", c("query", "-f", shQuote("%CHROM\\n"), shQuote(VCF_PATH)),
                   stdout = TRUE)
if (length(vcf_chr) != nrow(bim)) {
    stop("PLINK kept ", nrow(bim), " of the VCF's ", length(vcf_chr), " variants")
}
chr_name <- unique(data.table(code = bim$chr, name = vcf_chr))
if (anyDuplicated(chr_name$code)) stop("PLINK merged distinct VCF chromosome names")
chr_name <- setNames(chr_name$name, chr_name$code)
# Internal variant IDs: VCF IDs may be "." or repeated, and --ld-snp-list needs them unique.
bim[, id := paste0("v", .I)]
fwrite(bim, paste0(all_set, ".bim"), sep = "\t", col.names = FALSE)

plan <- ld_anchor_plan(bim$chr, bim$pos, max_dist, MAX_PAIRS)
anchor_args <- character(0)
if (!plan$exhaustive) {
    set.seed(SEED)
    anchors <- bim$id[sort(sample.int(nrow(bim), plan$n_anchors))]
    anchor_file <- file.path(tmp, "anchors.txt")
    writeLines(anchors, anchor_file)
    anchor_args <- c("--ld-snp-list", anchor_file)
}
message("INFO: ", nrow(bim), " SNPs, ", round(plan$total_pairs), " pairs within ",
        MAX_DIST_KB, " kb; ",
        if (plan$exhaustive) "all pairs used" else
            paste0(plan$n_anchors, " anchors (seed ", SEED, "), about ",
                   round(plan$n_anchors * plan$mean_partners), " pairs"))

per_chr <- SCOPE %in% c("per_chromosome", "both")
rows <- list()
for (i in seq_len(nrow(manifest))) {
    grp  <- manifest$group[i]
    stem <- manifest$list_file[i]
    keep_args <- character(0)
    if (grp != "All") {
        smp  <- readLines(file.path(SAMPLE_LISTS_DIR, paste0(stem, ".txt")))
        keep <- file.path(tmp, paste0(stem, ".keep"))
        fwrite(data.table("0", smp), keep, sep = "\t", col.names = FALSE)
        keep_args <- c("--keep", keep)
    }
    out <- file.path(tmp, stem)
    plink("--bfile", all_set, PLINK_CHR, keep_args, "--r2", anchor_args,
          "--ld-window", "999999", "--ld-window-kb", MAX_DIST_KB, "--ld-window-r2", "0",
          "--threads", "1", "--out", out)
    pairs <- fread(paste0(out, ".ld"), select = c("CHR_A", "BP_A", "SNP_A", "BP_B", "SNP_B", "R2"),
                   colClasses = c(CHR_A = "character", SNP_A = "character", SNP_B = "character"))
    unlink(paste0(out, ".ld"))
    st <- ld_pairs_to_stat(pairs, max_dist, dedupe = !plan$exhaustive)
    ld_write_stat(st$gw, file.path(STAT_GW_DIR, paste0(stem, ".stat.gz")))
    if (per_chr) {
        for (ch in unique(st$chr$chr)) {
            ld_write_stat(st$chr[chr == ch],
                          file.path(STAT_CHR_DIR, paste0(stem, "_", chr_name[[ch]], ".stat.gz")))
        }
    }
    rows[[grp]] <- data.table(group = grp, n_snps = nrow(bim), exhaustive = plan$exhaustive,
                              n_anchors = plan$n_anchors, n_pairs = sum(st$gw$n_pairs),
                              max_distance_kb = MAX_DIST_KB, seed = SEED)
    message("INFO: ", grp, ": ", sum(st$gw$n_pairs), " pairs")
}

unlink(tmp, recursive = TRUE)
fwrite(rbindlist(rows), OUT_SUMMARY, sep = "\t")
message("INFO: written ", OUT_SUMMARY)
