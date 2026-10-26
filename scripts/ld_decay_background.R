#!/usr/bin/env Rscript
# Background LD of unlinked markers, per LD-decay group (Breseghello & Sorrells 2006,
# as applied by Laido et al. 2014): r^2 between markers on DIFFERENT chromosomes, for
# a seeded sample of markers per chromosome (Vos et al. 2017 sampled 50), its mean and
# its critical value, the `quantile` of sqrt(r^2), squared. PLINK 1.9 --r2 inter-chr,
# with the chromosome set and missing-call handling of ld_decay_pairs.R, so the
# background and the decay curve come from one estimator. The contigs in
# LDdecay.exclude_chromosomes give no markers, as they give the curve no pairs.
# A dataset with a single chromosome has no unlinked pairs: its rows carry NA.
library(data.table)
source("/pipeline/scripts/R/lib/ld_decay.R")

args <- commandArgs(trailingOnly = TRUE)
################
VCF_PATH         <- args[1]  # filtered VCF, the one ld_decay_pairs.R reads
SAMPLE_LISTS_DIR <- args[2]  # ld_decay_prepare.R output: manifest.tsv + <group>.txt
MARKERS_PER_CHR  <- as.integer(args[3])
QUANTILE         <- as.numeric(args[4])
SEED             <- as.integer(args[5])
OUT_TABLE        <- args[6]
EXCLUDE          <- args[7]  # LDdecay.exclude_chromosomes, comma-separated, or "NULL"
################
stopifnot(length(args) == 7, MARKERS_PER_CHR >= 2, QUANTILE > 0, QUANTILE < 1)

manifest <- fread(file.path(SAMPLE_LISTS_DIR, "manifest.tsv"), colClasses = "character")
if (!"list_file" %in% names(manifest)) manifest[, list_file := group]
PLINK_CHR <- c("--chr-set", "95", "--allow-extra-chr")   # as in ld_decay_pairs.R
tmp <- ld_scratch_dir(SAMPLE_LISTS_DIR, "background_tmp")
pos <- fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\n'", shQuote(VCF_PATH)),
             header = FALSE, col.names = c("chr", "pos"),
             colClasses = c("character", "integer"), tmpdir = tmp)
excluded <- ld_excluded_contigs(pos$chr, EXCLUDE)
if (length(excluded)) {
    message("INFO: Left out (LDdecay.exclude_chromosomes): ", paste(excluded, collapse = ", "))
    pos <- pos[!chr %in% excluded]
}
n_chr <- uniqueN(pos$chr)
message("INFO: ", nrow(pos), " SNPs on ", n_chr, " chromosomes")

# At most MAX_MARKERS markers in all (ld_background_markers()), however fragmented the
# assembly: unlinked pairs grow with the square of the marker count.
MAX_MARKERS <- 2000L
bm <- ld_background_markers(pos$chr, pos$pos, MARKERS_PER_CHR, MAX_MARKERS, SEED)
picked <- bm$markers
per_chr <- bm$per_chr
if (per_chr < MARKERS_PER_CHR) {
    message("INFO: ", n_chr, " chromosomes: ", per_chr, " marker(s) per chromosome, ",
            nrow(picked), " in all (cap ", MAX_MARKERS, ")")
}
ranges <- file.path(tmp, "markers.range")
fwrite(picked[, .(chr, pos, pos, paste0(chr, ":", pos))], ranges,
       sep = "\t", col.names = FALSE)

# Read the VCF once: the sampled markers into a small PLINK binary set, which every
# group then reads (a VCF pass per group cost 101 passes per simulated replicate).
markers <- file.path(tmp, "markers")
if (n_chr >= 2) {
    status <- system2("plink", c("--vcf", VCF_PATH, "--const-fid", PLINK_CHR,
                                 "--vcf-half-call", "missing", "--extract", "range", ranges,
                                 "--make-bed", "--threads", "1", "--out", markers),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(paste0(markers, ".bed"))) {
        stop("plink could not extract the sampled markers (status ", status, ")")
    }
}

group_row <- function(grp, r2) {
    data.table(group = grp, n_markers = nrow(picked), n_unlinked_pairs = length(r2),
               critical_r2 = round(ld_critical_r2(r2, QUANTILE), 6),
               mean_r2 = if (length(r2)) round(mean(r2), 6) else NA_real_,
               quantile = QUANTILE, markers_per_chr = per_chr, seed = SEED)
}

rows <- list()
for (i in seq_len(nrow(manifest))) {
    grp  <- manifest$group[i]
    stem <- manifest$list_file[i]
    if (n_chr < 2) {
        message("WARNING: one chromosome, no unlinked pairs: critical r2 is NA for ", grp)
        rows[[grp]] <- group_row(grp, numeric(0))
        next
    }
    keep_args <- character(0)
    if (grp != "All") {
        smp <- readLines(file.path(SAMPLE_LISTS_DIR, paste0(stem, ".txt")))
        keep <- file.path(tmp, paste0(stem, ".keep"))
        fwrite(data.table("0", smp), keep, sep = "\t", col.names = FALSE)
        keep_args <- c("--keep", keep)
    }
    out <- file.path(tmp, stem)
    status <- system2("plink", c("--bfile", markers, PLINK_CHR, keep_args,
                                 "--r2", "inter-chr", "--ld-window-r2", "0",
                                 "--threads", "1", "--out", out),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(paste0(out, ".ld"))) {
        stop("plink --r2 inter-chr failed for group ", grp, " (status ", status, ")")
    }
    ld <- fread(paste0(out, ".ld"), colClasses = c(CHR_A = "character", CHR_B = "character"))
    ld[, R2 := suppressWarnings(as.numeric(R2))]   # PLINK writes "nan" for a monomorphic SNP
    r2 <- ld[CHR_A != CHR_B & is.finite(R2), R2]
    rows[[grp]] <- group_row(grp, r2)
    message("INFO: ", grp, ": ", length(r2), " unlinked pairs, critical r2 = ",
            signif(rows[[grp]]$critical_r2, 3))
}

unlink(tmp, recursive = TRUE)
fwrite(rbindlist(rows), OUT_TABLE, sep = "\t")
message("INFO: written ", OUT_TABLE)
