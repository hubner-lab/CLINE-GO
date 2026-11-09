#!/usr/bin/env Rscript
# Background LD of unlinked markers, per LD-decay group (Breseghello & Sorrells 2006,
# as applied by Laido et al. 2014): r^2 between markers on DIFFERENT chromosomes, for
# a seeded sample of markers per chromosome (Vos et al. 2017 sampled 50), its mean and
# its critical value, the `quantile` of sqrt(r^2), squared. PLINK 1.9 --r2 inter-chr,
# with the chromosome set and missing-call handling of ld_decay_pairs.R, so the
# background and the decay curve come from one estimator. The contigs in
# LDdecay.exclude_chromosomes give no markers, as they give the curve no pairs.
#
# The sample is drawn DRAWS times (LDdecay.background_draws), draw b on seed
# SEED + b - 1, so draw 1 is the single draw that DRAWS = 1 gives. OUT_TABLE has one
# row per group: n_unlinked_pairs, critical_r2 and mean_r2 are medians over draws
# (ld_summarise_draws()), critical_r2_lo and critical_r2_hi the 2.5 and 97.5 %
# quantiles of the per-draw critical r^2. When no chromosome has more SNPs than are
# drawn from it, every draw is the same and the interval is a point. DRAWS_TABLE (or
# "NULL") gets one row per group and draw. The groups' draws run in parallel, on the
# job's OMP_NUM_THREADS processes.
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
DRAWS            <- as.integer(args[7])  # LDdecay.background_draws
DRAWS_TABLE      <- args[8]  # one row per group and draw (or "NULL")
EXCLUDE          <- args[9]  # LDdecay.exclude_chromosomes, comma-separated, or "NULL"
################
stopifnot(length(args) == 9, MARKERS_PER_CHR >= 2, QUANTILE > 0, QUANTILE < 1,
          !is.na(DRAWS), DRAWS >= 1)
THREADS <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", "1")))
if (is.na(THREADS) || THREADS < 1L) THREADS <- 1L

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

# At most MAX_MARKERS markers per draw (ld_background_markers()), however fragmented the
# assembly: unlinked pairs grow with the square of the marker count.
MAX_MARKERS <- 2000L
draw_sets <- ld_background_draws(pos$chr, pos$pos, MARKERS_PER_CHR, MAX_MARKERS, SEED, DRAWS)
per_chr <- draw_sets[[1]]$per_chr
n_markers <- nrow(draw_sets[[1]]$markers)
if (per_chr < MARKERS_PER_CHR) {
    message("INFO: ", n_chr, " chromosomes: ", per_chr, " marker(s) per chromosome, ",
            n_markers, " in all (cap ", MAX_MARKERS, ")")
}
write_ranges <- function(markers, path) {
    fwrite(markers[, .(chr, pos, pos, paste0(chr, ":", pos))], path,
           sep = "\t", col.names = FALSE)
    path
}
draw_ranges <- vapply(seq_len(DRAWS), function(b)
    write_ranges(draw_sets[[b]]$markers, file.path(tmp, paste0("draw", b, ".range"))),
    character(1))
all_markers <- unique(rbindlist(lapply(draw_sets, `[[`, "markers")))
message("INFO: ", DRAWS, " draw(s) of ", n_markers, " markers on seeds ", SEED, "-",
        SEED + DRAWS - 1L, ": ", nrow(all_markers), " distinct markers")

# Read the VCF once: every draw's markers into one small PLINK binary set, which each
# group and draw then subsets (a VCF pass per group cost 101 passes per simulated
# replicate).
markers <- file.path(tmp, "markers")
if (n_chr >= 2) {
    status <- system2("plink", c("--vcf", VCF_PATH, "--const-fid", PLINK_CHR,
                                 "--vcf-half-call", "missing", "--extract", "range",
                                 write_ranges(all_markers, file.path(tmp, "markers.range")),
                                 "--make-bed", "--threads", "1", "--out", markers),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(paste0(markers, ".bed"))) {
        stop("plink could not extract the sampled markers (status ", status, ")")
    }
}

keep_args <- lapply(seq_len(nrow(manifest)), function(i) {
    if (manifest$group[i] == "All") return(character(0))
    smp <- readLines(file.path(SAMPLE_LISTS_DIR, paste0(manifest$list_file[i], ".txt")))
    keep <- file.path(tmp, paste0(manifest$list_file[i], ".keep"))
    fwrite(data.table("0", smp), keep, sep = "\t", col.names = FALSE)
    c("--keep", keep)
})

draw_row <- function(i, b, r2) {
    data.table(group = manifest$group[i], draw = b, seed = SEED + b - 1L,
               n_markers = n_markers, n_unlinked_pairs = length(r2),
               critical_r2 = round(ld_critical_r2(r2, QUANTILE), 6),
               mean_r2 = if (length(r2)) round(mean(r2), 6) else NA_real_)
}

# One group and draw: PLINK --r2 inter-chr on that draw's markers and that group's
# samples, the unlinked pairs' r^2 summarised and the .ld file removed.
run_draw <- function(i, b) {
    out <- file.path(tmp, paste0(manifest$list_file[i], "_draw", b))
    status <- system2("plink", c("--bfile", markers, PLINK_CHR, keep_args[[i]],
                                 "--extract", "range", draw_ranges[b],
                                 "--r2", "inter-chr", "--ld-window-r2", "0",
                                 "--threads", "1", "--out", out),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(paste0(out, ".ld"))) {
        stop("plink --r2 inter-chr failed for group ", manifest$group[i], ", draw ", b,
             " (status ", status, ")")
    }
    ld <- fread(paste0(out, ".ld"), colClasses = c(CHR_A = "character", CHR_B = "character"))
    unlink(paste0(out, c(".ld", ".log", ".nosex")))
    ld[, R2 := suppressWarnings(as.numeric(R2))]   # PLINK writes "nan" for a monomorphic SNP
    draw_row(i, b, ld[CHR_A != CHR_B & is.finite(R2), R2])
}

jobs <- CJ(i = seq_len(nrow(manifest)), b = seq_len(DRAWS))
if (n_chr < 2) {
    message("WARNING: one chromosome, no unlinked pairs: critical r2 is NA for ",
            paste(manifest$group, collapse = ", "))
    per_draw <- rbindlist(Map(function(i, b) draw_row(i, b, numeric(0)), jobs$i, jobs$b))
} else {
    message("INFO: ", nrow(jobs), " PLINK runs (", nrow(manifest), " group(s) x ", DRAWS,
            " draw(s)) on ", THREADS, " process(es)")
    res <- parallel::mclapply(seq_len(nrow(jobs)), function(j) run_draw(jobs$i[j], jobs$b[j]),
                              mc.cores = THREADS)
    failed <- vapply(res, function(r) inherits(r, "try-error") || !is.data.table(r), logical(1))
    if (any(failed)) {
        stop("background draw failed: ", as.character(res[[which(failed)[1]]]))
    }
    per_draw <- rbindlist(res)
}

summary_rows <- per_draw[, {
    s <- ld_summarise_draws(critical_r2, mean_r2, n_unlinked_pairs)
    .(n_markers = n_markers[1], n_unlinked_pairs = as.integer(round(s$n_unlinked_pairs)),
      critical_r2 = round(s$critical_r2, 6), mean_r2 = round(s$mean_r2, 6),
      quantile = QUANTILE, markers_per_chr = per_chr, seed = SEED, n_draws = s$n_draws,
      critical_r2_lo = round(s$critical_r2_lo, 6), critical_r2_hi = round(s$critical_r2_hi, 6))
}, by = group]
for (k in seq_len(nrow(summary_rows))) {
    r <- summary_rows[k]
    message("INFO: ", r$group, ": ", r$n_unlinked_pairs, " unlinked pairs per draw, critical r2 = ",
            signif(r$critical_r2, 3), " (", signif(r$critical_r2_lo, 3), "-",
            signif(r$critical_r2_hi, 3), " over ", r$n_draws, " draw(s))")
}

unlink(tmp, recursive = TRUE)
fwrite(summary_rows, OUT_TABLE, sep = "\t")
message("INFO: written ", OUT_TABLE)
if (DRAWS_TABLE != "NULL") {
    dir.create(dirname(DRAWS_TABLE), recursive = TRUE, showWarnings = FALSE)
    fwrite(per_draw, DRAWS_TABLE, sep = "\t")
    message("INFO: written ", DRAWS_TABLE)
}
