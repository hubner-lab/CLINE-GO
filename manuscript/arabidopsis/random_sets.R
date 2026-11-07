#!/usr/bin/env Rscript
# random_sets.R: random SNP sets matched on MAF to one real set, the draws of offset_eval.R's level 3.
#
#   Rscript manuscript/arabidopsis/random_sets.R <project> <real set> <draws> [--seed N]
#
# Phase 15, the user's choice on 2026-11-07: draws for the union set only, matched on MAF. Each draw
# holds the real set's SNP count, taken without replacement from the tested SNPs outside the real
# set, with the real set's count in every MAF bin of width 0.01, so a draw differs from the real set
# in association and not in allele frequency. Draw d uses seed N + d (default N = 20261107), so any
# one draw can be rebuilt alone.
#
# MAF comes from the matrix every offset script reads,
# <project>_results/_intermediate/climate_subset/lfmm_imp_full_climate.lfmm (imputed, climate-valid
# samples), its columns named chr:pos by the one *imp_maf*.vcfsnp under _work/, the file the offset
# rules take as W['vcfsnp_imp_full']. SNPs in that _work folder's .removed list stay out of the pool,
# as the offset scripts leave them out.
#
# Writes <project>_results/_intermediate/snp_sets/random.<real set>.<d>/ (selected_snps.tsv,
# snp_set.yaml): mode=maladaptation's `snp_sets: all` takes every stored set, and offset_eval.R reads
# random.<set>.<d> as a draw of <set>. Also writes ANALYSIS/Arabidopsis1135/random_sets/<real set>/
# maf_match.tsv, one row per draw. Refuses to touch an existing draw directory. Run from the repo
# root, in the image.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})
setDTthreads(8L)

args <- commandArgs(trailingOnly = TRUE)
usage <- "usage: random_sets.R <project> <real set> <draws> [--seed N]"
if (!length(args) %in% c(3L, 5L)) stop(usage)
project <- args[1]
real_set <- args[2]
n_draws <- as.integer(args[3])
seed <- 20261107L
if (length(args) == 5L) {
  if (args[4] != "--seed") stop(usage)
  seed <- as.integer(args[5])
}
stopifnot(!is.na(n_draws), n_draws >= 1L, !is.na(seed))
bin_width <- 0.01

results <- here(paste0(project, "_results"))
store <- file.path(results, "_intermediate", "snp_sets")
real_file <- file.path(store, real_set, "selected_snps.tsv")
matrix_file <- file.path(results, "_intermediate", "climate_subset", "lfmm_imp_full_climate.lfmm")
vcfsnp_file <- Sys.glob(file.path(results, "_work", "*", "*imp_maf*.vcfsnp"))
removed_file <- Sys.glob(file.path(results, "_work", "*", "*.removed"))
if (length(vcfsnp_file) != 1L) stop("Expected one *imp_maf*.vcfsnp under _work/, found ", length(vcfsnp_file))
if (length(removed_file) != 1L) stop("Expected one .removed under _work/, found ", length(removed_file))
stopifnot(file.exists(real_file), file.exists(matrix_file))
out_dir <- here("ANALYSIS", "Arabidopsis1135", "random_sets", real_set)
draw_names <- sprintf("random.%s.%d", real_set, seq_len(n_draws))
taken <- draw_names[dir.exists(file.path(store, draw_names))]
if (length(taken) > 0L) stop("Draw directories exist already: ", paste(taken, collapse = ", "))

snps <- fread(vcfsnp_file, header = FALSE, select = 1:2, col.names = c("chr", "pos"),
              colClasses = list(character = 1))
snps[, SNPID := paste0(chr, ":", pos)]
message("INFO: reading ", matrix_file)
genotypes <- as.matrix(fread(matrix_file, header = FALSE))
if (ncol(genotypes) != nrow(snps)) {
  stop(ncol(genotypes), " matrix columns against ", nrow(snps), " SNPs in ", vcfsnp_file)
}
if (any(genotypes == 9L)) stop("The matrix holds LFMM's missing code 9: not an imputed matrix")
n_samples <- nrow(genotypes)
frequency <- colMeans(genotypes) / 2
rm(genotypes)
snps[, maf := pmin(frequency, 1 - frequency)]
snps[, bin := as.integer(floor(maf / bin_width + 1e-9))]
message("INFO: ", nrow(snps), " SNPs by ", n_samples, " samples; MAF ", signif(min(snps$maf), 3),
        " to ", signif(max(snps$maf), 3))

removed <- if (file.size(removed_file) > 0) {
  fread(removed_file, header = FALSE, select = 1:2)[, paste0(V1, ":", V2)]
} else character()
real <- fread(real_file, colClasses = c(SNPID = "character"))$SNPID
missing_real <- setdiff(real, snps$SNPID)
if (length(missing_real) > 0L) stop(length(missing_real), " SNPs of ", real_set, " are not in the matrix")
target <- snps[SNPID %in% real, .N, by = bin]
pool <- snps[!SNPID %in% c(real, removed)]
short <- target[pool[, .(available = .N), by = bin], on = "bin"][N > available]
if (nrow(short) > 0L) stop("MAF bins with fewer pool SNPs than the real set needs: ", paste(short$bin, collapse = ", "))
message("INFO: ", real_set, " holds ", length(real), " SNPs in ", nrow(target), " MAF bins; pool ",
        nrow(pool), " SNPs")

real_maf <- snps[SNPID %in% real, maf]
quality <- rbindlist(lapply(seq_len(n_draws), function(d) {
  set.seed(seed + d)
  drawn <- pool[target, on = "bin", .SD[sample.int(.N, i.N)], by = .EACHI]
  drawn <- drawn[order(match(chr, unique(snps$chr)), pos)]
  stopifnot(nrow(drawn) == length(real), !anyDuplicated(drawn$SNPID), !any(drawn$SNPID %in% real))
  set_dir <- file.path(store, draw_names[d])
  dir.create(set_dir)
  fwrite(drawn[, .(SNPID, chr, pos, min_pvalue = NA_real_)], file.path(set_dir, "selected_snps.tsv"),
         sep = "\t", na = "NA")
  writeLines(c(
    paste0("name: ", draw_names[d]),
    "origin: manual",
    "source_module: GEA",
    "regime: snp",
    "strategy: random_maf_matched",
    "clumping_distance: ''",
    paste0("matched_to: ", real_set),
    paste0("draw: ", d),
    paste0("seed: ", seed + d),
    paste0("maf_bin_width: ", bin_width),
    paste0("pool: tested SNPs outside ", real_set, " and the removed list (", nrow(pool), ")"),
    "script: manuscript/arabidopsis/random_sets.R"
  ), file.path(set_dir, "snp_set.yaml"))
  ks <- suppressWarnings(ks.test(drawn$maf, real_maf)$statistic)
  quantiles <- quantile(drawn$maf, c(0.25, 0.5, 0.75), names = FALSE)
  data.table(set = draw_names[d], seed = seed + d, n_snps = nrow(drawn), maf_mean = mean(drawn$maf),
             maf_q1 = quantiles[1], maf_median = quantiles[2], maf_q3 = quantiles[3],
             ks_distance_to_real = unname(ks), chromosomes = uniqueN(drawn$chr))
}))
real_quantiles <- quantile(real_maf, c(0.25, 0.5, 0.75), names = FALSE)
quality <- rbind(data.table(set = real_set, seed = NA_integer_, n_snps = length(real_maf),
                            maf_mean = mean(real_maf), maf_q1 = real_quantiles[1],
                            maf_median = real_quantiles[2], maf_q3 = real_quantiles[3],
                            ks_distance_to_real = 0, chromosomes = uniqueN(snps[SNPID %in% real, chr])),
                 quality)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(quality, file.path(out_dir, "maf_match.tsv"), sep = "\t")
print(quality[, .(set, n_snps, maf_mean = round(maf_mean, 4), maf_median = round(maf_median, 4),
                  ks = round(ks_distance_to_real, 4), chromosomes)], nrows = 100)
message("Wrote ", n_draws, " draws to ", store, " and ", file.path(out_dir, "maf_match.tsv"))
