#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_set_profile.R -- Job C of the landscapes pack: what the per-method SNP sets pick, per landscape, on all
# 1,800 replicates. work/mvp_land/lfmm_probe/set_profile.R (stratified sample of 90) with three changes: every
# replicate instead of a draw, 24 cores instead of 48, a stop if a landscape does not hold 600 replicates.
# Per replicate x set (LFMM, EMMAX, RDA, the causal loci): composition of the set (causal temp / sal) and each SNP's
# deme-level allele-frequency cline: |r| with temperature, with salinity, with neutral PC1, and its sd.
# PC1 here = first principal component of the deme x SNP allele frequencies of the SNPs neutral for both traits;
# r2_pc1_temp / r2_pc1_sal here = Pearson r^2 of that PC1 with the predictor across demes.
# Usage: Rscript benchmarks/mvp_land_set_profile.R <out.tsv>
# =============================================================================
suppressMessages({library(data.table); library(parallel)})
args <- commandArgs(TRUE); OUT <- args[1]
stopifnot(length(args) == 1L)
setDTthreads(1)
PS <- fread("/pipeline/benchmarks/mvp_eval/landscape_compare/per_seed.tsv", colClasses = c(seed = "character"))
SEEDS <- unique(PS[family == "TP", .(landscape, seed, arch)])
stopifnot(all(SEEDS[, .N, by = landscape]$N == 600L), nrow(SEEDS) == 1800L)
SETS <- c(LFMM = "gea.solo_lfmm", EMMAX = "gea.solo_emmax", RDA = "gea.solo_rda", truth = "truth")
tree <- function(L, s) if (L == "SS-Clines")
  sprintf("/pipeline/_archive/mvp_ssclines/results/MVP%shd_results", s) else sprintf("/pipeline/MVP%shd_results", s)

one <- function(i) {
  L <- SEEDS$landscape[i]; s <- SEEDS$seed[i]; D <- sprintf("/pipeline/data/mvp/MVP%s", s)
  meta <- fread(file.path(D, sprintf("MVP%s_metadata.tsv", s)), colClasses = c(site = "character", sample = "character"))
  env  <- fread(file.path(D, sprintf("MVP%s_env_present.tsv", s)), colClasses = c(site = "character"))
  vcf  <- fread(cmd = sprintf("grep -v '^##' %s", file.path(D, sprintf("MVP%s.vcf", s))), colClasses = "character")
  ids  <- paste0(vcf[["#CHROM"]], ":", vcf$POS)
  G <- as.matrix(vcf[, -(1:9)]); G <- (substr(G, 1, 1) == "1") + (substr(G, 3, 3) == "1")
  site <- meta$site[match(colnames(vcf)[-(1:9)], meta$sample)]
  stopifnot(!anyNA(site))
  demes <- sort(unique(site))
  F <- sapply(demes, function(d) rowMeans(G[, site == d, drop = FALSE]) / 2)   # SNP x deme
  e <- env[match(demes, site)]
  tr <- fread(file.path(D, "loci_freq.tsv"))
  stopifnot(identical(paste0(tr$chr, ":", tr$pos), ids))
  bg <- which(tr$causal_temp == "neutral" & tr$causal_sal == "neutral")
  bgF <- F[bg, ]; bgF <- bgF[apply(bgF, 1, sd) > 0, ]
  pc1 <- prcomp(t(bgF))$x[, 1]
  ok <- apply(F, 1, sd) > 0
  rt <- rs <- rp <- rep(NA_real_, nrow(F))
  rt[ok] <- abs(cor(t(F[ok, ]), e$bio_1)); rs[ok] <- abs(cor(t(F[ok, ]), e$bio_2)); rp[ok] <- abs(cor(t(F[ok, ]), pc1))
  sdf <- apply(F, 1, sd)
  T <- tree(L, s)
  rbindlist(lapply(names(SETS), function(nm) {
    f <- file.path(T, "_intermediate/snp_sets", SETS[[nm]], "selected_snps.tsv")
    x <- fread(f, colClasses = c(SNPID = "character")); j <- match(x$SNPID, ids)
    stopifnot(!anyNA(j))
    data.table(landscape = L, seed = s, arch = SEEDS$arch[i], set = nm, n = length(j),
      sh_ctemp = mean(tr$causal_temp[j] == "causal"), sh_csal = mean(tr$causal_sal[j] == "causal"),
      sh_causal = mean(tr$causal[j]),
      r_temp = median(rt[j], na.rm = TRUE), r_sal = median(rs[j], na.rm = TRUE),
      r_pc1 = median(rp[j], na.rm = TRUE), sd_freq = median(sdf[j]),
      r2_pc1_temp = cor(pc1, e$bio_1)^2, r2_pc1_sal = cor(pc1, e$bio_2)^2, n_bg = nrow(bgF))
  }))
}
RES <- mclapply(seq_len(nrow(SEEDS)), function(i) tryCatch(one(i), error = function(e) {
  message("FAIL ", SEEDS$seed[i], ": ", conditionMessage(e)); NULL }), mc.cores = 24)
R <- rbindlist(RES)
fwrite(R, OUT, sep = "\t")
message("seeds done: ", uniqueN(R$seed), " of ", nrow(SEEDS))
print(R[, lapply(.SD, median), by = .(landscape, set),
        .SDcols = c("n", "sh_ctemp", "sh_csal", "r_temp", "r_sal", "r_pc1", "sd_freq", "r2_pc1_temp", "r2_pc1_sal")][order(landscape, set)], digits = 3)
if (uniqueN(R$seed) != 1800L) stop("only ", uniqueN(R$seed), " replicates profiled")
