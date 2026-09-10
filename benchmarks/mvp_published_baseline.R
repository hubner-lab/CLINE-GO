#!/usr/bin/env Rscript
# mvp_published_baseline.R -- extract Lotterhos 2023's OWN published detection metrics for the
# 14 benchmark seeds into one lookup table, so comparing a future CLINE-GO run against the
# original study is a join, not a re-derivation.
#
# Source: data/mvp/selection/summary_20220428_20220726.csv (BCO-DMO doi:10.26008/1912/bco-dmo.889769.1),
# one row per seed. Column definitions: data/mvp/selection/summary_metadata.md.
#
# Output: benchmarks/mvp_published_baseline.tsv
#
# ---------------------------------------------------------------------------------------
# COMPARABILITY: their `*_neutSNPs` variant == our convention.
#
# They define *_neutSNPs as "including only causal loci and neutral loci not affected by
# selection (any non-causal loci that arose on the half of the genome affected by selection
# was excluded)". Our lib_detection.R counts ONLY background_neutral as false positives
# (`is_bg <- ranked %in% truth$bg_keys`); linked_neutral loci are ranked but excluded from the
# precision denominator. Verified numerically equal on seed 1232548: removing linked-neutral
# rows entirely and skipping them in the denominator both give AP = 0.2027 / 0.3938.
#
# => ALWAYS compare our AUC-PR against *_neutSNPs, never against *_allSNPs.
# ---------------------------------------------------------------------------------------
#
# KNOWN DEFECT IN THE PUBLISHED TABLE -- see the `lfmm_temp_aucpr_suspect` column.
# `LEA3_2_lfmm2_AUCPR_temp_neutSNPs` is NOT reproducible from the shipped temperature p-value
# column (`LEA3.2_lfmm2_mlog10P_tempenv`): re-scoring those p-values gives 0.2027 for seed
# 1232548 against a published 0.0783. It IS reproducible from the SALINITY p-value column
# (0.0791), within the same +-0.001..0.007 offset that every other published cell shows.
# Consistent with a column swap in that one metric. Not confirmed against the authors'
# MVP-NonClinalAF source, so treat as "verify before citing" -- and for comparison purposes
# use the re-scored value from benchmarks/mvp_score_published.R instead of this cell.

suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
SUMM <- file.path(ROOT, "data/mvp/selection/summary_20220428_20220726.csv")
MAN  <- file.path(ROOT, "benchmarks/mvp_seeds.tsv")
OUT  <- file.path(ROOT, "benchmarks/mvp_published_baseline.tsv")

s <- fread(SUMM)
m <- fread(MAN)
# meanFst exists in BOTH tables (identical values) — take it from the summary side only, so
# the merge produces no .x/.y suffixes.
d <- merge(m[, .(seed, arm, arch_level, r2_pc1_temp, r2_pc1_sal,
                 n_snps, n_causal_maf01, n_causal_temp, n_causal_sal, K_authors, k_best)],
           s, by = "seed")
stopifnot(nrow(d) == nrow(m))

out <- d[, .(
  seed, arm, arch_level, k_best,
  K_authors,                      # their K -- what their LFMM/PCA used
  meanFst, r2_pc1_temp, r2_pc1_sal,
  n_snps, n_causal_temp, n_causal_sal, n_causal_any = n_causal_maf01,

  # ---- LFMM (LEA3.2 lfmm2) -- the same method our pipeline runs ----
  pub_lfmm_aucpr_temp_neut = LEA3_2_lfmm2_AUCPR_temp_neutSNPs,
  pub_lfmm_aucpr_sal_neut  = LEA3_2_lfmm2_AUCPR_sal_neutSNPs,
  pub_lfmm_aucpr_temp_all  = LEA3_2_lfmm2_AUCPR_temp_allSNPs,
  pub_lfmm_aucpr_sal_all   = LEA3_2_lfmm2_AUCPR_sal_allSNPs,
  pub_lfmm_tpr_temp        = LEA3_2_lfmm2_TPR_temp,
  pub_lfmm_tpr_sal         = LEA3_2_lfmm2_TPR_sal,
  pub_lfmm_fdr_temp_neut   = LEA3_2_lfmm2_FDR_neutSNPs_temp,
  pub_lfmm_fdr_sal_neut    = LEA3_2_lfmm2_FDR_neutSNPs_sal,

  # ---- RDA (rdadapt), with and without Condition(PC1+PC2) ----
  pub_rda_aucpr_neut       = RDA_AUCPR_neutSNPs,
  pub_rda_aucpr_neut_corr  = RDA_AUCPR_neutSNPs_corr,
  pub_rda_aucpr_all        = RDA_AUCPR_allSNPs,
  pub_rda_tpr              = RDA_TPR,
  pub_rda_tpr_corr         = RDA_TPR_corr,
  pub_rda_fdr_neut         = RDA_FDR_neutSNPs,
  pub_rda_fdr_neut_corr    = RDA_FDR_neutSNPs_corr,

  # ---- Kendall correlation of deme allele freq ~ deme env (their simplest baseline) ----
  pub_cor_aucpr_temp_neut  = cor_AUCPR_temp_neutSNPs,
  pub_cor_aucpr_sal_neut   = cor_AUCPR_sal_neutSNPs,
  pub_cor_tpr_temp         = cor_TPR_temp,
  pub_cor_tpr_sal          = cor_TPR_sal,

  # ---- their phenotype GWAS -- reference for OUR mode=gwas arm ----
  pub_gwas_tpr_temp        = gwas_TPR_temp,
  pub_gwas_tpr_sal         = gwas_TPR_sal,
  pub_gwas_fdr_temp        = gwas_FDR_temp_neutbase,
  pub_gwas_fdr_sal         = gwas_FDR_sal_neutbase
)]
out[, lfmm_temp_aucpr_suspect := TRUE]   # see header note; applies to every seed

setorder(out, arm, -meanFst)
fwrite(out, OUT, sep = "\t")
message("INFO: wrote ", OUT, " -- ", nrow(out), " seeds x ", ncol(out), " columns")
print(out[, .(seed, arm, K_authors, pub_lfmm_aucpr_sal_neut, pub_rda_aucpr_neut,
              pub_lfmm_tpr_temp, pub_lfmm_tpr_sal, pub_gwas_tpr_temp, pub_gwas_fdr_temp)])
