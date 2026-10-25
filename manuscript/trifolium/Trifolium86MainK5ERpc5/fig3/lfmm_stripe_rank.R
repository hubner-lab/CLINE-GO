#!/usr/bin/env Rscript
# Which LFMM Manhattans go into the LFMM supplement (plan ~/.claude/plans/fancy-scribbling-garden.md, Phase 5; user
# review: "K3 remove at all. Analyse k5, keep only several Manhattan plots which has the stripes"). Project
# Trifolium86MainK5PC5 (LFMM at K 5, 5 PCs; LFMM is not run in Trifolium86MainK5ERpc5).
# Rule: per predictor, count its LFMM-significant SNPs (pipeline sig table, bonf 0.05) that sit in a repeated imputed
#   genotype vector, i.e. pattern_n >= STRIPE in snp_stripe_imputation.tsv, with STRIPE = 10 as defined at
#   ANALYSIS/Trifolium86Main/imputation_check/stripe_imputation.py:54. Rank by that count (descending); ties by the
#   predictor's total significant SNPs (descending), then by name. bio_11 = bio_8 and bio_19 = bio_16 at every plant,
#   so the duplicates are dropped before ranking (14 distinct predictors, as in the rest of the paper). Top 4 drawn.
# Writes manuscript_fig3/draft/lf_k5_stripe_rank.tsv; prints the top predictors as one comma-separated line (last line).
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/lfmm_stripe_rank.R
suppressPackageStartupMessages(library(data.table))
STRIPE <- 10   # stripe_imputation.py:54
N_TOP  <- 4
DUP    <- c("bio_11", "bio_19")
ST  <- "/pipeline/ANALYSIS/Trifolium86Main/imputation_check/tables_Trifolium86MainK5PC5/snp_stripe_imputation.tsv"
SIG <- "/pipeline/Trifolium86MainK5PC5_results/GEA/tables/methods/LFMM/LFMM_pvalues_K5_sig_snps_bonf_0.05.tsv"
OUT <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/manuscript_fig3/draft/lf_k5_stripe_rank.tsv"

st  <- fread(ST, select = c("SNPID", "pattern_n", "sig_LFMM"), colClasses = list(character = "SNPID"))
sig <- fread(SIG, select = c("SNPID", "trait"), colClasses = list(character = "SNPID"))
# staleness guard: the stripe table must describe the sig table on disk now
stopifnot(setequal(st[sig_LFMM == 1, SNPID], unique(sig$SNPID)))
sig <- merge(sig[!trait %in% DUP], st[, .(SNPID, pattern_n)], by = "SNPID")
rk <- sig[, .(n_sig = .N, n_sig_in_stripe = sum(pattern_n >= STRIPE)), by = .(predictor = trait)]
rk[, share_in_stripe := round(n_sig_in_stripe / n_sig, 3)]
setorder(rk, -n_sig_in_stripe, -n_sig, predictor)
rk[, `:=`(rank = .I, drawn = .I <= N_TOP,
          rule = sprintf("sig SNPs with pattern_n >= %d, desc; ties: n_sig desc, name", STRIPE))]
stopifnot(nrow(rk) == 14)
fwrite(rk, OUT, sep = "\t")
print(rk)
cat(paste(rk[drawn == TRUE, predictor], collapse = ","), "\n", sep = "")
