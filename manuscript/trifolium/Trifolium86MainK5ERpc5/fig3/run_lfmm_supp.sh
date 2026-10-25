#!/usr/bin/env bash
# LFMM supplement at K 5 (plan ~/.claude/plans/fancy-scribbling-garden.md, Phase 5): rank the 14 distinct predictors by
# LFMM-significant SNPs in repeated imputed vectors (fig3/lfmm_stripe_rank.R), then draw the top 4 Manhattans with the
# derived pipeline script fig3/plot_manhattan_ms.R (same inputs as the pipeline's manhattan rule for LFMM K 5).
# Output: manuscript_fig3/draft/lf_k5_manhattan.{svg,png} + _counts.tsv, lf_k5_stripe_rank.tsv
# Run INSIDE the image: bash /pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/run_lfmm_supp.sh
set -euo pipefail
A=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5
M=/pipeline/Trifolium86MainK5PC5_results/GEA/tables/methods/LFMM
TOP=$(Rscript $A/fig3/lfmm_stripe_rank.R | tail -1)
[ -n "$TOP" ] || { echo "ERROR: no predictors ranked"; exit 1; }
echo "=== drawing $TOP"
Rscript $A/fig3/plot_manhattan_ms.R $M/LFMM_pvalues_K5.tsv bonf_0.05 5 LFMM "$TOP" \
    $A/manuscript_fig3/draft/lf_k5_manhattan $M/LFMM_pvalues_K5_sig_snps_bonf_0.05.tsv 180x150
echo "=== done"
