#!/usr/bin/env bash
# The 15 candidate regions in topr, one figure per region, for Fig 3 v2 (plan ~/.claude/plans/fancy-scribbling-garden.md,
# Phase 3). Same calls and inputs as regions/figures/run_regionplots.sh (the region-QC run), drawn by
# fig3/plot_regionplot_topr_ms.R at 180 x 79 mm. Optional first argument: one region order number (test run).
# Output: manuscript_fig3/draft/rg_topr/<order>/regionplot_custom_*.{svg,png}, copied to draft/rgt_<NN>_<region>.{svg,png}
# Run INSIDE the image: bash /pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/run_regions_topr.sh [order]
set -euo pipefail
A=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5
F=$A/regions/figures
D=$A/manuscript_fig3/draft
M=/pipeline/Trifolium86MainK5ERpc5_results/GEA/tables/methods
ASSOC="EMMAX:bonf_0.05:$M/EMMAX/EMMAX_pvalues_K5.tsv,RDA:qval_0.05:$M/RDA/RDA_pvalues_K5.tsv"
ONLY=${1:-}
tail -n +2 "$F/input/regionplot_calls.tsv" | while IFS=$'\t' read -r order region custom traits vlines nser; do
    [ -n "$ONLY" ] && [ "$order" != "$ONLY" ] && continue
    out="$D/rg_topr/$order/"
    mkdir -p "$out"
    echo "=== $order $region $custom n_series=$nser"
    Rscript "$A/fig3/plot_regionplot_topr_ms.R" NULL "$F/input/topr_labelled.tsv" "$ASSOC" 0 all "$out" \
        "$custom" "$traits" NULL snp "$vlines" 180x79
    stem=$(ls "$out"regionplot_custom_*.svg | head -1); stem=${stem%.svg}
    [ -s "$stem.svg" ] && [ -s "$stem.png" ] || { echo "ERROR: no figure for $region"; exit 1; }
    nn=$(printf '%02d' "$order")
    cp "$stem.svg" "$D/rgt_${nn}_${region}.svg"; cp "$stem.png" "$D/rgt_${nn}_${region}.png"
done
echo "=== done"
