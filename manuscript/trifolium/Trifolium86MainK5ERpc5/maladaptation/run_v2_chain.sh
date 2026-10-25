#!/usr/bin/env bash
# v2 (PC5) offset chain: the v1 scripts in ANALYSIS/Trifolium86MainK5/maladaptationER/ run against the PC5 project via
# env AN / PROJ (plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 4). R steps only; the Python steps
# (random_consensus_test.py, kendall_check.py) run on the host between gf/geo nulls and build_paper_pack.R.
# Run inside the image from /pipeline: bash ANALYSIS/Trifolium86MainK5ERpc5/maladaptation/run_v2_chain.sh <stage>
#   stage nulls: project_gf_random.sh arm3, build_site_maps.R, gf_random_null.R, geo_rda_random_null.R, offset_attribution.R
#   stage pack:  build_paper_pack.R, fig_k_choice.R
set -euo pipefail
cd /pipeline
export AN=ANALYSIS/Trifolium86MainK5ERpc5/maladaptation PROJ=Trifolium86MainK5ERpc5 CORES=${CORES:-160}
V1=ANALYSIS/Trifolium86MainK5/maladaptationER
case "${1:?stage: nulls | pack}" in
  nulls)
    bash $V1/project_gf_random.sh arm3
    Rscript $V1/build_site_maps.R
    Rscript $V1/gf_random_null.R
    Rscript $V1/geo_rda_random_null.R
    Rscript $V1/offset_attribution.R ;;
  pack)
    Rscript $V1/build_paper_pack.R
    Rscript $AN/fig_k_choice.R ;;
esac
echo "stage $1 done"
