#!/usr/bin/env bash
# Offset figures recoloured to minou for the Fig 3 review (plan scalable-doodling-allen.md, Phase 2). Writes ONLY into
# fig3/offset_minou/ (symlinked inputs, copied paper_pack/ + site_maps/), so the approved v2 pack is untouched.
# Phase 5 of fancy-scribbling-garden.md (2026-10-16): site labels = short names 1..11 (SITE_NAMES), S6 = panel (a)
# only (KCHOICE_PANELS=a). S4 (fig_s_gf_bio5) is still drawn but dropped from the pack. Defaults (no env) are
# byte-identical to v1 / v2 (checked against sha256 manifests).
# Run inside the image: bash /pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/run_offset_minou.sh
set -euo pipefail
cd /pipeline
export AN=ANALYSIS/Trifolium86MainK5ERpc5/fig3/offset_minou PROJ=Trifolium86MainK5ERpc5 PALETTE=minou
export SITE_NAMES=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.tsv KCHOICE_PANELS=a DELTA_CORR=1 NO_LINEAR=1 \
    MAIN_MS=1 SUPP_MS=1   # Phase 5 (2026-10-17): Fig 5 + S15-S19 also at final mm size as *_ms.{svg,png}
export MAIN6_MS=1   # 2026-10-18: manuscript Figure 5 = offset maps (a-c) over fig_main_ms (d-f), fig5_offset_ms.{svg,png}
export SUPP2_MS=1   # 2026-10-25 (supplement build): s_k_choice, s_site_rankings, s_predictor_change, s_random_sets
V1=ANALYSIS/Trifolium86MainK5/maladaptationER
Rscript $V1/build_site_maps.R
Rscript $V1/build_paper_pack.R
Rscript ANALYSIS/Trifolium86MainK5ERpc5/maladaptation/fig_k_choice.R
Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/plot_random_climate_ms.R   # 2026-10-25: s_random_climate (B8)
echo done
