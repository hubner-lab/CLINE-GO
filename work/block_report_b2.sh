#!/usr/bin/env bash
# Step 12 figure set for ssclines_b2. Sequential; each script gets 8 cores.
set -u
COHORT=ssclines_b2
TAG=ssclines_ncline_ns
FIGS=/pipeline/benchmarks/mvp_eval/figures_${COHORT}
DK="nix shell nixpkgs#docker-client -c docker"
cd /mnt/data/eugene/ADAPTOGENE
mkdir -p "benchmarks/mvp_eval/figures_${COHORT}"
fail=0
for s in mvp_oracle_stats.R mvp_main_figure_v2.R mvp_block_report.R \
         mvp_arch_panel.R mvp_method_panel.R mvp_dist_panel.R mvp_regime_panel.R \
         mvp_slope_panel.R mvp_single_method_panels.R mvp_main_figure.R \
         mvp_oracle_figure.R mvp_method_figure.R mvp_absolute_figures.R \
         mvp_manuscript_figures.R; do
    echo "=== $s  $(date -Is)"
    $DK run --rm --name "mvp-fig-${COHORT}-${s%.R}" --user "$(id -u):$(id -g)" \
        -e USER=adaptogene -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=48g \
        -e OFFSET_DIR=offset12_${COHORT} -e MVP_ARM=primary_ssclines \
        -e MVP_ADDED=${TAG} -e MVP_N_EXPECT=120 -e MVP_N_PER_ARCH=40 \
        -e FIG_OUT=$FIGS -e STATS_DIR=$FIGS \
        -v /mnt/data/eugene/ADAPTOGENE:/pipeline adaptogene:latest \
        Rscript /pipeline/benchmarks/$s > "work/blockrep_b2_${s%.R}.log" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then echo "    FAILED rc=$rc  (work/blockrep_b2_${s%.R}.log)"; fail=$((fail+1));
    else echo "    ok"; fi
done
echo "=== figure set done, $fail failure(s)  $(date -Is)"
ls benchmarks/mvp_eval/figures_${COHORT} | wc -l
