#!/usr/bin/env bash
# Step 12 figure set for one SS-Clines block (or the pooled corpus).
#
#   work/block_report.sh <cohort> <added_tags_csv> <n_expect> <n_per_arch> [offset_dir_name]
#   work/block_report.sh ssclines_b3 ssclines_nequal_mconst 120 40
#   work/block_report.sh ssclines_pooled tag1,tag2,tag3,tag4,tag5 600 200 offset12_ssclines_pooled
#
# Generalised from work/block_report_b2.sh. Two differences from that copy:
#   - image is cline-go:latest (the drivers' IMAGE default), not adaptogene:latest
#   - PIPELINE_ROOT is exported into the container
#
# MVP_ADDED is MANDATORY: the manifest holds 120 rows per block under ONE
# arm == "primary_ssclines", so filtering on arm alone silently pools every block.
set -u
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

COHORT="${1:?cohort}"
TAGS="${2:?added tags csv}"
NEXP="${3:?n_expect}"
NARCH="${4:?n_per_arch}"
ODIR="${5:-offset12_${COHORT}}"

FIGS_HOST="benchmarks/mvp_eval/figures_${COHORT}"
FIGS=/pipeline/${FIGS_HOST}
DK=(nix shell nixpkgs#docker-client -c docker)
mkdir -p "$FIGS_HOST" work/blockrep_${COHORT}

echo "=== step 12 figure set: cohort=$COHORT offset_dir=$ODIR n=$NEXP per_arch=$NARCH  $(date -Is)"
fail=0
for s in mvp_oracle_stats.R mvp_main_figure_v2.R mvp_block_report.R \
         mvp_arch_panel.R mvp_method_panel.R mvp_dist_panel.R mvp_regime_panel.R \
         mvp_slope_panel.R mvp_single_method_panels.R mvp_main_figure.R \
         mvp_oracle_figure.R mvp_method_figure.R mvp_absolute_figures.R \
         mvp_manuscript_figures.R; do
    echo "--- $s  $(date -Is)"
    "${DK[@]}" run --rm --name "mvp-fig-${COHORT}-${s%.R}" --user "$(id -u):$(id -g)" \
        -e USER=adaptogene -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 \
        --cpus=8 --memory=48g \
        -e OFFSET_DIR="$ODIR" -e MVP_ARM=primary_ssclines \
        -e MVP_ADDED="$TAGS" -e MVP_N_EXPECT="$NEXP" -e MVP_N_PER_ARCH="$NARCH" \
        -e FIG_OUT="$FIGS" -e STATS_DIR="$FIGS" \
        -v "$PWD":/pipeline cline-go:latest \
        Rscript /pipeline/benchmarks/$s > "work/blockrep_${COHORT}/${s%.R}.log" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "    FAILED rc=$rc  (work/blockrep_${COHORT}/${s%.R}.log)"
        tail -5 "work/blockrep_${COHORT}/${s%.R}.log" | sed 's/^/      /'
        fail=$((fail+1))
    else
        echo "    ok"
    fi
done
echo "=== figure set done, $fail failure(s)  $(date -Is)"
echo "files: $(ls "$FIGS_HOST" | wc -l)"

# 4 of 29 manuscript figures skip by design (B6, B6b, D8, T5 -- detection arm sweep07).
echo "--- manuscript-figure skips (expect B6, B6b, D8, T5):"
grep -iE "skip" "work/blockrep_${COHORT}/mvp_manuscript_figures.log" 2>/dev/null | head -8

# Cross-arm pooling scan: a TSV carrying more distinct seeds than this corpus should have
# means a script absorbed another arm. +2 allows the degenerate controls several scripts keep.
echo "--- cross-arm pooling scan (must print nothing):"
for f in "$FIGS_HOST"/*.tsv; do
    [[ -f "$f" ]] || continue
    head -1 "$f" | grep -q 'seed' || continue
    n=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="seed")c=i;next}c{print $c}' "$f" | sort -u | wc -l)
    (( n > NEXP + 2 )) && echo "POOLING  $f  distinct_seeds=$n"
done
echo "--- scan done"

# The SHORT report: exactly four figures, no prose.
echo "--- short report figures:"
for stem in MAIN_FIGURE_simulation panelA_composition panelB_match_oracle G3_dist_transposed; do
    hit=$(ls "$FIGS_HOST"/${stem}.* 2>/dev/null | head -2 | paste -sd' ')
    echo "  ${stem}: ${hit:-MISSING}"
done
exit $fail
