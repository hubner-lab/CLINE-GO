#!/usr/bin/env bash
# Phase 5 (SS-Clines re-analysis): manuscript simulation figures on the rebuilt panels.
#
#   step 1  benchmarks/mvp_j16_gallery.R    -> benchmarks/mvp_eval/figures_ssclines_gallery13/
#           the journal-16 compute layer re-pointed at offset13 (per-seed detection + offset
#           tables, stats); regression-tied to mvp_oracle_stats.R's offset13 tables
#   step 1b benchmarks/mvp_panel_rank.R     -> benchmarks/mvp_eval/figures_ssclines_offset13/rank/
#           within-replicate rank of the six panels per engine (rank_long, rank_summary; V1-V5
#           exploration figures); journal 18 Step 18
#   step 1c benchmarks/mvp_rank_reliability.R -> benchmarks/mvp_eval/figures_ssclines_offset13/rank/
#           split-half reliability, engine agreement and variance partition of the
#           within-replicate panel rank (numbers for the manuscript text only, no figure)
#   step 2  benchmarks/mvp_ms_sim_figures.R -> benchmarks/mvp_eval/figures_ssclines_ms/
#           manuscript figures (main: detection plane + within-replicate rank), captions.md,
#           numbers.tsv, index.html + zip; served on :8099 via work/journal/sim_figures
#
# Every arm input is passed explicitly: the scripts' defaults name older arms that still
# exist on disk (offset12_ssclines_pooled, figures_ssclines_pooled_j15, detection600).
#
#   bash work/p5_figures.sh [gallery|rank|reliability|ms|all]  logs: work/p5_logs/<script>.log
#
# Image cline-go:latest (as work/j16_figures.sh). --name is required on this host.
set -u
cd /mnt/data/eugene/ADAPTOGENE || exit 1
WHAT="${1:-all}"
DK=(nix shell nixpkgs#docker-client -c docker)
TAGS=ssclines_nvar_mvar,ssclines_ncline_ns,ssclines_nequal_mconst,ssclines_ncline_ctredge,ssclines_nequal_mbreaks
EV=/pipeline/benchmarks/mvp_eval
mkdir -p work/p5_logs

run() {   # run <script> <cpus> <mem> [extra -e args...]
    local s="$1" cpus="$2" mem="$3"; shift 3
    echo "--- $s  $(date -Is)"
    "${DK[@]}" run --rm --name "mvp-p5-${s%.R}" --user "$(id -u):$(id -g)" \
        -e USER=adaptogene -e PIPELINE_ROOT=/pipeline --cpus="$cpus" --memory="$mem" \
        -e MVP_ARM=primary_ssclines -e MVP_ADDED="$TAGS" -e MVP_N_EXPECT=600 \
        -e MVP_SUBTITLE=0 "$@" \
        -v "$PWD":/pipeline cline-go:latest \
        Rscript "/pipeline/benchmarks/$s" > "work/p5_logs/${s%.R}.log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "    FAILED rc=$rc (work/p5_logs/${s%.R}.log)"; tail -5 "work/p5_logs/${s%.R}.log"
        exit $rc
    fi
    echo "    ok; warnings: $(grep -ci 'warning' "work/p5_logs/${s%.R}.log")"
}

[[ "$WHAT" == gallery || "$WHAT" == all ]] && \
    run mvp_j16_gallery.R 8 48g -e OPENBLAS_NUM_THREADS=4 \
        -e OFFSET_DIR=offset13 -e J15_DIR=$EV/figures_ssclines_offset13 \
        -e AUC_FILE=$EV/remeasure600/rdaunc_rda2x/rank_metrics.tsv \
        -e FIG_OUT=$EV/figures_ssclines_gallery13
[[ "$WHAT" == rank || "$WHAT" == all ]] && \
    run mvp_panel_rank.R 4 16g -e OPENBLAS_NUM_THREADS=1 \
        -e OFFSET_DIR=offset13 -e FIG_OUT=$EV/figures_ssclines_offset13/rank
[[ "$WHAT" == reliability || "$WHAT" == all ]] && \
    run mvp_rank_reliability.R 8 64g -e OPENBLAS_NUM_THREADS=1 \
        -e OFFSET_DIR=offset13 -e FIG_OUT=$EV/figures_ssclines_offset13/rank
[[ "$WHAT" == ms || "$WHAT" == all ]] && \
    run mvp_ms_sim_figures.R 4 16g -e OPENBLAS_NUM_THREADS=2 \
        -e GALLERY_DIR=$EV/figures_ssclines_gallery13 -e OFFSET_DIR=offset13 \
        -e RANK_DIR=$EV/figures_ssclines_offset13/rank -e SIZE_DIR=$EV/figures_ssclines_offset13/size_b1 \
        -e RDA_FIG_DIR=$EV/figures_ssclines_rda -e FIG_OUT=$EV/figures_ssclines_ms

# Cross-arm pooling scan: no emitted table may carry more distinct seeds than the corpus.
echo "--- cross-arm pooling scan (must print nothing):"
for f in benchmarks/mvp_eval/figures_ssclines_gallery13/*.tsv benchmarks/mvp_eval/figures_ssclines_ms/*.tsv; do
    [[ -e "$f" ]] || continue
    col=$(head -1 "$f" | tr '\t' '\n' | grep -nx seed | cut -d: -f1)
    [[ -z "$col" ]] && continue
    n=$(tail -n +2 "$f" | cut -f"$col" | sort -u | wc -l)
    (( n > 600 )) && echo "    $f: $n seeds"
done
echo "=== done $(date -Is)"
