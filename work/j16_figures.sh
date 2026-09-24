#!/usr/bin/env bash
# Journal-16 gallery: labelled presentation variants (A detection view, B detection
# scalar, C offset view, D split) for the pooled SS-Clines corpus, 600 replicates.
#
#   step 1  benchmarks/mvp_detection_600.R  -> benchmarks/mvp_eval/detection600/
#           AUC-PR per seed x single method; hard gate: legacy seed 1231288 must
#           reproduce journal 07 (redundancy07) before the corpus is scored
#   step 2  benchmarks/mvp_j16_gallery.R    -> benchmarks/mvp_eval/figures_ssclines_j16_gallery/
#           26 figures + same-stem TSVs; regression-tied to journal 15 delta tables
#
#   bash work/j16_figures.sh [detection|gallery|all]     logs: work/j16_logs/<script>.log
#
# Image cline-go:latest (as work/block_report.sh). --name is required on this host.
set -u
cd /mnt/data/eugene/ADAPTOGENE || exit 1
WHAT="${1:-all}"
DK=(nix shell nixpkgs#docker-client -c docker)
TAGS=ssclines_nvar_mvar,ssclines_ncline_ns,ssclines_nequal_mconst,ssclines_ncline_ctredge,ssclines_nequal_mbreaks
mkdir -p work/j16_logs

run() {   # run <script> <cpus> <mem> [extra -e args...]
    local s="$1" cpus="$2" mem="$3"; shift 3
    echo "--- $s  $(date -Is)"
    "${DK[@]}" run --rm --name "mvp-j16-${s%.R}" --user "$(id -u):$(id -g)" \
        -e USER=adaptogene -e PIPELINE_ROOT=/pipeline --cpus="$cpus" --memory="$mem" \
        -e MVP_ARM=primary_ssclines -e MVP_ADDED="$TAGS" -e MVP_N_EXPECT=600 \
        -e MVP_SUBTITLE=0 "$@" \
        -v "$PWD":/pipeline cline-go:latest \
        Rscript "/pipeline/benchmarks/$s" > "work/j16_logs/${s%.R}.log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "    FAILED rc=$rc (work/j16_logs/${s%.R}.log)"; tail -5 "work/j16_logs/${s%.R}.log"
        exit $rc
    fi
    echo "    ok; warnings: $(grep -ci 'warning' "work/j16_logs/${s%.R}.log")"
}

[[ "$WHAT" == detection || "$WHAT" == all ]] && \
    run mvp_detection_600.R 16 64g -e OPENBLAS_NUM_THREADS=1 -e NCORES=16
[[ "$WHAT" == gallery || "$WHAT" == all ]] && \
    run mvp_j16_gallery.R 8 48g -e OPENBLAS_NUM_THREADS=4

# Cross-arm pooling scan: no emitted table may carry more distinct seeds than the corpus.
echo "--- cross-arm pooling scan (must print nothing):"
for f in benchmarks/mvp_eval/detection600/*.tsv benchmarks/mvp_eval/figures_ssclines_j16_gallery/*.tsv; do
    col=$(head -1 "$f" | tr '\t' '\n' | grep -nx seed | cut -d: -f1)
    [[ -z "$col" ]] && continue
    n=$(tail -n +2 "$f" | cut -f"$col" | sort -u | wc -l)
    (( n > 600 )) && echo "    $f: $n seeds"
done
echo "=== done $(date -Is)"
