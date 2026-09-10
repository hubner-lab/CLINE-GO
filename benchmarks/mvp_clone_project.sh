#!/usr/bin/env bash
# mvp_clone_project.sh -- make a lightweight clone of a project results tree so that several
# gardens of the SAME seed can run concurrently.
#
# WHY THIS IS NEEDED. All gardens of a seed share one {PROJECT}_results/ tree and therefore one
# Snakemake lock, so they can only run one at a time. Measured consequence: a lane uses ~2 of its
# 8 CPUs, because a garden's DAG is only ~4 jobs wide in its tail (one rda_offset per marker
# panel), and 14 lanes reach just ~27 of 120 cores.
#
# WHY IT IS SAFE. Every work/intermediate filename derives from the VCF BASENAME
# (common.smk: VCF_BASE = get_vcf_basename(VCF_RAW)); only OUTDIR/LOGDIR use project_name. So a
# clone with a different project_name finds byte-identical inputs at identical paths.
#
# WHAT IS CLONED. A WHITELIST of read-only inputs, hardlinked (instant, ~0 disk). Nothing the
# maladaptation DAG writes is cloned, so every output lands on a fresh inode and can never
# overwrite the original through a shared hardlink:
#   _work/                      filtered/imputed matrices, PCA projections, vcfsnp, removed
#   _intermediate/samples/      sample lists and metadata
#   _intermediate/climate_subset/  row-subset LFMM matrices
#   _intermediate/snp_sets/     the marker panels
#   _intermediate/gradient_forest/**/adaptive_model.qs   fitted GF models (read-only here)
#   climate/tables/present/, climate/rasters/present/    staged present climate
#   PreStructure/tables/        Q-matrices used by the piemap rules
#
# Usage:  mvp_clone_project.sh <SEED> <CLONE_TAG>      e.g. 1232548 c1
#         -> creates MVP{SEED}{CLONE_TAG}_results/ and MVP{SEED}{CLONE_TAG}_logs/
set -euo pipefail

ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
SEED="$1"; TAG="$2"
SRC="$ROOT/MVP${SEED}_results"
DST="$ROOT/MVP${SEED}${TAG}_results"

[[ -d "$SRC" ]] || { echo "FATAL: no source tree $SRC" >&2; exit 1; }

rm -rf "$DST"
mkdir -p "$DST" "$ROOT/MVP${SEED}${TAG}_logs"

link_dir() {   # link_dir <relative path>
    local rel="$1"
    [[ -d "$SRC/$rel" ]] || return 0
    mkdir -p "$(dirname "$DST/$rel")"
    cp -al "$SRC/$rel" "$DST/$rel"
}

link_dir "_work"
link_dir "_intermediate/samples"
link_dir "_intermediate/climate_subset"
link_dir "_intermediate/snp_sets"
link_dir "_intermediate/annotation"
link_dir "climate/tables/present"
link_dir "climate/rasters/present"
link_dir "PreStructure/tables"

# Gradient Forest models: reused read-only across gardens, so the clone must have them or every
# garden would refit a forest (3105 s on the `all` panel).
if [[ -d "$SRC/_intermediate/gradient_forest" ]]; then
    mkdir -p "$DST/_intermediate/gradient_forest"
    for d in "$SRC/_intermediate/gradient_forest"/*/; do
        [[ -d "$d" ]] || continue
        b=$(basename "$d")
        mkdir -p "$DST/_intermediate/gradient_forest/$b"
        # models only -- NOT offset_raster.tif, which this DAG rewrites per garden
        for f in "$d"adaptive_model.qs "$d"random_model.qs; do
            [[ -f "$f" ]] && cp -al "$f" "$DST/_intermediate/gradient_forest/$b/"
        done
    done
fi

echo "cloned MVP${SEED} -> MVP${SEED}${TAG} ($(du -sh --apparent-size "$DST" 2>/dev/null | cut -f1) apparent, hardlinked)"
