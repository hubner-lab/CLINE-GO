#!/usr/bin/env bash
# mvp_chain_remaining.sh -- run the remaining panel sweeps back to back.
#
# Exists because each stage was being launched by hand and the machine sat idle
# between them. The stages MUST be sequential: they target the same 32 result
# trees, and two Snakemake processes on one tree collide on its lock.
#
#   step 2  rand_best1 / rand_solo1 / rand_union1   (139-223 SNPs)  x 112 gardens
#   step 3  neutral_all / all       (5.6k / 11.2k SNPs median)      x 112 gardens
#
# Step 3 is the expensive one and is RAM-bound, not core-bound: its Gradient
# Forest models were 1.6 GB each on the ORIGINAL replicates, and Group 1's are
# 2.4x larger (median `all` = 24k SNPs, max 26k). Hence fewer concurrent seeds
# and a bigger per-seed memory cap than step 2.
#
# SCOPE. Written when the manifest held exactly one cohort, so it took every seed
# in it. That is now wrong by default: rewriting a finished seed's sweep config
# and re-entering its tree re-triggers cached Gradient Forest models that cost
# ~50 min each to rebuild. Pass the cohort's seeds explicitly.
#
# Usage:  mvp_chain_remaining.sh [SEEDS_CSV|all]
set -uo pipefail
ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
cd "$ROOT"
DOCKER=(nix shell nixpkgs#docker-client -c docker)
GATE="$ROOT/benchmarks/mvp_eval/gate"
SEEDS_ARG="${1:-all}"
if [[ "$SEEDS_ARG" == "all" ]]; then
    SEEDS=$(awk -F'\t' 'NR>1{print $1}' benchmarks/mvp_seeds.tsv)
    echo "scope: ALL manifest seeds ($(echo "$SEEDS" | wc -w))" >&2
else
    SEEDS=${SEEDS_ARG//,/ }
    echo "scope: $(echo "$SEEDS" | wc -w) requested seed(s)"
fi
SEEDS_CSV=$(echo "$SEEDS" | tr ' \n' ',,' | sed 's/,*$//;s/^,*//')

run_stage() {           # run_stage <suffix> <sets> <seed_jobs> <cores> <mem>
    local sfx="$1" sets="$2" jobs="$3" cores="$4" mem="$5"
    echo "=== stage $sfx : $sets ($(date -Is)) ==="
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -v "$ROOT:/pipeline" -w /pipeline \
        -e PIPELINE_ROOT=/pipeline cline-go:latest \
        Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds="$SEEDS_CSV" \
            --gardens_subdir=gardens --suffix="$sfx" --sets="$sets" 2>&1 | tail -1

    # Clear any lock left by an aborted run; a stale lock fails a seed in seconds
    # and the driver would skip it for the whole stage.
    for s in $SEEDS; do
        "${DOCKER[@]}" run --user "$(id -u):$(id -g)" -e USER=cline-go --rm \
            -w /pipeline -v "$ROOT:/pipeline" cline-go:latest \
            snakemake --unlock -s Snakefile --config mode=maladaptation \
            --configfile "config_MVP${s}${sfx}.yaml" >/dev/null 2>&1 &
        while (( $(jobs -rp | wc -l) >= 12 )); do sleep 1; done
    done
    wait

    date -Is > "$GATE/chain${sfx}.start"
    PIPELINE_ROOT="$ROOT" CONFIG_SUFFIX="$sfx" BLAS_THREADS=2 HOST_RAM_CEILING_GB=800 \
        bash benchmarks/mvp_garden_sweep.sh "$SEEDS_CSV" "$jobs" "$cores" "$mem" \
        > "$GATE/chain${sfx}.log" 2>&1
    echo "=== stage $sfx done ($(date -Is)): DONE=$(grep -cE '\] DONE ' "$GATE/chain${sfx}.log") FAILED=$(grep -cE 'FAILED' "$GATE/chain${sfx}.log") ==="
}

run_stage _floor "rand_best1,rand_solo1,rand_union1" 12 6 50g
run_stage _big   "neutral_all,all"                    6 8 140g

echo "=== chain complete $(date -Is) ==="
