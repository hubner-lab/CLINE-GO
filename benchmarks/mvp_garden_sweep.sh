#!/usr/bin/env bash
# mvp_garden_sweep.sh -- clone-free garden sweep. One Snakemake invocation per seed.
#
# Replaces mvp_garden_run_par.sh, which ran 2 invocations per garden (84 per seed)
# across 4 clone lanes, rewriting a single env table between each. That shape existed
# only because (a) future_table was a scalar param and (b) RDA-corrected was not a
# registered method. Both are fixed, so a seed is now ONE invocation covering every
# garden and all four offset methods.
#
# What this means operationally:
#   - no clone lanes            -> no mvp_clone_project.sh, and no MVP{seed}p{1..4}_results
#   - no per-garden harvest     -> outputs land at their final scenario paths directly
#   - no .complete marker dance -> Snakemake's own DAG is the resume mechanism
#   - no mvp_unstick.sh         -> the bare-`wait` hang is gone (lane PIDs are waited on)
#
# The HARVESTED LAYOUT IS UNCHANGED. eval_offset_lind.R reads exclusively from
# offset09/gardens/{seed}/{garden}/{method}__{panel}.tsv, so this copies into that
# exact shape once per seed after the run, instead of 42 times during it.
#
# PREREQUISITES per seed, in order:
#   Rscript benchmarks/mvp_write_garden_env.R --seed=SEED     # 42 scenario env tables
#   Rscript benchmarks/mvp_write_sweep_config.R --seeds=SEED  # config_MVP{seed}_sweep.yaml
#   models for all 7 panels must exist in the base tree, else the run rebuilds them:
#     ls MVP{seed}_results/_intermediate/gradient_forest/*/adaptive_model.qs
#
# Usage:  mvp_garden_sweep.sh [SEEDS_CSV|all] [SEED_JOBS] [SNAKE_CORES] [MEM_PER_SEED]
set -uo pipefail

ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
SEEDS_ARG="${1:-all}"
SEED_JOBS="${2:-6}"          # seeds running concurrently
SNAKE_CORES="${3:-28}"       # DAG is 28 offset jobs wide per seed
MEM_PER_SEED="${4:-80g}"
BLAS_THREADS="${BLAS_THREADS:-4}"
HOST_RAM_CEILING_GB="${HOST_RAM_CEILING_GB:-900}"

MANIFEST="$ROOT/benchmarks/mvp_seeds.tsv"
OUT="${OUT:-$ROOT/benchmarks/mvp_eval/offset09}"
LOGDIR="${LOGDIR:-$OUT/sweeplogs}"
mkdir -p "$LOGDIR"
DOCKER=(nix shell nixpkgs#docker-client -c docker)

if [[ "$SEEDS_ARG" == "all" ]]; then
    mapfile -t SEEDS < <(awk 'NR>1 && NF {print $1}' "$MANIFEST")
else
    IFS=',' read -r -a SEEDS <<< "$SEEDS_ARG"
fi
echo "seeds=${#SEEDS[@]} concurrent=$SEED_JOBS snake -c$SNAKE_CORES mem/seed=$MEM_PER_SEED"

# Host-RAM guardrail. Shared machine: roughly a third of it belongs to other users,
# so the ceiling is on TOTAL host usage, not on this sweep's share.
watchdog() {
    while true; do
        u=$(free -g | awk '/^Mem:/ {print $3}')
        if [[ -n "$u" ]] && (( u > HOST_RAM_CEILING_GB )); then
            v=$("${DOCKER[@]}" ps --filter "name=gsweep_" --format '{{.Names}}' | head -1)
            echo "[watchdog] host RAM ${u}G > ${HOST_RAM_CEILING_GB}G -- stopping ${v:-<none>}" \
                | tee -a "$LOGDIR/watchdog.log"
            [[ -n "$v" ]] && "${DOCKER[@]}" stop "$v" >/dev/null 2>&1
        fi
        sleep 30
    done
}
watchdog & WD=$!
# Kill the watchdog explicitly on exit. The old driver ended with a bare `wait`, which
# also waited on this infinite loop and hung the driver forever after its last lane.
trap 'kill "$WD" 2>/dev/null' EXIT

run_seed() {
    local seed="$1"
    local proj="MVP${seed}"
    # CONFIG_SUFFIX selects the garden phase: _p1 (the 42-garden journal-09 set),
    # _p2 (the remaining 70), or _sweep (all 112). Phases write to disjoint
    # per-scenario paths and merge in the harvest, so running p1 then p2 costs the
    # same total work as running all 112 at once -- but yields results after p1.
    local cfg="config_${proj}${CONFIG_SUFFIX:-_sweep}.yaml"
    local res="$ROOT/${proj}_results"
    local log="$LOGDIR/${proj}.log"
    : > "$log"

    if [[ ! -f "$ROOT/$cfg" ]]; then
        echo "[$proj] FATAL: missing $cfg -- run mvp_write_sweep_config.R first" | tee -a "$log"
        return 1
    fi

    echo "[$proj] start $(date -Is)" | tee -a "$log"
    local t0; t0=$(date +%s)

    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm --name "gsweep_${proj}" \
        -e USER=cline-go -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" \
        -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --memory="$MEM_PER_SEED" \
        -w /pipeline -v "$ROOT:/pipeline" "$IMAGE" \
        snakemake -c"$SNAKE_CORES" -s Snakefile \
            --config mode=maladaptation --configfile "$cfg" \
            --scheduler greedy --rerun-incomplete >> "$log" 2>&1
    local rc=$?
    if (( rc != 0 )); then
        echo "[$proj] FAILED rc=$rc after $(($(date +%s) - t0))s" | tee -a "$log"
        return 1
    fi

    # Harvest into the layout eval_offset_lind.R reads. Once per seed, not per garden.
    local n=0
    local harvested="$LOGDIR/.harvested_${seed}"
    : > "$harvested"
    for mdir in "$res/Maladaptation/tables"/*; do
        [[ -d "$mdir" ]] || continue
        local method; method=$(basename "$mdir")
        # The harvested vocabulary predates the method registry: eval_offset_lind.R's
        # METHOD_LABEL keys the corrected RDA as `rda_corrected`. The pipeline now calls
        # that method `rda_offset_corrected`. Map on the way out so the harvested layout
        # stays exactly what the scorer already reads -- renaming it there instead would
        # invalidate the 11 seeds harvested by the old driver.
        local token="$method"
        [[ "$token" == "rda_offset_corrected" ]] && token="rda_corrected"
        for pdir in "$mdir"/*_nospatial; do
            [[ -d "$pdir" ]] || continue
            local panel; panel=$(basename "$pdir" _nospatial)
            for sdir in "$pdir"/*/; do
                [[ -f "$sdir/genetic_offset_site.tsv" ]] || continue
                local garden; garden=$(basename "$sdir")
                mkdir -p "$OUT/gardens/${seed}/${garden}"
                cp "$sdir/genetic_offset_site.tsv" \
                   "$OUT/gardens/${seed}/${garden}/${token}__${panel}.tsv"
                echo "$garden" >> "$harvested"
                n=$((n + 1))
            done
        done
    done
    # Mark complete ONLY the gardens this run actually harvested. Globbing every
    # directory under gardens/{seed} would also stamp the empty leftovers that the
    # old per-garden driver mkdir'd for gardens it never finished, inflating the
    # completeness count (they hold no TSVs, so scoring ignores them either way).
    sort -u "$harvested" | while read -r garden; do
        [[ -n "$garden" ]] && touch "$OUT/gardens/${seed}/${garden}/.complete"
    done
    rm -f "$harvested"

    echo "[$proj] DONE $(($(date +%s) - t0))s harvested=$n $(date -Is)" | tee -a "$log"
}

PIDS=()
# Throttle on the seed PIDs only. `jobs -rp` would also count the RAM watchdog,
# which never exits -- that silently costs one slot, so N seeds would run N-1 wide.
alive() { local n=0 p; for p in "${PIDS[@]}"; do kill -0 "$p" 2>/dev/null && n=$((n+1)); done; echo "$n"; }
for seed in "${SEEDS[@]}"; do
    while (( $(alive) >= SEED_JOBS )); do sleep 5; done
    run_seed "$seed" & PIDS+=($!)
done
# Wait only on the seed jobs -- never on $WD, which never exits.
for p in "${PIDS[@]}"; do wait "$p"; done

echo; echo "=== summary ==="
tot=0
for seed in "${SEEDS[@]}"; do
    n=$(ls -d "$OUT/gardens/${seed}"/*/.complete 2>/dev/null | wc -l)
    tot=$((tot + n)); echo "MVP${seed}  gardens complete: $n"
done
echo "total gardens: $tot"
