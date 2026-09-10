#!/usr/bin/env bash
# mvp_garden_run_par.sh -- garden sweep with garden-level parallelism (journal 09, phase 2).
#
# Same work as mvp_garden_run.sh, but each seed gets N CLONE lanes (see mvp_clone_project.sh), so
# N gardens of one seed run at once instead of one. Rationale, measured: a single garden's DAG is
# only ~4 jobs wide in its tail, so one lane per seed used ~2 of its 8 CPUs and 14 lanes reached
# 27 of 120 cores. With N=4 the sweep runs 56 lanes at ~2 cores each.
#
# Gardens already harvested (offset09/gardens/{seed}/{garden}/.complete) are skipped, so this
# resumes whatever the serial driver already produced.
#
# Usage:  mvp_garden_run_par.sh [SEEDS_CSV|all] [CLONES_PER_SEED] [CPUS_PER_LANE] [SNAKE_CORES]
set -uo pipefail

ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
SEEDS_ARG="${1:-all}"
CLONES="${2:-4}"
CLONE_ONLY="${CLONE_ONLY:-4}"
CPUS_PER_LANE="${3:-2}"
SNAKE_CORES="${4:-4}"
MEM_PER_LANE="${MEM_PER_LANE:-20g}"
BLAS_THREADS="${BLAS_THREADS:-2}"
HOST_RAM_CEILING_GB="${HOST_RAM_CEILING_GB:-900}"

MANIFEST="$ROOT/benchmarks/mvp_seeds.tsv"
OUT="${OUT:-$ROOT/benchmarks/mvp_eval/offset09}"
LOGDIR="${LOGDIR:-$OUT/runlogs}"
mkdir -p "$LOGDIR"
DOCKER=(nix shell nixpkgs#docker-client -c docker)

if [[ "$SEEDS_ARG" == "all" ]]; then
    mapfile -t SEEDS < <(awk 'NR>1 && NF {print $1}' "$MANIFEST")
else
    IFS=',' read -r -a SEEDS <<< "$SEEDS_ARG"
fi
echo "seeds=${#SEEDS[@]} clones=$CLONES cpus/lane=$CPUS_PER_LANE snake -c$SNAKE_CORES mem=$MEM_PER_LANE"
echo "total lanes = $(( ${#SEEDS[@]} * CLONES ))"

watchdog() {
    while true; do
        u=$(free -g | awk '/^Mem:/ {print $3}')
        if [[ -n "$u" ]] && (( u > HOST_RAM_CEILING_GB )); then
            v=$("${DOCKER[@]}" ps --filter "name=gp09_" --format '{{.Names}}' | head -1)
            echo "[watchdog] host RAM ${u}G > ${HOST_RAM_CEILING_GB}G -- stopping ${v:-<none>}" \
                | tee -a "$LOGDIR/watchdog.log"
            [[ -n "$v" ]] && "${DOCKER[@]}" stop "$v" >/dev/null 2>&1
        fi
        sleep 30
    done
}
watchdog & WD=$!
trap 'kill "$WD" 2>/dev/null' EXIT

snake() {   # snake <name> <configfile> <log> <cores>
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm --name "$1" \
        -e USER=cline-go -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" \
        -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --cpus="$CPUS_PER_LANE" --memory="$MEM_PER_LANE" \
        -v "$ROOT:/pipeline" "$IMAGE" \
        snakemake -c"$4" -s Snakefile -R stage_custom_climate_future \
            --config mode=maladaptation --configfile "$2" \
            --scheduler greedy --rerun-incomplete >> "$3" 2>&1
}

run_lane() {                     # run_lane <seed> <clone_index>
    local seed="$1" k="$2"
    local tag="p${k}"
    local proj="MVP${seed}${tag}"
    local res="$ROOT/${proj}_results"
    local gardens_f="$OUT/gardens_${seed}.tsv"
    local presfile="$ROOT/data/mvp/MVP${seed}/MVP${seed}_env_present.tsv"
    local envfile="$ROOT/data/mvp/MVP${seed}/${proj}_env_garden.tsv"
    local log="$LOGDIR/${proj}.log"
    : > "$log"

    "$ROOT/benchmarks/mvp_clone_project.sh" "$seed" "$tag" >> "$log" 2>&1 || {
        echo "[$proj] FATAL: clone failed" | tee -a "$log"; return 1; }

    # every CLONES-th garden, offset by this clone's index
    mapfile -t GARDENS < <(awk -F'\t' -v n="$CLONES" -v k="$k" \
        'NR>1 {i++; if ((i-1) % n == (k-1)) print $1}' "$gardens_f")
    echo "[$proj] ${#GARDENS[@]} gardens $(date -Is)" | tee -a "$log"

    for garden in "${GARDENS[@]}"; do
        local dest="$OUT/gardens/${seed}/${garden}"
        [[ -f "$dest/.complete" ]] && continue
        mkdir -p "$dest"

        read -r OPT1 OPT2 < <(awk -F'\t' -v g="$garden" 'NR>1 && $1==g {print $4"\t"$5; exit}' "$gardens_f")
        [[ -z "${OPT1:-}" ]] && { echo "[$proj] FATAL: garden $garden not found" | tee -a "$log"; return 1; }
        awk -F'\t' -v o1="$OPT1" -v o2="$OPT2" 'BEGIN{OFS="\t"}
             NR==1 {print "site","bio_1","bio_2"; next} {print $1, o1, o2}' "$presfile" > "$envfile"

        local t0 rc; t0=$(date +%s)
        snake "gp09_${proj}" "config_${proj}_g.yaml" "$log" "$SNAKE_CORES"; rc=$?
        (( rc != 0 )) && { snake "gp09_${proj}" "config_${proj}_g.yaml" "$log" 2; rc=$?; }
        (( rc != 0 )) && { echo "[$proj] FAILED $garden pass A" | tee -a "$log"; return 1; }
        for m in gradient_forest geometric_offset rda_offset; do
            for d in "$res/Maladaptation/tables/$m"/*_nospatial; do
                [[ -d "$d" ]] || continue
                cp "$d/genetic_offset_site.tsv" "$dest/${m}__$(basename "$d" _nospatial).tsv" 2>/dev/null
            done
        done

        snake "gp09c_${proj}" "config_${proj}_grdac.yaml" "$log" "$SNAKE_CORES"; rc=$?
        (( rc != 0 )) && { snake "gp09c_${proj}" "config_${proj}_grdac.yaml" "$log" 2; rc=$?; }
        (( rc != 0 )) && { echo "[$proj] FAILED $garden pass B" | tee -a "$log"; return 1; }
        for d in "$res/Maladaptation/tables/rda_offset"/*_nospatial; do
            [[ -d "$d" ]] || continue
            cp "$d/genetic_offset_site.tsv" "$dest/rda_corrected__$(basename "$d" _nospatial).tsv" 2>/dev/null
        done

        touch "$dest/.complete"
        echo "[$proj] $garden $(($(date +%s) - t0))s files=$(ls "$dest"/*.tsv 2>/dev/null | wc -l)" | tee -a "$log"
    done
    echo "[$proj] DONE $(date -Is)" | tee -a "$log"
}

for seed in "${SEEDS[@]}"; do
    run_lane "$seed" "$CLONE_ONLY" &
done
wait

echo; echo "=== summary ==="
tot=0
for seed in "${SEEDS[@]}"; do
    n=$(ls -d "$OUT/gardens/${seed}"/*/.complete 2>/dev/null | wc -l)
    tot=$((tot + n)); echo "MVP${seed}  gardens complete: $n"
done
echo "total gardens: $tot"
