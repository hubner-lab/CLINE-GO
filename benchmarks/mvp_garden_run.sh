#!/usr/bin/env bash
# mvp_garden_run.sh -- run the Lind & Lotterhos garden sweep (journal 09).
#
# For every seed, for every garden: write that garden's environment as the future-climate table,
# run mode=maladaptation, and HARVEST the per-site offsets before the next garden overwrites them.
# Two passes per garden, because the pipeline registry has one rda_offset entry and the two RDA
# implementations differ only by condition_pcs:
#     pass A (config_MVP{seed}_g.yaml)      GFoffset + LFMM2offset + RDA-uncorrected
#     pass B (config_MVP{seed}_grdac.yaml)  RDA-corrected only
#
# WHY -R stage_custom_climate_future: the future-climate TABLE reaches that rule as a PARAM, not an
# input, so Snakemake cannot see that the garden changed and would report "nothing to be done".
# Forcing that one rule re-fires every offset below it while correctly leaving the Gradient Forest
# models alone -- they are built from the PRESENT climate and the SNP set, neither of which moves
# between gardens, so each is fitted once per (seed, set) and reused across all gardens.
#
# Parallelism: seeds are independent (own results tree, own Snakemake lock) and run concurrently.
# Gardens are sequential WITHIN a seed -- they share one project directory, so they cannot overlap.
#
# RAM: per-container --memory cap, a host watchdog that stops OUR largest container above the host
# ceiling, and a per-garden OOMKilled check with one retry at reduced cores.
#
# Usage:  benchmarks/mvp_garden_run.sh [SEEDS_CSV|all] [PARALLEL] [CPUS_PER_SEED] [GARDENS]
#   GARDENS: 'all' | an integer N (first N gardens) | a comma list of garden_ids
set -uo pipefail

ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
SEEDS_ARG="${1:-all}"
PARALLEL="${2:-14}"
CPUS_PER_SEED="${3:-8}"
GARDENS_ARG="${4:-all}"
MEM_PER_SEED="${MEM_PER_SEED:-40g}"
SNAKE_CORES="${SNAKE_CORES:-6}"
BLAS_THREADS="${BLAS_THREADS:-6}"
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

echo "Seeds (${#SEEDS[@]}): ${SEEDS[*]}"
echo "gardens=$GARDENS_ARG  parallel=$PARALLEL  cpus/seed=$CPUS_PER_SEED  mem/seed=$MEM_PER_SEED"

# ---- host RAM watchdog ------------------------------------------------------
watchdog() {
    while true; do
        used_gb=$(free -g | awk '/^Mem:/ {print $3}')
        if [[ -n "$used_gb" ]] && (( used_gb > HOST_RAM_CEILING_GB )); then
            victim=$("${DOCKER[@]}" ps --filter "name=g09_" --format '{{.Names}}' | head -1)
            echo "[watchdog] host RAM ${used_gb}G > ${HOST_RAM_CEILING_GB}G -- stopping ${victim:-<none>}" \
                | tee -a "$LOGDIR/watchdog.log"
            [[ -n "$victim" ]] && "${DOCKER[@]}" stop "$victim" >/dev/null 2>&1
        fi
        sleep 30
    done
}
watchdog & WATCHDOG_PID=$!
trap 'kill "$WATCHDOG_PID" 2>/dev/null' EXIT

# ---- one snakemake invocation ----------------------------------------------
snake() {   # snake <container_name> <configfile> <logfile> <cores>
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm --name "$1" \
        -e USER=cline-go -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" \
        -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --cpus="$CPUS_PER_SEED" --memory="$MEM_PER_SEED" \
        -v "$ROOT:/pipeline" "$IMAGE" \
        snakemake -c"$4" -s Snakefile -R stage_custom_climate_future \
            --config mode=maladaptation --configfile "$2" \
            --scheduler greedy --rerun-incomplete \
        >> "$3" 2>&1
}

run_seed() {
    local seed="$1"
    local proj="MVP${seed}"
    local res="$ROOT/${proj}_results"
    local gardens_f="$OUT/gardens_${seed}.tsv"
    local envfile="$ROOT/data/mvp/${proj}/${proj}_env_garden.tsv"
    local presfile="$ROOT/data/mvp/${proj}/${proj}_env_present.tsv"
    local log="$LOGDIR/${proj}.log"
    : > "$log"

    [[ -f "$gardens_f" ]] || { echo "[$proj] FATAL: missing $gardens_f" | tee -a "$log"; return 1; }
    for cfg in "config_${proj}_g.yaml" "config_${proj}_grdac.yaml"; do
        [[ -f "$ROOT/$cfg" ]] || { echo "[$proj] FATAL: missing $cfg" | tee -a "$log"; return 1; }
    done

    # garden list
    local -a GARDENS
    if [[ "$GARDENS_ARG" == "all" ]]; then
        mapfile -t GARDENS < <(awk -F'\t' 'NR>1 {print $1}' "$gardens_f")
    elif [[ "$GARDENS_ARG" =~ ^[0-9]+$ ]]; then
        mapfile -t GARDENS < <(awk -F'\t' -v n="$GARDENS_ARG" 'NR>1 && NR<=n+1 {print $1}' "$gardens_f")
    else
        IFS=',' read -r -a GARDENS <<< "$GARDENS_ARG"
    fi
    echo "[$proj] ${#GARDENS[@]} gardens  start $(date -Is)" | tee -a "$log"

    # a stale lock at seed start is a corpse: this driver is the only writer of MVP* projects
    if [[ -d "$res/.snakemake/locks" ]] && [[ -n "$(ls -A "$res/.snakemake/locks" 2>/dev/null)" ]]; then
        echo "[$proj] clearing stale snakemake lock" | tee -a "$log"
        "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=cline-go \
            -v "$ROOT:/pipeline" "$IMAGE" \
            snakemake -s Snakefile --unlock --config mode=maladaptation \
                --configfile "config_${proj}_g.yaml" >> "$log" 2>&1
    fi

    local g_i=0
    for garden in "${GARDENS[@]}"; do
        g_i=$((g_i + 1))
        local dest="$OUT/gardens/${seed}/${garden}"
        if [[ -f "$dest/.complete" ]]; then continue; fi
        mkdir -p "$dest"

        # ---- write this garden's future-climate table: EVERY deme's future = the garden's env
        awk -F'\t' -v g="$garden" 'NR>1 && $1==g {print $4"\t"$5; exit}' "$gardens_f" \
            > "$dest/.opt" || true
        read -r OPT1 OPT2 < "$dest/.opt"
        if [[ -z "${OPT1:-}" ]]; then
            echo "[$proj] FATAL: garden $garden not found in $gardens_f" | tee -a "$log"; return 1
        fi
        awk -F'\t' -v o1="$OPT1" -v o2="$OPT2" 'BEGIN{OFS="\t"}
             NR==1 {print "site","bio_1","bio_2"; next} {print $1, o1, o2}' "$presfile" > "$envfile"

        local t0 rc
        t0=$(date +%s)

        # ---- pass A: GF + geometric + RDA-uncorrected
        snake "g09_${proj}" "config_${proj}_g.yaml" "$log" "$SNAKE_CORES"; rc=$?
        if (( rc != 0 )); then
            echo "[$proj] garden $garden pass A rc=$rc -- retry at -c4" | tee -a "$log"
            snake "g09_${proj}" "config_${proj}_g.yaml" "$log" 4; rc=$?
        fi
        if (( rc != 0 )); then
            echo "[$proj] FAILED garden $garden pass A" | tee -a "$log"; return 1
        fi
        for m in gradient_forest geometric_offset rda_offset; do
            for d in "$res/Maladaptation/tables/$m"/*_nospatial; do
                [[ -d "$d" ]] || continue
                local set_name; set_name=$(basename "$d" _nospatial)
                cp "$d/genetic_offset_site.tsv" "$dest/${m}__${set_name}.tsv" 2>/dev/null
            done
        done

        # ---- pass B: RDA-corrected
        snake "g09c_${proj}" "config_${proj}_grdac.yaml" "$log" "$SNAKE_CORES"; rc=$?
        if (( rc != 0 )); then
            echo "[$proj] garden $garden pass B rc=$rc -- retry at -c4" | tee -a "$log"
            snake "g09c_${proj}" "config_${proj}_grdac.yaml" "$log" 4; rc=$?
        fi
        if (( rc != 0 )); then
            echo "[$proj] FAILED garden $garden pass B" | tee -a "$log"; return 1
        fi
        for d in "$res/Maladaptation/tables/rda_offset"/*_nospatial; do
            [[ -d "$d" ]] || continue
            local set_name; set_name=$(basename "$d" _nospatial)
            cp "$d/genetic_offset_site.tsv" "$dest/rda_corrected__${set_name}.tsv" 2>/dev/null
        done

        rm -f "$dest/.opt"
        local n_files; n_files=$(ls "$dest"/*.tsv 2>/dev/null | wc -l)
        touch "$dest/.complete"
        echo "[$proj] garden $g_i/${#GARDENS[@]} $garden  $(($(date +%s) - t0))s  files=$n_files" \
            | tee -a "$log"
    done

    echo "[$proj] DONE $(date -Is)" | tee -a "$log"
    return 0
}

fail=0; running=0
for seed in "${SEEDS[@]}"; do
    run_seed "$seed" &
    running=$((running + 1))
    if (( running >= PARALLEL )); then wait -n || fail=$((fail + 1)); running=$((running - 1)); fi
done
while (( running > 0 )); do wait -n || fail=$((fail + 1)); running=$((running - 1)); done

echo; echo "=== summary ==="
for seed in "${SEEDS[@]}"; do
    n=$(ls -d "$OUT/gardens/${seed}"/*/ 2>/dev/null | wc -l)
    f=$(ls "$OUT/gardens/${seed}"/*/*.tsv 2>/dev/null | wc -l)
    if grep -q "^\[MVP${seed}\] DONE" "$LOGDIR/MVP${seed}.log" 2>/dev/null; then
        echo "MVP${seed}  OK      gardens=$n  offset files=$f"
    else
        echo "MVP${seed}  FAILED  gardens=$n  offset files=$f  -- see $LOGDIR/MVP${seed}.log"
    fi
done
echo "failed lanes: $fail"
exit $(( fail > 0 ))
