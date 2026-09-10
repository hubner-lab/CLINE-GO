#!/usr/bin/env bash
# mvp_offset_run.sh -- run mode=maladaptation for the offset arm, one project per seed.
#
# Each seed is independent (its own {PROJECT}_results/ tree and its own Snakemake lock), so
# seeds run CONCURRENTLY. Modes are never run concurrently WITHIN a seed -- there is only
# one mode here, so that constraint is satisfied by construction.
#
# Docker client: /snap/bin/docker cannot start on this host (setuid blocked under
# nosuid+NoNewPrivs), so every call goes through `nix shell nixpkgs#docker-client -c docker`
# -- see CLAUDE.local.md. `-e USER=cline-go` is mandatory: under --user $(id -u):$(id -g)
# there is no /etc/passwd entry inside the container and snakemake's getpass.getuser() dies
# with KeyError: getpwuid(). OPENBLAS/OMP thread caps are mandatory too: containers see all
# 192 host cores and uncapped BLAS runs far slower while eating tens of GB.
#
# Resource budget: 120 CPUs / 600 GB across all our containers, and a host-level guardrail
# at 900 GB total server RAM (this is a shared machine). The watchdog below samples host
# usage and, if it crosses the ceiling, stops OUR largest container first -- never anyone
# else's.
#
# Usage:  benchmarks/mvp_offset_run.sh [SEEDS_CSV|all] [PARALLEL] [CPUS_PER_SEED]
set -uo pipefail

ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
SEEDS_ARG="${1:-all}"
PARALLEL="${2:-7}"
CPUS_PER_SEED="${3:-16}"
MEM_PER_SEED="${MEM_PER_SEED:-80g}"
SNAKE_CORES="${SNAKE_CORES:-8}"
BLAS_THREADS="${BLAS_THREADS:-8}"
CFG_SUFFIX="${CFG_SUFFIX:-_off}"
HOST_RAM_CEILING_GB="${HOST_RAM_CEILING_GB:-900}"
# Rules to force-rerun. The future-climate TABLE reaches stage_custom_climate_future as a
# PARAM, not an input, so Snakemake cannot see that the scenario changed -- editing the
# scenario and re-running without this flag reports "nothing to be done". Forcing that one
# rule re-fires every offset downstream of it while correctly leaving the GF adaptive models
# alone: they are built from the PRESENT climate and the SNP set, neither of which moved.
RERUN_RULE="${RERUN_RULE:-stage_custom_climate_future}"

MANIFEST="$ROOT/benchmarks/mvp_seeds.tsv"
LOGDIR="${LOGDIR:-$ROOT/benchmarks/mvp_eval/offset08/runlogs}"
mkdir -p "$LOGDIR"

DOCKER=(nix shell nixpkgs#docker-client -c docker)

if [[ "$SEEDS_ARG" == "all" ]]; then
    mapfile -t SEEDS < <(awk 'NR>1 && NF {print $1}' "$MANIFEST")
else
    IFS=',' read -r -a SEEDS <<< "$SEEDS_ARG"
fi
echo "Seeds (${#SEEDS[@]}): ${SEEDS[*]}"
echo "parallel=$PARALLEL  cpus/seed=$CPUS_PER_SEED  mem/seed=$MEM_PER_SEED  snake -c$SNAKE_CORES"

# ---- host RAM watchdog ------------------------------------------------------
# Samples total host usage every 30 s. Above the ceiling it stops the largest-memory
# container whose name starts with off08_ (ours), one per sample, and logs what it did.
watchdog() {
    while true; do
        used_gb=$(free -g | awk '/^Mem:/ {print $3}')
        if [[ -n "$used_gb" ]] && (( used_gb > HOST_RAM_CEILING_GB )); then
            victim=$("${DOCKER[@]}" ps --filter "name=off08_" --format '{{.Names}}' | head -1)
            echo "[watchdog] host RAM ${used_gb}G > ${HOST_RAM_CEILING_GB}G -- stopping ${victim:-<none>}" \
                | tee -a "$LOGDIR/watchdog.log"
            [[ -n "$victim" ]] && "${DOCKER[@]}" stop "$victim" >/dev/null 2>&1
        fi
        sleep 30
    done
}
watchdog & WATCHDOG_PID=$!
trap 'kill "$WATCHDOG_PID" 2>/dev/null' EXIT

run_seed() {
    local seed="$1"
    local proj="MVP${seed}"
    local cfg="config_${proj}${CFG_SUFFIX}.yaml"
    local log="$LOGDIR/${proj}.log"
    : > "$log"

    if [[ ! -f "$ROOT/$cfg" ]]; then
        echo "[$proj] FATAL: missing $cfg" | tee -a "$log"; return 1
    fi

    # A lock present at seed start is a corpse: this driver is the only thing that runs
    # snakemake against MVP* projects and it never starts two runs for one seed.
    local res="$ROOT/${proj}_results"
    if [[ -d "$res/.snakemake/locks" ]] && [[ -n "$(ls -A "$res/.snakemake/locks" 2>/dev/null)" ]]; then
        echo "[$proj] clearing stale snakemake lock" | tee -a "$log"
        "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=cline-go \
            -v "$ROOT:/pipeline" "$IMAGE" \
            snakemake -s Snakefile --unlock --config mode=maladaptation \
                --configfile "$cfg" >> "$log" 2>&1
    fi

    echo "[$proj] start $(date -Is)" | tee -a "$log"
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm --name "off08_${proj}" \
        -e USER=cline-go -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" \
        -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --cpus="$CPUS_PER_SEED" --memory="$MEM_PER_SEED" \
        -v "$ROOT:/pipeline" "$IMAGE" \
        snakemake -c"$SNAKE_CORES" -s Snakefile \
            ${RERUN_RULE:+-R "$RERUN_RULE"} \
            --config mode=maladaptation --configfile "$cfg" \
            --scheduler greedy --rerun-incomplete \
        >> "$log" 2>&1
    local rc=$?
    echo "[$proj] exit=$rc $(date -Is)" | tee -a "$log"
    if (( rc == 0 )); then echo "[$proj] DONE" >> "$log"; fi
    return $rc
}

fail=0
running=0
for seed in "${SEEDS[@]}"; do
    run_seed "$seed" &
    running=$((running + 1))
    if (( running >= PARALLEL )); then wait -n || fail=$((fail + 1)); running=$((running - 1)); fi
done
while (( running > 0 )); do wait -n || fail=$((fail + 1)); running=$((running - 1)); done

echo
echo "=== summary ==="
for seed in "${SEEDS[@]}"; do
    if grep -q "^\[MVP${seed}\] DONE" "$LOGDIR/MVP${seed}.log" 2>/dev/null; then
        n=$(ls -d "$ROOT/MVP${seed}_results/Maladaptation/tables/"*/*_nospatial 2>/dev/null | wc -l)
        echo "MVP${seed}  OK   ($n method x set outputs)"
    else
        echo "MVP${seed}  FAILED -- see $LOGDIR/MVP${seed}.log"
    fi
done
echo "failed seeds: $fail"
exit $(( fail > 0 ))
