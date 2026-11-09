#!/usr/bin/env bash
# run_mode.sh: one pipeline mode through the clinego launcher, with its time and memory recorded.
#
#   manuscript/arabidopsis/run_mode.sh <config> <mode> [snakemake args]
#   manuscript/arabidopsis/run_mode.sh --wait <exit file>
#
# The launcher runs its container with --rm and this kernel (5.15) has no memory.peak, so peak memory
# exists only if it is sampled during the run. Every 15 s this script reads the container's cgroup
# (anon, memory.current) and the RSS and command of each process above 1 GB, so a peak can be tied to
# a script. It also keeps work/ram_guard.sh watching the host for this one container.
#
# The run is detached, so ending the calling session does not end it. The caller blocks on the exit
# file, and returns the run's exit code; if the caller dies, `--wait <exit file>` waits again.
#
# Writes under ANALYSIS/Arabidopsis1135/:
#   runlogs/<config stem>_<mode>_<stamp>.log    the launcher's output
#   runlogs/<config stem>_<mode>_<stamp>.exit   the exit code, written last
#   resources/<config stem>_<mode>_<stamp>.tsv  container memory, one row per sample
#   resources/<config stem>_<mode>_<stamp>_procs.tsv  processes above 1 GB
#   resources/modes.tsv                         one row per finished run
#
# ARAB_CORES and ARAB_MEMORY override --cores 120 and --memory 500g.
#
# The image is the code that runs: clinego executes the pipeline baked at /pipeline, not this
# checkout, and --image skips clinego's check that the two agree. So each config is pinned to an
# image below, built with `clinego build --tag` from the commit its runs need. A pipeline change
# reaches a project only through a new image and a new line here.
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd -P)
SELF=$ROOT/manuscript/arabidopsis/run_mode.sh
[ -x "$SELF" ] || { echo "run_mode.sh: $SELF is not executable" >&2; exit 2; }
cd "$ROOT" || exit 2

CORES=${ARAB_CORES:-120}
MEMORY=${ARAB_MEMORY:-500g}
INTERVAL=15
PROC_MIN_KB=1048576
RUNLOGS=ANALYSIS/Arabidopsis1135/runlogs
RESOURCES=ANALYSIS/Arabidopsis1135/resources
USAGE="usage: run_mode.sh <config> <mode> [snakemake args] | run_mode.sh --wait <exit file>"

wait_for() {
    local exit_file=$1
    until [ -s "$exit_file" ]; do sleep 20; done
    local code
    code=$(cat "$exit_file")
    tail -n 5 "${exit_file%.exit}.log"
    echo "run_mode.sh: exit $code ($exit_file)"
    exit "$code"
}

if [ "${1:-}" = "--wait" ]; then
    wait_for "${2:?$USAGE}"
fi

CONFIG=${1:?$USAGE}
MODE=${2:?$USAGE}
shift 2
[ -f "$CONFIG" ] || { echo "run_mode.sh: no config $CONFIG in $ROOT" >&2; exit 2; }
STEM=$(basename "$CONFIG" .yaml)
NAME=clinego-$STEM-$MODE
case "$STEM" in
    # 6e7efbb: the critical r2 from LDdecay.background_draws draws, on top of ce4c704's relatedness
    # on SNPs passing Filter.maf and snp_miss (plan shiny-seeking-scott). 1700e49 ran processing and
    # prestructure; the two images differ only in ld_decay_background, write_summary.R and the app.
    config_Arabidopsis1135maf01) IMAGE=cline-go:6e7efbb ;;
    *) IMAGE=cline-go:6cdedcc ;;
esac

if [ "${ARAB_DETACHED:-}" != 1 ]; then
    mkdir -p "$RUNLOGS" "$RESOURCES"
    STAMP=$(date +%Y%m%d_%H%M%S)
    ARAB_DETACHED=1 ARAB_STAMP=$STAMP setsid nohup "$SELF" "$CONFIG" "$MODE" "$@" \
        > "$RUNLOGS/${STEM}_${MODE}_${STAMP}.wrapper.log" 2>&1 &
    disown
    echo "run_mode.sh: $NAME detached; exit file $RUNLOGS/${STEM}_${MODE}_${STAMP}.exit"
    wait_for "$RUNLOGS/${STEM}_${MODE}_${STAMP}.exit"
fi

# ------------------------------------------------------------------ the detached run
STAMP=${ARAB_STAMP:?}
BASE=${STEM}_${MODE}_${STAMP}
LOG=$RUNLOGS/$BASE.log
TRACE=$RESOURCES/$BASE.tsv
PROCS=$RESOURCES/${BASE}_procs.tsv
MODES=$RESOURCES/modes.tsv
GUARD_LOG=$RUNLOGS/ram_guard.log

docker_cli() { nix shell nixpkgs#docker-client -c docker "$@"; }

if docker_cli ps --format '{{.Names}}' | grep -qx "$NAME"; then
    echo "run_mode.sh: $NAME is already running" >> "$LOG"
    echo 3 > "$RUNLOGS/$BASE.exit"
    exit 3
fi

printf 'time\tanon_bytes\tcurrent_bytes\thost_used_gb\n' > "$TRACE"
printf 'time\tpid\trss_kb\tcommand\n' > "$PROCS"

NAME_PREFIX=$NAME LOG=$GUARD_LOG bash work/ram_guard.sh > /dev/null 2>&1 &
GUARD=$!

extra=()
[ $# -gt 0 ] && extra=(-- "$@")
start=$(date +%s)
nix shell nixpkgs#docker-client -c ./clinego run "$MODE" --config "$CONFIG" --cores "$CORES" \
    --memory "$MEMORY" --image "$IMAGE" ${extra[@]+"${extra[@]}"} > "$LOG" 2>&1 &
RUN=$!

scope=""
peak_anon=0
peak_current=0
while kill -0 "$RUN" 2>/dev/null; do
    if [ -z "$scope" ]; then
        id=$(docker_cli inspect --format '{{.Id}}' "$NAME" 2>/dev/null) || id=""
        [ -n "$id" ] && scope=/sys/fs/cgroup/system.slice/docker-$id.scope
    fi
    if [ -n "$scope" ] && [ -r "$scope/memory.current" ]; then
        now=$(date +%s)
        current=$(cat "$scope/memory.current" 2>/dev/null) || current=0
        anon=$(awk '$1 == "anon" {print $2}' "$scope/memory.stat" 2>/dev/null) || anon=0
        current=${current:-0}
        anon=${anon:-0}
        used=$(free -g | awk '/^Mem:/ {print $3}')
        printf '%s\t%s\t%s\t%s\n' "$now" "$anon" "$current" "$used" >> "$TRACE"
        [ "$anon" -gt "$peak_anon" ] && peak_anon=$anon
        [ "$current" -gt "$peak_current" ] && peak_current=$current
        while read -r pid; do
            rss=$(awk '/^VmRSS:/ {print $2}' "/proc/$pid/status" 2>/dev/null)
            if [ -n "$rss" ] && [ "$rss" -ge "$PROC_MIN_KB" ]; then
                command=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null | cut -c1-200)
                printf '%s\t%s\t%s\t%s\n' "$now" "$pid" "$rss" "$command" >> "$PROCS"
            fi
        done < "$scope/cgroup.procs" 2>/dev/null
    fi
    sleep "$INTERVAL"
done

wait "$RUN"
status=$?
end=$(date +%s)
kill "$GUARD" 2>/dev/null

[ -f "$MODES" ] || printf 'config\tmode\tstart\tend\twall_seconds\tpeak_anon_gb\tpeak_current_gb\texit\tsnakemake_args\n' > "$MODES"
awk -v a="$peak_anon" -v c="$peak_current" 'BEGIN {printf "%.1f\t%.1f", a / 1073741824, c / 1073741824}' > "$RUNLOGS/$BASE.peaks"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$STEM" "$MODE" "$(date -d "@$start" -Is)" "$(date -d "@$end" -Is)" \
    "$((end - start))" "$(cat "$RUNLOGS/$BASE.peaks")" "$status" "$*" >> "$MODES"
rm -f "$RUNLOGS/$BASE.peaks"
echo "$status" > "$RUNLOGS/$BASE.exit"
