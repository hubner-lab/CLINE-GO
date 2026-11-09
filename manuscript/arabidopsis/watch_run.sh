#!/usr/bin/env bash
# watch_run.sh: stall watchdog for one detached run_mode.sh run.
#
#   manuscript/arabidopsis/watch_run.sh <exit file> <mode> [quiet minutes, default 45]
#
# Start it with the Bash tool's run_in_background beside the run_mode.sh call, so that a stall
# wakes the session instead of the run's end notification arriving hours late.
#
# Ends 0 when the run's exit file appears (and prints the run's exit code).
# Ends 2 on a stall: the run log unchanged for the quiet limit while the run's container gained
#   under 60 CPU seconds in that time. A stage that works without writing (sNMF, RDA forward
#   selection) is not a stall; its window restarts and one line says so.
# Ends 3 when the container is gone for 3 checks and no exit file follows.
#
# The log's age comes from `date -r`: `stat` on this host's PATH is LDhat's binary. The container
# id comes from `docker ps --no-trunc`, never `docker inspect`, whose full output holds the
# container's environment. WATCH_INTERVAL (seconds, default 60) sets the check period.
set -uo pipefail

EXIT_FILE=${1:?usage: watch_run.sh <exit file> <mode> [quiet minutes]}
MODE=${2:?usage: watch_run.sh <exit file> <mode> [quiet minutes]}
QUIET=$(( ${3:-45} * 60 ))
INTERVAL=${WATCH_INTERVAL:-60}
LOG=${EXIT_FILE%.exit}.log
BASE=$(basename "$EXIT_FILE" .exit)        # <config stem>_<mode>_<date>_<time>
STEM=${BASE%_"$MODE"_*_*}
NAME=clinego-$STEM-$MODE                    # run_mode.sh's container name

docker_cli() { nix shell nixpkgs#docker-client -c docker "$@"; }

# usage_usec of the run's container; empty when it is not running.
cpu_of() {
    local id
    id=$(docker_cli ps --no-trunc --format '{{.ID}} {{.Names}}' 2>/dev/null |
        awk -v n="$NAME" '$2 == n {print $1}')
    [ -n "$id" ] || return 0
    awk '$1 == "usage_usec" {print $2}' "/sys/fs/cgroup/system.slice/docker-$id.scope/cpu.stat" 2>/dev/null
}

echo "watch_run.sh: $NAME, log $LOG, quiet limit $((QUIET / 60)) min"
seen=""
mark_t=$(date +%s)
mark_cpu=""
gone=0
while :; do
    if [ -s "$EXIT_FILE" ]; then
        echo "watch_run.sh: run ended, exit $(cat "$EXIT_FILE")"
        exit 0
    fi
    now=$(date +%s)
    cpu=$(cpu_of)
    if [ -z "$cpu" ]; then
        gone=$((gone + 1))
        if [ "$gone" -ge 3 ]; then
            [ -s "$EXIT_FILE" ] && continue
            echo "watch_run.sh: no container $NAME for $gone checks and no exit file"
            exit 3
        fi
    else
        gone=0
        [ -n "$mark_cpu" ] || mark_cpu=$cpu
    fi
    mtime=$(date -r "$LOG" +%s 2>/dev/null || echo 0)
    if [ "$mtime" != "$seen" ]; then
        seen=$mtime
        mark_t=$now
        mark_cpu=$cpu
    elif [ $((now - mark_t)) -ge "$QUIET" ]; then
        gained=$(( (${cpu:-0} - ${mark_cpu:-0}) / 1000000 ))
        if [ "$gained" -lt 60 ]; then
            echo "watch_run.sh: STALL: log quiet $(( (now - mtime) / 60 )) min, container gained ${gained} CPU s in the last $(( (now - mark_t) / 60 )) min"
            tail -n 5 "$LOG"
            exit 2
        fi
        echo "$(date -Is) silent but busy: ${gained} CPU s in $(( (now - mark_t) / 60 )) min"
        mark_t=$now
        mark_cpu=$cpu
    fi
    sleep "$INTERVAL"
done
