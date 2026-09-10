#!/usr/bin/env bash
# =============================================================================
# ram_watchdog.sh -- kill the biggest RAM consumer when host memory crosses a limit.
#
#   benchmarks/ram_watchdog.sh [LIMIT_GB] [SCOPE] [INTERVAL_S]
#
#   LIMIT_GB    breach threshold for USED memory (default 920, host total ~1007)
#   SCOPE       which containers may be killed (default: sweep)
#                 sweep  -- only containers this sweep launched (snakemake + benchmark
#                           Rscripts started from benchmarks/).
#                 mine   -- any docker container running under THIS uid, which adds the
#                           user's own Shiny/dev containers.
#   INTERVAL_S  poll interval (default 15)
#
# THIS IS A SHARED MACHINE -- THE HARD SAFETY RULE
# ------------------------------------------------
# At the time of writing, ~390 GB of the host's 1007 GB belongs to ANOTHER user: eleven
# `minimap2-nd` processes from /mnt/data/guyh/programs/NextDenovo/ running as `nobody`.
# That is a long-running genome assembly and killing it would destroy days of someone
# else's work.
#
# Two structural guarantees, not conventions:
#   1. Only DOCKER CONTAINERS are ever eligible. Bare host processes are never touched,
#      full stop -- guyh's assembly is bare processes and therefore unreachable by this
#      script under any scope, threshold, or argument.
#   2. A container is eligible only if its configured user matches this uid.
# The threshold exists to stop OUR footprint from OOM-killing the machine, not to reclaim
# memory from whoever else is on it.
#
# WHAT "USED" MEANS HERE
# ----------------------
# Breach is `MemTotal - MemAvailable >= LIMIT_GB`, read from /proc/meminfo. MemAvailable is
# the kernel's own estimate of what can be handed out without swapping, so it already
# discounts reclaimable page cache -- this is the number that actually predicts an OOM kill.
# Using the `used` column from `free` would under-count, and using `free` would over-count by
# treating all cache as lost.
#
# BEHAVIOUR ON BREACH
# -------------------
# Stop the single largest-RSS eligible container, wait, re-measure. Repeat up to
# MAX_KILLS_PER_CYCLE times per poll so a runaway cannot take the whole fleet down in one
# sweep of the loop. `docker stop` (SIGTERM, 10s grace) is used rather than kill -9 so
# snakemake gets the chance to unlink partial outputs; the sweep driver's --rerun-incomplete
# handles whatever it does not.
#
# Never eligible: the docker daemon, this script, the Claude session, PID 1, and any
# container not matching SCOPE.
# =============================================================================
set -uo pipefail

LIMIT_GB="${1:-920}"
SCOPE="${2:-sweep}"
INTERVAL_S="${3:-15}"
MAX_KILLS_PER_CYCLE="${MAX_KILLS_PER_CYCLE:-2}"

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
LOG="$PIPELINE_ROOT/benchmarks/mvp_eval/ram_watchdog.log"
DOCKER=(nix shell nixpkgs#docker-client -c docker)

mkdir -p "$(dirname "$LOG")"
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

MY_UID="$(id -u)"

case "$SCOPE" in
    sweep|mine) ;;
    *) say "FATAL: SCOPE must be 'sweep' or 'mine', got '$SCOPE'"; exit 2 ;;
esac

used_gb() {
    awk '/MemTotal:/{t=$2} /MemAvailable:/{a=$2} END{printf "%d", (t-a)/1048576}' /proc/meminfo
}
total_gb() { awk '/MemTotal:/{printf "%d", $2/1048576}' /proc/meminfo; }

# Eligible containers, largest resident first: "<id> <gb> <command>"
candidates() {
    "${DOCKER[@]}" stats --no-stream --format '{{.ID}}\t{{.MemUsage}}' 2>/dev/null \
    | while IFS=$'\t' read -r id mem; do
        raw="${mem%% /*}"
        gb=$(awk -v v="$raw" 'BEGIN{
            u=v; sub(/[0-9.]+/,"",u); n=v; sub(/[A-Za-z]+$/,"",n);
            if (u=="GiB") printf "%.1f", n;
            else if (u=="MiB") printf "%.1f", n/1024;
            else if (u=="KiB") printf "%.1f", n/1048576;
            else printf "0"}')
        cmd=$("${DOCKER[@]}" ps --filter "id=$id" --format '{{.Command}}' 2>/dev/null)

        # GUARD 1 (ownership): the container must have been started as this uid. Anything
        # started by another user, or with no explicit --user, is skipped. Combined with
        # the fact that only containers are ever considered, another user's work cannot be
        # reached from here.
        cuser=$("${DOCKER[@]}" inspect --format '{{.Config.User}}' "$id" 2>/dev/null)
        [[ "${cuser%%:*}" == "$MY_UID" ]] || continue

        # GUARD 2 (scope): what this sweep started -- snakemake runs and benchmark
        # Rscripts. A bare `Rscript /pipeline/scripts/...` is a user-launched job (the
        # Shiny app, dev.R) and is excluded unless SCOPE=mine.
        if [[ "$SCOPE" == "sweep" ]]; then
            [[ "$cmd" == *snakemake* || "$cmd" == *"/pipeline/benchmarks/"* ]] || continue
        fi
        printf "%s\t%s\t%s\n" "$id" "$gb" "$cmd"
      done | sort -k2 -gr
}

say "watchdog START limit=${LIMIT_GB}GB scope=${SCOPE} interval=${INTERVAL_S}s host_total=$(total_gb)GB"
say "eligible now:"; candidates | awk -F'\t' '{printf "    %s  %sGB  %.60s\n",$1,$2,$3}' | tee -a "$LOG"

while true; do
    u=$(used_gb)
    if (( u >= LIMIT_GB )); then
        say "BREACH used=${u}GB >= ${LIMIT_GB}GB"
        killed=0
        while (( killed < MAX_KILLS_PER_CYCLE )); do
            line=$(candidates | head -1)
            if [[ -z "$line" ]]; then say "  no eligible container to stop (scope=$SCOPE)"; break; fi
            id=$(cut -f1 <<<"$line"); gb=$(cut -f2 <<<"$line"); cmd=$(cut -f3 <<<"$line")
            say "  stopping $id (${gb}GB) :: ${cmd:0:80}"
            "${DOCKER[@]}" stop -t 10 "$id" >/dev/null 2>&1 \
                && say "  stopped $id" || say "  FAILED to stop $id"
            killed=$((killed+1))
            sleep 5
            u=$(used_gb)
            say "  used now ${u}GB"
            (( u < LIMIT_GB )) && break
        done
    fi
    sleep "$INTERVAL_S"
done
