#!/usr/bin/env bash
# mvp_unstick.sh -- external supervisor for the garden sweep.
#
# WHY: benchmarks/mvp_garden_run_par.sh ends with a bare `wait`, which also waits for its own
# RAM-watchdog subshell -- an infinite loop. So the driver never exits after its last lane
# finishes, and benchmarks/mvp_overnight.sh (which blocks on the driver PID) never starts the
# next pass. Measured cost: ~70 minutes of an idle machine tonight.
#
# The driver script cannot be patched while it is executing (bash re-reads a running script by
# byte offset), so this supervisor fixes the symptom from outside: when no lane containers have
# been alive for GRACE consecutive checks but a driver process still exists, kill the driver.
# The orchestrator then proceeds to the next pass, the gap-fill, and scoring on its own.
#
# Safe by construction: it only ever kills mvp_garden_run_par.sh processes, and only when zero
# gp09_* containers exist -- i.e. when there is provably no lane work in flight.
set -uo pipefail

DOCKER=(nix shell nixpkgs#docker-client -c docker)
LOG="${LOG:-/mnt/data/eugene/CLINE-GO/benchmarks/mvp_eval/offset09/runlogs/unstick.log}"
INTERVAL="${INTERVAL:-120}"
GRACE="${GRACE:-3}"          # consecutive empty checks before acting (~6 min)

idle=0
while true; do
    n_containers=$("${DOCKER[@]}" ps --filter "name=gp09" -q 2>/dev/null | wc -l)
    n_driver=$(pgrep -fc "mvp_garden_run_par.sh" 2>/dev/null || echo 0)

    if [[ "$n_containers" -eq 0 && "$n_driver" -gt 0 ]]; then
        idle=$((idle + 1))
        echo "[$(date +%H:%M:%S)] no lane containers, driver alive ($n_driver) -- idle=$idle/$GRACE" >> "$LOG"
        if [[ "$idle" -ge "$GRACE" ]]; then
            echo "[$(date +%H:%M:%S)] killing stuck driver so the orchestrator can advance" >> "$LOG"
            pkill -9 -f "mvp_garden_run_par.sh"
            idle=0
        fi
    else
        idle=0
    fi

    # stop supervising once the orchestrator has finished everything
    if grep -q "OVERNIGHT RUN COMPLETE" \
        /mnt/data/eugene/CLINE-GO/benchmarks/mvp_eval/offset09/runlogs/overnight.log 2>/dev/null; then
        echo "[$(date +%H:%M:%S)] orchestrator reported completion -- supervisor exiting" >> "$LOG"
        exit 0
    fi
    sleep "$INTERVAL"
done
