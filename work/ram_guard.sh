#!/usr/bin/env bash
# Host RAM guard for the SS-Clines run.
#
# WHY IT EXISTS. The only host-level watchdog in this corpus is inside mvp_garden_sweep.sh
# (:56-66), so it covers the GARDEN sweep only. The GEA sweep (mvp_run_sweep.sh) has none:
# its sole guard is the per-container `--memory=40g` cap, and 28 concurrent GEA containers
# cap out at 1120 GB -- above the 900 GB budget for this host. This process covers every
# stage, not just the garden sweep.
#
# It is also strictly better behaved than the built-in one, which stops
# `docker ps --filter name=gsweep_ | head -1` -- an ARBITRARY container, and one that can
# belong to the other half's driver. This picks the largest-RSS victim and logs it first.
#
# SAFETY. Only containers running IMAGE are ever candidates, and the scoring/tables
# containers are exempt: they are small (7 GB), they hold the offset12 serialisation lock,
# and killing one would corrupt a block's close-out. Every action is logged with the full
# candidate list so a wrong victim is auditable after the fact.
#
# Consequence of a stop: mvp_run_sweep.sh's snake() sees a nonzero exit, marks the seed
# FAILED at that mode, and the driver reports it as INCOMPLETE at the end. The seed is
# re-runnable -- Snakemake's DAG resumes, and the driver clears the stale lock at seed start.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

IMAGE=cline-go:latest
SOFT=${SOFT:-780}          # GB host used -- warn
HARD=${HARD:-880}          # GB host used -- stop the largest of our containers
INTERVAL=${INTERVAL:-20}
LOG=work/ram_guard.log
DK=(nix shell nixpkgs#docker-client -c docker)

log() { echo "$(date -Is) $*" | tee -a "$LOG"; }

log "GUARD START soft=${SOFT}G hard=${HARD}G interval=${INTERVAL}s image=$IMAGE"
warned=0
while true; do
    used=$(free -g | awk '/^Mem:/ {print $3}')
    [[ -z "$used" ]] && { sleep "$INTERVAL"; continue; }

    if (( used >= HARD )); then
        # Candidates: running containers of our image, excluding the scoring/tables//prep
        # singletons. Sorted by resident size, biggest first.
        mapfile -t cand < <("${DK[@]}" ps --filter "ancestor=$IMAGE" \
                              --format '{{.ID}}\t{{.Names}}' 2>/dev/null \
                            | grep -vE 'mvp-(score|tables|fitness|identity|snpsets|gcfg|genv)')
        if (( ${#cand[@]} == 0 )); then
            log "HARD ${used}G but no eligible container to stop -- nothing done"
            sleep "$INTERVAL"; continue
        fi
        ids=$(printf '%s\n' "${cand[@]}" | cut -f1 | paste -sd' ')
        victim=$("${DK[@]}" stats --no-stream --format '{{.ID}}\t{{.Name}}\t{{.MemUsage}}' $ids 2>/dev/null \
                 | awk -F'\t' '{split($3,a," "); v=a[1];
                                if (v ~ /GiB/) {sub(/GiB/,"",v)}
                                else if (v ~ /MiB/) {sub(/MiB/,"",v); v=v/1024}
                                else if (v ~ /KiB/) {sub(/KiB/,"",v); v=v/1048576}
                                printf "%.2f\t%s\t%s\n", v, $1, $2}' \
                 | sort -rn | head -1)
        vsize=$(cut -f1 <<<"$victim"); vid=$(cut -f2 <<<"$victim"); vname=$(cut -f3 <<<"$victim")
        log "HARD THRESHOLD ${used}G >= ${HARD}G -- candidates=${#cand[@]} stopping ${vname} (${vid}, ${vsize} GiB)"
        printf '%s\n' "${cand[@]}" >> "$LOG"
        "${DK[@]}" stop "$vid" >/dev/null 2>&1 \
            && log "STOPPED ${vname} (${vid}) -- its seed will report INCOMPLETE and is re-runnable" \
            || log "FAILED to stop ${vname} (${vid})"
        sleep 60          # let the memory actually come back before re-evaluating
        continue
    fi

    if (( used >= SOFT )); then
        if (( warned == 0 )); then
            log "SOFT THRESHOLD ${used}G >= ${SOFT}G -- $("${DK[@]}" ps -q 2>/dev/null | wc -l) containers running"
            warned=1
        fi
    else
        (( warned == 1 )) && log "recovered: ${used}G < ${SOFT}G"
        warned=0
    fi
    sleep "$INTERVAL"
done
