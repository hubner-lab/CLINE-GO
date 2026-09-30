#!/usr/bin/env bash
# Phase 2b autoscaler for the gated driver (work/phase2b_run_gated.sh). Every 120 s:
#   seed_jobs  = seeds driver C holds in flight + free CPU room, where room is measured, not
#                modelled: (CPU_BUDGET - 5-min load average) / 3 new seeds, capped at +6 per step
#                so a burst of fresh launches (which ramp up slowly) cannot overshoot. The first
#                version costed containers at 3.5 CPUs upstream / 2 in gea and left ~60 CPUs idle
#                (load 119 at 30 gea + 34 upstream). Seeds waiting for a gea slot cost nothing,
#                so seed_jobs rises as they queue.
#   gea_slots  = memory brake only: if the latest memtrace tick exceeds MEM_BRAKE_GIB, slots are
#                lowered to (gea containers - 1) so no further seed enters gea; restored to
#                GEA_SLOTS_MAX once the total is back under MEM_OK_GIB. Otherwise left alone, so
#                a hand-set value is respected.
# Exits once driver C has finished and no mvp-sweep container is left. Log: onefit_run/autoscale.log
#
#   setsid nohup work/phase2b_autoscale.sh > /dev/null 2>&1 &
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE
R=benchmarks/mvp_eval/onefit_run
G=benchmarks/mvp_eval/params_onefit/.gea_gate
LOG=$R/autoscale.log
CPU_BUDGET="${CPU_BUDGET:-180}"
MEM_BRAKE_GIB="${MEM_BRAKE_GIB:-800}"
MEM_OK_GIB="${MEM_OK_GIB:-650}"
GEA_SLOTS_MAX="${GEA_SLOTS_MAX:-40}"
DOCKER=$(nix shell nixpkgs#docker-client -c sh -c 'command -v docker')
braked=0

while :; do
    names=$("$DOCKER" ps --format '{{.Names}}' 2>/dev/null | grep '^mvp-sweep-' || true)
    n_gea=$(grep -c -- '-gea$' <<< "$names" || true)
    n_up=$(grep -vc -- '-gea$' <<< "$names" || true); [[ -z "$names" ]] && n_up=0
    started=$(grep -oE '^\[MVP[0-9]+\] k_best=' $R/driver_C.log | sort -u | wc -l)
    ended=$(grep -oE '^\[MVP[0-9]+\] (DONE|FAILED|FATAL)' $R/driver_C.log | grep -oE 'MVP[0-9]+' | sort -u | wc -l)
    inflight=$((started - ended))
    load1=$(cut -d' ' -f2 /proc/loadavg)   # 5-min: the 1-min swung 163-211 around a 180 target
    room=$(awk -v b="$CPU_BUDGET" -v l="$load1" 'BEGIN{r=(b - l)/3; if (r<0) r=0; if (r>6) r=6; printf "%d", r}')
    jobs=$((inflight + room))
    # Over budget: set the limit BELOW the current count, so seeds that finish are not replaced
    # (at jobs == inflight every finished seed is backfilled and the load never comes down).
    awk -v b="$CPU_BUDGET" -v l="$load1" 'BEGIN{exit !(l > b + 10)}' && jobs=$((inflight - 3))
    (( jobs < 10 )) && jobs=10; (( jobs > 90 )) && jobs=90
    echo "$jobs" > $G/seed_jobs

    last_ts=$(tail -n 1 $R/memtrace.tsv | cut -f1)
    mem=$(awk -F'\t' -v t="$last_ts" '$1==t {s+=$3} END{printf "%d", s}' $R/memtrace.tsv)
    if (( mem > MEM_BRAKE_GIB )); then
        s=$(( n_gea > 1 ? n_gea - 1 : 1 )); echo "$s" > $G/gea_slots; braked=1
    elif (( braked && mem < MEM_OK_GIB )); then
        echo "$GEA_SLOTS_MAX" > $G/gea_slots; braked=0
    fi
    echo -e "$(date -Is)\tload1=$load1\tgea=$n_gea\tupstream=$n_up\tC_inflight=$inflight\tseed_jobs=$jobs\tmem_gib=$mem\tgea_slots=$(cat $G/gea_slots)\tbraked=$braked" >> $LOG

    if ! pgrep -f "work/phase2b_run_gated.sh" > /dev/null && [[ -z "$names" ]]; then break; fi
    sleep 120
done
