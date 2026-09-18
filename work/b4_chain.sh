#!/usr/bin/env bash
# Stage D + F -- b4 from step 9 through the short report. Auto-chained with hard gates.
#
# Runs AFTER work/b4_prep.sh (steps 6-8) reports done. Launched detached; each stage stops
# the chain on a gate failure and leaves its logs. b5's chain is independent and keeps going.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

COHORT=ssclines_b4
TAG=ssclines_ncline_ctredge
PREV=ssclines_b3

die() { echo "CHAIN STOP: $*  $(date -Is)"; exit 1; }

# ---- wait for prep (steps 6-8)
echo "=== b4 chain: waiting for b4 prep  $(date -Is)"
for i in $(seq 1 2880); do
    grep -q "b4 prep done\|=== b4 prep done" work/b4_prep.nohup 2>/dev/null && break
    grep -q "GATE FAIL" work/b4_prep.nohup 2>/dev/null && die "b4 prep gate failed"
    sleep 30
done
grep -q "prep done" work/b4_prep.nohup 2>/dev/null || die "b4 prep never finished"
echo "=== b4 prep done, starting garden sweep  $(date -Is)"

# ---- step 9, wave 1
nohup bash work/garden_wave.sh $COHORT h1 > work/${COHORT}_garden_h1.nohup 2>&1 &
W1=$!

# Wave 2 needs TWO conditions, not one.
#
#  (a) wave 1's FIRST DONE line -- the first point past a completed Gradient Forest fit,
#      which is the RSS peak. A fixed clock would sample the wrong phase.
#  (b) host RAM below RAM_GATE. This block's garden sweep overlaps b5's GEA sweep, which
#      alone sits at ~420 GB across 28 containers. 32 garden seeds on top of that would
#      project past the 900 GB budget and the guard would start stopping containers.
#      Waiting here costs wall time; getting it wrong costs killed work.
RAM_GATE=${RAM_GATE:-600}
echo "--- wave 2 gate: wave-1 first DONE, and host RAM < ${RAM_GATE}G  $(date -Is)"
seen_done=0
for i in $(seq 1 2880); do
    if (( seen_done == 0 )); then
        grep -q "\] DONE " "logs_garden_${COHORT}.h1.log" 2>/dev/null && {
            seen_done=1; echo "    wave-1 first DONE seen  $(date -Is)"; }
        kill -0 $W1 2>/dev/null || { echo "WARNING: wave 1 exited before any DONE"; seen_done=1; }
    fi
    if (( seen_done == 1 )); then
        u=$(free -g | awk '/^Mem:/ {print $3}')
        (( u < RAM_GATE )) && { echo "    RAM gate open: ${u}G < ${RAM_GATE}G"; break; }
        (( i % 20 == 0 )) && echo "    holding wave 2: host ${u}G >= ${RAM_GATE}G  $(date -Is)"
    fi
    sleep 30
done
echo "--- host at wave-2 launch: $(free -g | awk '/^Mem:/ {print $3"G used"}')  $(date -Is)"
nohup bash work/garden_wave.sh $COHORT h2 > work/${COHORT}_garden_h2.nohup 2>&1 &
W2=$!

wait $W1 $W2
echo "=== b4 garden sweep both waves finished  $(date -Is)"

# ---- step 9 gate: 8 panels x 4 methods x 112 gardens = 3584 files per seed,
#      minus 448 per panel dropped at step 6.
echo "--- step 9 gate"
short=0
while read -r s; do
    n=$(find "benchmarks/mvp_eval/offset12/gardens/$s" -name '*.tsv' 2>/dev/null | wc -l)
    # grep -c exits 1 on no match but still prints 0 -- do NOT add `|| echo 0`, that yields "0 0"
    d=$(grep -c -P "^MVP${s}\t" panel_underfilled_${COHORT}.tsv 2>/dev/null)
    d=${d:-0}
    want=$(( 3584 - 448 * d ))
    (( n != want )) && { echo "  SHORT MVP$s: $n (want $want, $d panel(s) dropped)"; short=1; }
done < work/seeds_b4.txt
(( short == 1 )) && die "step 9 gate"
echo "step 9 gate passed"

# ---- steps 10, 11, snapshot (takes the shared offset12 scoring lock)
bash work/block_close.sh $COHORT $PREV > work/b4_close.nohup 2>&1
grep -q "close-out done" work/b4_close.nohup || { tail -20 work/b4_close.nohup; die "b4 close-out"; }
echo "=== b4 close-out done  $(date -Is)"

# ---- step 12 + short report
bash work/block_report.sh $COHORT $TAG 120 40 > work/b4_report.nohup 2>&1
tail -25 work/b4_report.nohup
echo "=== b4 CHAIN COMPLETE  $(date -Is)"
