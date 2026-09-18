#!/usr/bin/env bash
# Stage G + H -- b5 from the end of the GEA sweep through the short report.
#
# Waits for BOTH GEA waves, then steps 6-8 (which themselves wait for b4's garden_qc
# snapshot -- the one corpus-wide serialisation point in step 7), then the garden sweep,
# then close-out and the report.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

COHORT=ssclines_b5
TAG=ssclines_nequal_mbreaks
PREV=ssclines_b4

die() { echo "CHAIN STOP: $*  $(date -Is)"; exit 1; }

# ---- wait for both GEA waves
echo "=== b5 chain: waiting for both GEA waves  $(date -Is)"
for i in $(seq 1 5760); do
    d1=$(grep -c "\] DONE " logs_gea_ssclines_b5.h1.log 2>/dev/null); d1=${d1:-0}
    d2=$(grep -c "\] DONE " logs_gea_ssclines_b5.h2.log 2>/dev/null); d2=${d2:-0}
    (( d1 + d2 >= 120 )) && break
    sleep 60
done
d1=$(grep -c "\] DONE " logs_gea_ssclines_b5.h1.log 2>/dev/null); d1=${d1:-0}
d2=$(grep -c "\] DONE " logs_gea_ssclines_b5.h2.log 2>/dev/null); d2=${d2:-0}
echo "GEA DONE lines: h1=$d1 h2=$d2 total=$((d1+d2))"
(( d1 + d2 >= 120 )) || die "GEA incomplete ($((d1+d2))/120) -- see INCOMPLETE lines in the wave logs"

np=$(ls -d benchmarks/mvp_eval/params/MVP*/c1 2>/dev/null | wc -l)
echo "step 5 gate: params/MVP*/c1 dirs on disk = $np"

# ---- steps 6, 7, 8 (waits internally for b4's garden_qc snapshot)
bash work/block_prep.sh $COHORT $PREV > work/b5_prep.nohup 2>&1
grep -q "prep done" work/b5_prep.nohup || { tail -20 work/b5_prep.nohup; die "b5 prep"; }
echo "=== b5 prep done  $(date -Is)"

# ---- step 9, two waves.
#
# WIDTH, NOT DEPTH. mvp_garden_sweep.sh:89 passes NO --cpus, only --memory, so cores were
# never capped -- the limiters are `snakemake -cN` per seed and how many seeds run at once.
# Intra-seed parallelism scales badly (b3 -c8: 1.51 h/seed vs b1 -c2: 2.05 h/seed -- 4x the
# cores for 1.36x the speed), so the machine goes into JOBS. b4 ran 32 concurrent for 21 h
# WHILE sharing the box with 28 GEA containers and never crossed the guard's 780 G soft
# threshold. b5's garden runs alone.
GJOBS=${GJOBS:-40}          # per wave -> 80 concurrent seeds
GSCORES=${GSCORES:-4}       # snakemake -c per seed
BASE=$(free -g | awk '/^Mem:/ {print $3}')
echo "--- baseline before wave 1: ${BASE}G used; JOBS=$GJOBS/wave snake -c$GSCORES  $(date -Is)"
nohup bash work/garden_wave.sh $COHORT h1 $GJOBS $GSCORES > work/${COHORT}_garden_h1.nohup 2>&1 &
W1=$!

# Wave-2 gate is a PROJECTION, not a fixed number. A fixed 600 G gate is wrong in both
# directions: it blocks a cheap wave 1 needlessly, and it waves through an expensive one
# straight into the guard's 880 G hard stop. Measure what wave 1 actually costs, double it,
# and only open if the projected total clears the 780 G soft threshold.
PLAN_CEIL=${PLAN_CEIL:-780}
seen_done=0
for i in $(seq 1 2880); do
    if (( seen_done == 0 )); then
        grep -q "\] DONE " "logs_garden_${COHORT}.h1.log" 2>/dev/null && {
            seen_done=1; echo "    wave-1 first DONE seen  $(date -Is)"; }
        kill -0 $W1 2>/dev/null || { echo "WARNING: wave 1 exited before any DONE"; seen_done=1; }
    fi
    if (( seen_done == 1 )); then
        u=$(free -g | awk '/^Mem:/ {print $3}')
        cost=$(( u - BASE )); (( cost < 0 )) && cost=0
        proj=$(( u + cost ))
        if (( proj < PLAN_CEIL )); then
            echo "    RAM gate open: wave1 cost=${cost}G, now=${u}G, projected=${proj}G < ${PLAN_CEIL}G"
            break
        fi
        (( i % 20 == 0 )) && echo "    holding wave 2: projected ${proj}G >= ${PLAN_CEIL}G (now ${u}G, wave1 cost ${cost}G)  $(date -Is)"
    fi
    sleep 30
done
echo "--- host at wave-2 launch: $(free -g | awk '/^Mem:/ {print $3"G used"}')  $(date -Is)"
nohup bash work/garden_wave.sh $COHORT h2 $GJOBS $GSCORES > work/${COHORT}_garden_h2.nohup 2>&1 &
W2=$!
wait $W1 $W2
echo "=== b5 garden sweep both waves finished  $(date -Is)"

# ---- step 9 gate
echo "--- step 9 gate"
short=0
while read -r s; do
    n=$(find "benchmarks/mvp_eval/offset12/gardens/$s" -name '*.tsv' 2>/dev/null | wc -l)
    d=$(grep -c -P "^MVP${s}\t" panel_underfilled_${COHORT}.tsv 2>/dev/null); d=${d:-0}
    want=$(( 3584 - 448 * d ))
    (( n != want )) && { echo "  SHORT MVP$s: $n (want $want, $d panel(s) dropped)"; short=1; }
done < work/seeds_b5.txt
(( short == 1 )) && die "step 9 gate"
echo "step 9 gate passed"

# ---- steps 10, 11, snapshot
bash work/block_close.sh $COHORT $PREV > work/b5_close.nohup 2>&1
grep -q "close-out done" work/b5_close.nohup || { tail -20 work/b5_close.nohup; die "b5 close-out"; }
echo "=== b5 close-out done  $(date -Is)"

# ---- step 12 + short report
bash work/block_report.sh $COHORT $TAG 120 40 > work/b5_report.nohup 2>&1
tail -25 work/b5_report.nohup
echo "=== b5 CHAIN COMPLETE -- all five blocks scored; pooled report is the remaining step  $(date -Is)"
