#!/usr/bin/env bash
# mvp_autoscore.sh -- score each replicate as soon as mvp_run_sweep.sh finishes it.
#
#   benchmarks/mvp_autoscore.sh [SWEEP_LOG] [CELLS]
#
# Polls the sweep log for "[MVPxxxx] DONE" markers and runs mvp_score_sweep.sh on each seed
# once, in completion order. Scoring is cheap (minutes, single-digit GB) next to the fits, so
# overlapping it with the remaining seeds costs nothing and means the cross-seed aggregation
# can start the moment the last replicate lands.
#
# Exits when every seed in the manifest has been scored, or when the sweep driver is gone and
# no unscored DONE markers remain.
set -uo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
LOG="${1:-$PIPELINE_ROOT/benchmarks/mvp_eval/sweep_run4.log}"
CELLS="${2:-5}"
MANIFEST="$PIPELINE_ROOT/benchmarks/mvp_seeds.tsv"
STATE="$PIPELINE_ROOT/benchmarks/mvp_eval/autoscore.log"

total=$(awk 'NR>1 && NF' "$MANIFEST" | wc -l)
declare -A done_seeds
echo "$(date '+%F %T') autoscore START (expecting $total seeds)" | tee -a "$STATE"

while true; do
    while read -r s; do
        [[ -n "$s" && -z "${done_seeds[$s]:-}" ]] || continue
        done_seeds[$s]=1
        echo "$(date '+%F %T') scoring MVP$s" | tee -a "$STATE"
        "$PIPELINE_ROOT/benchmarks/mvp_score_sweep.sh" "$s" 1 "$CELLS" >> "$STATE" 2>&1 \
            && echo "$(date '+%F %T') scored MVP$s ($((${#done_seeds[@]}))/$total)" | tee -a "$STATE" \
            || echo "$(date '+%F %T') SCORING FAILED MVP$s" | tee -a "$STATE"
    done < <(grep -oE '^\[MVP[0-9]+\] DONE' "$LOG" 2>/dev/null | sed 's/^\[MVP//; s/\] DONE//')

    (( ${#done_seeds[@]} >= total )) && { echo "$(date '+%F %T') autoscore COMPLETE" | tee -a "$STATE"; break; }
    if ! pgrep -f "benchmarks/mvp_run_sweep.sh" >/dev/null 2>&1; then
        sleep 20
        n=$(grep -cE '^\[MVP[0-9]+\] DONE' "$LOG" 2>/dev/null || echo 0)
        (( n <= ${#done_seeds[@]} )) && {
            echo "$(date '+%F %T') autoscore STOPPING: sweep driver gone, ${#done_seeds[@]}/$total scored" | tee -a "$STATE"
            break; }
    fi
    sleep 30
done
