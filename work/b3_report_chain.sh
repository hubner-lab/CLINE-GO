#!/usr/bin/env bash
# Stage E -- b3 step 12 + short report, once b3's close-out lands.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1
echo "=== b3 report chain: waiting for close-out  $(date -Is)"
for i in $(seq 1 2880); do
    grep -q "b3 close-out done" work/b3_close.nohup 2>/dev/null && break
    grep -q "GATE FAIL" work/b3_close.nohup 2>/dev/null && { echo "CHAIN STOP: b3 close-out gate failed"; exit 1; }
    sleep 30
done
grep -q "b3 close-out done" work/b3_close.nohup || { echo "CHAIN STOP: close-out never finished"; exit 1; }
bash work/block_report.sh ssclines_b3 ssclines_nequal_mconst 120 40
echo "=== b3 REPORT COMPLETE  $(date -Is)"
