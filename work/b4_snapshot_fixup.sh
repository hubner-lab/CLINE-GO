#!/usr/bin/env bash
# b4's block_close.sh is ALREADY RUNNING and cannot be edited (bash re-reads a running script
# by byte offset). When it reaches its snapshot step it will cp offset12/garden_qc.tsv and
# offset12/fitness_identity.tsv -- but b5's prep now runs concurrently and overwrites both
# (garden_qc is a fixed-name whole-file rewrite; the identity table grows 480 -> 600).
# Pristine copies of b4's correct content were taken before b5 was unblocked; this restores
# them afterwards and verifies the result.
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE || exit 1
SNAP=benchmarks/mvp_eval/offset12_ssclines_b4
LOG=work/b4_snapshot_fixup.log
log() { echo "$(date -Is) $*" | tee -a "$LOG"; }

log "waiting for b4 close-out to finish"
for i in $(seq 1 2880); do
    grep -q "close-out done" work/b4_close.nohup 2>/dev/null && break
    grep -q "GATE FAIL" work/b4_close.nohup 2>/dev/null && { log "b4 close-out failed; leaving snapshot alone"; exit 1; }
    sleep 30
done
grep -q "close-out done" work/b4_close.nohup 2>/dev/null || { log "timed out"; exit 1; }
sleep 5

# garden_qc: must contain b4's 120 seeds and nothing else.
b4seeds=$(cut -f1 work/seeds_b4.txt | sort -u | wc -l)
foreign=$(awk -F'\t' 'NR>1 {print $1}' "$SNAP/garden_qc.tsv" 2>/dev/null | sort -u \
          | grep -vxF -f <(sort -u work/seeds_b4.txt) | wc -l)
if (( foreign > 0 )); then
    log "garden_qc.tsv held $foreign non-b4 seed(s) -- restoring pristine"
    cp -p "$SNAP/garden_qc.b4_pristine.tsv" "$SNAP/garden_qc.tsv"
else
    log "garden_qc.tsv already correct ($b4seeds b4 seeds, 0 foreign)"
fi

# accumulated identity: b4's snapshot must show 480 rows (b1+b2+b3+b4), not 600.
rows=$(( $(wc -l < "$SNAP/fitness_identity.accumulated.tsv" 2>/dev/null || echo 1) - 1 ))
if (( rows != 480 )); then
    log "fitness_identity.accumulated.tsv had $rows rows (want 480) -- restoring pristine"
    cp -p "$SNAP/fitness_identity.accumulated.b4_pristine.tsv" "$SNAP/fitness_identity.accumulated.tsv"
else
    log "fitness_identity.accumulated.tsv already correct (480 rows)"
fi
log "fixup done; snapshot contents:"; ls -1 "$SNAP" | tee -a "$LOG"
