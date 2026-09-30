#!/usr/bin/env bash
# Phase 2b memory sampler: every 15 s, the resident memory of every mvp-sweep-* container
# (docker stats, page cache excluded), one row per container per tick. Container names carry
# the seed and mode (mvp-sweep-MVP{seed}-{mode}), so per-mode peaks fall out of the table.
# Exits once no mvp-sweep-* container and no driver has been seen for 10 minutes.
#
#   setsid nohup work/phase2b_memtrace.sh > /dev/null 2>&1 &
set -uo pipefail
OUT=/mnt/data/eugene/ADAPTOGENE/benchmarks/mvp_eval/onefit_run/memtrace.tsv
DOCKER=$(nix shell nixpkgs#docker-client -c sh -c 'command -v docker')
[[ -s "$OUT" ]] || printf "ts\tcontainer\tmem_gib\tcap_gib\n" > "$OUT"
idle=0
while (( idle < 40 )); do
    ts=$(date -Is)
    rows=$("$DOCKER" stats --no-stream --format '{{.Name}}\t{{.MemUsage}}' 2>/dev/null | grep '^mvp-sweep-')
    if [[ -n "$rows" ]]; then
        idle=0
        awk -F'\t' -v t="$ts" '
            function gib(x,  v, u) { v = x + 0; u = x; gsub(/[0-9.]/, "", u)
                if (u ~ /^Ki/) return v / 1048576; if (u ~ /^Mi/) return v / 1024
                if (u ~ /^Gi/) return v; if (u ~ /^Ti/) return v * 1024; return v / 1073741824 }
            { split($2, m, " / "); printf "%s\t%s\t%.3f\t%.1f\n", t, $1, gib(m[1]), gib(m[2]) }' \
            <<< "$rows" >> "$OUT"
    elif pgrep -f "mvp_run_sweep.sh" > /dev/null; then
        idle=0
    else
        idle=$((idle + 1))
    fi
    sleep 15
done
