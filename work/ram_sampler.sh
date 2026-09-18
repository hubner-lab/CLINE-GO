#!/usr/bin/env bash
# Periodic host/container RSS trace. The guard only logs THRESHOLD CROSSINGS, so after b4's
# garden sweep there was no record of the actual peak -- which is the number needed to size
# the next stage's concurrency. This fills that gap: one line every 5 min, no notifications.
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE || exit 1
while true; do
    u=$(free -g | awk '/^Mem:/ {print $3}')
    n=$(nix shell nixpkgs#docker-client -c docker ps -q 2>/dev/null | wc -l)
    g=$(nix shell nixpkgs#docker-client -c docker ps --format '{{.Names}}' 2>/dev/null | grep -c gsweep)
    echo "$(date -Is) used=${u}G containers=${n} gsweep=${g}" >> work/ram_trace.log
    sleep 300
done
