#!/usr/bin/env bash
# diffpoll.sh <container-name-regex> <out> <max_wait_s> -- poll `docker diff` every 2 s on every
# running container whose name matches, until none has matched for 30 s after the first sighting.
set -uo pipefail
re=$1; out=$2; maxw=${3:-3600}; : > "$out"
D=(nix shell --offline nixpkgs#docker-client -c docker)
polls=0; seen=0; idle=0; t0=$(date +%s)
while :; do
  names=$("${D[@]}" ps --format '{{.Names}}' | grep -E "$re" || true)
  if [[ -n "$names" ]]; then seen=1; idle=0
    for n in $names; do "${D[@]}" diff "$n" 2>/dev/null | sed "s|^|$n |" >> "$out"; done; polls=$((polls+1))
  else idle=$((idle+2)); fi
  (( seen && idle >= 30 )) && break
  (( $(date +%s) - t0 > maxw )) && break
  sleep 2
done
sort -u "$out" -o "$out"
echo "polls: $polls; distinct layer entries: $(awk '{print $2,$3}' "$out" | sort -u | wc -l)"
awk '{print $2,$3}' "$out" | sort -u
