#!/usr/bin/env bash
# Sequential pipeline modes for config_Trifolium86MainKsearch.yaml (repo root); stop at first failure.
# usage: run_chain.sh <mode>...
set -u
CFG=config_Trifolium86MainKsearch.yaml
LOGDIR=/mnt/data/eugene/ADAPTOGENE/ANALYSIS/Trifolium86MainKsearch/logs
cd /mnt/data/eugene/ADAPTOGENE || exit 1
for M in "$@"; do
  echo "$(date +%T) START $M"
  nix shell nixpkgs#docker-client -c docker run --name adaptogene-run-Trifolium86MainKsearch \
    --user "$(id -u):$(id -g)" --rm -e USER=adaptogene -e OPENBLAS_NUM_THREADS=16 \
    --cpus=120 --memory=400g -v /mnt/data/eugene/ADAPTOGENE:/pipeline adaptogene:latest \
    snakemake -c90 -s Snakefile --config mode="$M" --configfile "$CFG" --scheduler greedy \
    > "$LOGDIR/$M.log" 2>&1
  rc=$?
  echo "$(date +%T) END $M rc=$rc"
  [ "$rc" -eq 0 ] || { echo "STOP at $M"; exit "$rc"; }
done
echo "ALL DONE"
