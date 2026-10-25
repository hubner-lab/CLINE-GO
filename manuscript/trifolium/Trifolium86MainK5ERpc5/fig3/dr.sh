#!/usr/bin/env bash
# One build_fig3.R call in the image (host helper). Usage: dr.sh <only-prefix|all> ; log in fig3/logs/build_<prefix>.log
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE
P=${1:?prefix}; ARG=""; [ "$P" != all ] && ARG="only=$P"
nix shell nixpkgs#docker-client -c docker run --rm --name "fig3_build_${P//[^a-z0-9]/}" --user "$(id -u):$(id -g)" \
  -e USER=adaptogene -e OPENBLAS_NUM_THREADS=16 -v /mnt/data/eugene/ADAPTOGENE:/pipeline -w /pipeline adaptogene:latest \
  Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/build_fig3.R $ARG > "ANALYSIS/Trifolium86MainK5ERpc5/fig3/logs/build_${P}.log" 2>&1
echo "exit $?"; tail -5 "ANALYSIS/Trifolium86MainK5ERpc5/fig3/logs/build_${P}.log"
