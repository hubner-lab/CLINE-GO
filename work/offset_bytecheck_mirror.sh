#!/usr/bin/env bash
# Build one scratch mirror for the offset byte-identity gate (2026-10-18, plan let-s-implement-this-plan-mighty-bubble.md).
# usage (host, repo root): bash work/offset_bytecheck_mirror.sh <default|minou> <pre|post>
# default = mirror of ANALYSIS/Trifolium86MainK5/maladaptationER (v1), minou = mirror of fig3/offset_minou.
# Inputs are relative symlinks to the real targets; paper_pack/ and site_maps/ are copied (snapshot subfolders excluded).
set -euo pipefail
cd /mnt/data/eugene/ADAPTOGENE
env=$1; tag=$2
M=work/offset_bytecheck_2026-10-18/${env}_${tag}
[ -e "$M" ] && { echo "exists: $M"; exit 1; }
mkdir -p "$M"
case $env in
  default) SRC=ANALYSIS/Trifolium86MainK5/maladaptationER; ARMS="$SRC/arm3 $SRC/arm3_ax3" ;;
  minou)   SRC=ANALYSIS/Trifolium86MainK5ERpc5/fig3/offset_minou
           ARMS="ANALYSIS/Trifolium86MainK5ERpc5/maladaptation/arm3 ANALYSIS/Trifolium86MainK5ERpc5/maladaptation/arm3_ax3" ;;
  *) echo "env?"; exit 1 ;;
esac
for a in $ARMS; do ln -s "$(realpath --relative-to="$M" "$a")" "$M/$(basename "$a")"; done
for s in paper_pack site_maps; do rsync -a --exclude '_20*/' "$SRC/$s/" "$M/$s/"; done
echo "$M"
