#!/usr/bin/env bash
# 12_fetch_associated_snps.sh -- copy the run's associated SNPs (GEA/tables/selected_snps.tsv, the 30 SNPs of the
# EMMAX-RDA agreement set) from lab, for the distance column of the genes table (13_table_s9.py).
#
# Lab is treated as read-only: scp copies the file back, the ssh command only prints its md5 sum.
#
# Usage: bash 12_fetch_associated_snps.sh <data_dir>
# Writes: <data_dir>/lab/{selected_snps.tsv, selected_snps_md5.txt}
set -euo pipefail

DATA=${1:?usage: 12_fetch_associated_snps.sh <data_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"

OUT="$DATA/lab"
mkdir -p "$OUT"
SRC="$LAB_RESULTS/GEA/tables/selected_snps.tsv"

echo "[12] copying the associated SNPs"
scp -p -o BatchMode=yes "$LAB_HOST:$SRC" "$OUT/"
ssh -o BatchMode=yes "$LAB_HOST" "md5sum $SRC" > "$OUT/selected_snps_md5.txt"
local_md5=$(md5sum "$OUT/selected_snps.tsv" | cut -d' ' -f1)
[[ "$local_md5" == "$(cut -d' ' -f1 "$OUT/selected_snps_md5.txt")" ]] || { echo "[12] md5 differs from lab" >&2; exit 1; }
echo "[12] done (md5 $local_md5)"
