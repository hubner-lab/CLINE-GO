#!/usr/bin/env bash
# 01_fetch_lab_inputs.sh -- copy the pipeline region tables and the BRAKER GTF from lab, and translate the
# BRAKER transcripts on lab with gffread (protein FASTA streamed to stdout, saved locally).
#
# Lab is treated as read-only: scp copies files back, every ssh command only reads and prints to stdout
# (gffread -o /dev/null discards its GFF output; the assembly's .fai index already exists, so neither
# gffread nor samtools writes one).
#
# Usage: bash 01_fetch_lab_inputs.sh <data_dir>
# Writes: <data_dir>/lab/{genes_per_region.tsv, regions_combined.tsv, Trifolium_genes.gtf,
#         braker_all_proteins.faa, gffread.stderr, lab_md5.txt, lab_tool_versions.txt}
set -euo pipefail

DATA=${1:?usage: 01_fetch_lab_inputs.sh <data_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"

OUT="$DATA/lab"
mkdir -p "$OUT"
SSH=(ssh -o BatchMode=yes "$LAB_HOST")

echo "[01] copying pipeline tables and the BRAKER GTF"
scp -p -o BatchMode=yes "$LAB_HOST:$LAB_RESULTS/GEA/tables/genes_per_region.tsv" "$OUT/"
scp -p -o BatchMode=yes "$LAB_HOST:$LAB_RESULTS/GEA/tables/regions_combined.tsv" "$OUT/"
scp -p -o BatchMode=yes "$LAB_HOST:$LAB_GTF" "$OUT/"

echo "[01] translating all BRAKER transcripts on lab (gffread -y to stdout)"
# -S: '*' for in-frame stops (gffread drops the terminal stop itself); no -V/-J, so no model is discarded.
"${SSH[@]}" "$LAB_GFFREAD -g $LAB_FASTA -S -y /dev/stdout -o /dev/null $LAB_GTF" \
  > "$OUT/braker_all_proteins.faa" 2> "$OUT/gffread.stderr"
n_prot=$(grep -c '^>' "$OUT/braker_all_proteins.faa")
echo "[01] proteins translated: $n_prot"

echo "[01] recording md5 sums and tool versions on lab"
"${SSH[@]}" "md5sum $LAB_RESULTS/GEA/tables/genes_per_region.tsv $LAB_RESULTS/GEA/tables/regions_combined.tsv $LAB_GTF $LAB_FASTA ${LAB_FASTA}.fai" \
  > "$OUT/lab_md5.txt"
{
  printf 'gffread\t%s\n' "$("${SSH[@]}" "$LAB_GFFREAD --version" 2>&1 | head -1)"
  printf 'samtools\t%s\n' "$("${SSH[@]}" "$LAB_SAMTOOLS --version" | head -1)"
  printf 'lab_host\t%s\n' "$("${SSH[@]}" hostname)"
  printf 'fetched\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/lab_tool_versions.txt"
cat "$OUT/lab_tool_versions.txt"
echo "[01] done"
