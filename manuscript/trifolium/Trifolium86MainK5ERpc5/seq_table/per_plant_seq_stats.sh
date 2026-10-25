#!/usr/bin/env bash
# Per-plant sequencing statistics of the Trifolium israeliticum plants that were called
# (Methods, "Plant material, sequencing and variant calling"; the per-plant sheet of Table S3).
# Written 2026-10-25 (user: reads, depth and missing calls per plant, so "low-coverage" has a number).
#
# Input: one BAM per called plant (bwa-mem2 alignment, Picard FixMateInformation, MarkDuplicates with
# duplicates marked, not removed, and AddOrReplaceReadGroups; trim_align.sh and process.sh on the lab
# server), and the framework's per-sample missingness table of the run Trifolium86MainK5ERpc5.
# For each BAM, samtools (1.17 on the lab server) gives:
#   read_pairs     flagstat "read1" (QC-passed), i.e. read pairs after adapter trimming
#   mapped_pct     flagstat "primary mapped %"
#   duplicate_pct  flagstat "primary duplicates" / "primary" x 100
#   mean_depth     samtools coverage over the six pseudochromosomes Chr_1-Chr_6, the length-weighted
#                  mean of "meandepth" (samtools coverage's default read filter: unmapped, secondary,
#                  QC-failed and duplicate reads excluded; no mapping- or base-quality threshold)
# Missing calls: F_MISS from plink --missing over the 7,163,793 SNPs the framework read
#   (<run>_results/_intermediate/samples/samples_missing_stats.tsv).
#
# Usage (on the lab server): per_plant_seq_stats.sh <bam_dir> <missing_stats.tsv> <out_dir> [jobs]
# Writes <out_dir>/per_plant_seq_stats.tsv and <out_dir>/per_bam/<plant>.{flagstat.tsv,coverage.tsv}.
set -euo pipefail

BAM_DIR=${1:?bam_dir}
MISS=${2:?missing_stats.tsv}
OUT=${3:?out_dir}
JOBS=${4:-12}
SAMTOOLS=$(command -v samtools) || { echo "samtools not found" >&2; exit 1; }
mkdir -p "$OUT/per_bam"
"$SAMTOOLS" --version | head -2 > "$OUT/samtools_version.txt"

one_bam() {
    bam=$1; out=$2; st=$3
    plant=$(basename "$bam" | grep -oE '^[0-9]{5}-[0-9]+')
    [ -s "$out/per_bam/$plant.coverage.tsv" ] && [ -s "$out/per_bam/$plant.flagstat.tsv" ] && return 0
    "$st" flagstat -@ 2 -O tsv "$bam" > "$out/per_bam/$plant.flagstat.tsv.tmp"
    "$st" coverage "$bam" > "$out/per_bam/$plant.coverage.tsv.tmp"
    mv "$out/per_bam/$plant.flagstat.tsv.tmp" "$out/per_bam/$plant.flagstat.tsv"
    mv "$out/per_bam/$plant.coverage.tsv.tmp" "$out/per_bam/$plant.coverage.tsv"
}
export -f one_bam

ls "$BAM_DIR"/*.sorted.fixmate.deDup.RG.bam > "$OUT/bams.list"
n=$(wc -l < "$OUT/bams.list")
[ "$n" -gt 0 ] || { echo "no BAMs in $BAM_DIR" >&2; exit 1; }
echo "$(date) start: $n BAMs, $JOBS at a time" >&2
xargs -a "$OUT/bams.list" -P "$JOBS" -I{} bash -c 'one_bam "$@"' _ {} "$OUT" "$SAMTOOLS"
echo "$(date) samtools done" >&2

# Assemble the table (awk only, so the script needs nothing beyond samtools and coreutils).
{
    printf "plant\tsite\tread_pairs\tmapped_pct\tduplicate_pct\tmean_depth_chr1_6\tmissing_calls_pct\n"
    for f in "$OUT"/per_bam/*.flagstat.tsv; do
        plant=$(basename "$f" .flagstat.tsv)
        site=${plant%%-*}
        pairs=$(awk -F'\t' '$3 == "read1" {print $1}' "$f")
        mapped=$(awk -F'\t' '$3 == "primary mapped %" {sub(/%$/, "", $1); print $1}' "$f")
        dup=$(awk -F'\t' '$3 == "primary" {p = $1} $3 == "primary duplicates" {d = $1} END {printf "%.2f", 100 * d / p}' "$f")
        depth=$(awk -F'\t' '$1 ~ /^Chr_[1-6]$/ {len = $3 - $2 + 1; s += $7 * len; L += len} END {printf "%.2f", s / L}' "$OUT/per_bam/$plant.coverage.tsv")
        miss=$(awk -v p="$plant" 'NR > 1 && $2 == p {printf "%.2f", 100 * $6}' "$MISS")
        [ -n "$miss" ] || { echo "no missingness row for $plant" >&2; exit 1; }
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$plant" "$site" "$pairs" "$mapped" "$dup" "$depth" "$miss"
    done
} > "$OUT/per_plant_seq_stats.tsv"
rows=$(($(wc -l < "$OUT/per_plant_seq_stats.tsv") - 1))
[ "$rows" -eq "$n" ] || { echo "table has $rows rows for $n BAMs" >&2; exit 1; }
echo "$(date) wrote $OUT/per_plant_seq_stats.tsv ($rows plants)" >&2
