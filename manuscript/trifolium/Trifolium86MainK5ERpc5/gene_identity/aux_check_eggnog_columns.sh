#!/usr/bin/env bash
# aux_check_eggnog_columns.sh -- side check, not part of the identity procedure: compares the eggNOG columns of
# the pipeline table genes_per_region.tsv with the attributes of the same transcripts in the eggNOG GFF on lab
# (Trifolium_eggnog.gff, read with grep on lab and printed to stdout).
#
# Usage: bash aux_check_eggnog_columns.sh <data_dir>     (after 01 and 02)
# Writes: <data_dir>/aux_eggnog/{ids.txt, eggnog_subset.gff, column_check.tsv}
set -euo pipefail

DATA=${1:?usage: aux_check_eggnog_columns.sh <data_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"
OUT="$DATA/aux_eggnog"
mkdir -p "$OUT"
awk -F'\t' 'NR > 1 {print "ID=" $3 ";"}' "$DATA/work/genes.tsv" > "$OUT/ids.txt"
ssh -o BatchMode=yes "$LAB_HOST" "grep -F -f - ${LAB_GENOME_DIR}/Trifolium_eggnog.gff" < "$OUT/ids.txt" > "$OUT/eggnog_subset.gff"

python3 -I - "$OUT/eggnog_subset.gff" "$DATA/lab/genes_per_region.tsv" "$OUT/column_check.tsv" <<'EOF'
import csv, sys
gff, table, out = sys.argv[1:4]
eg = {}
for line in open(gff):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 9:
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        eg[a.get("ID")] = a
rows = {r["gene_id"]: r for r in csv.DictReader(open(table), delimiter="\t")}
cols = ["em_target", "em_desc", "em_COG_cat", "em_OGs", "em_EC", "em_KEGG_ko", "em_BRITE", "em_PFAMs"]
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["column", "genes_compared", "genes_differing", "example_gene", "pipeline_value", "eggnog_gff_value"])
    for c in cols:
        diff = [(g, r[c], eg[g].get(c, "")) for g, r in sorted(rows.items()) if g in eg and r[c] != eg[g].get(c, "")]
        ex = diff[0] if diff else ("", "", "")
        w.writerow([c, sum(1 for g in rows if g in eg), len(diff), ex[0], ex[1][:60], ex[2][:60]])
        print(f"[aux] {c}: {len(diff)} of {sum(1 for g in rows if g in eg)} genes differ")
EOF
