#!/usr/bin/env bash
# 09_tblastn.sh -- TBLASTN check of possibly truncated gene models among the named genes.
#
# Trigger: a named gene (named_genes_spec.tsv) whose best Arabidopsis hit (UP000006548, 07) covers less than
# TBLASTN_SCOV_MAX of the Arabidopsis protein (subject cover = merged HSP length / protein length).
# For each such gene, the Arabidopsis protein (representative entry of the best locus) is searched with
# TBLASTN (default parameters, E-value cutoff from config.sh) against the gene's region +/- WINDOW_FLANK bp,
# extracted on lab with samtools faidx (printed to stdout; nothing is written on lab).
# 10_named_genes_check.py reads the HSP table and reports whether HSPs fall outside the gene model.
#
# Usage: bash 09_tblastn.sh <data_dir>
# Writes: <data_dir>/tblastn/{triggers.tsv, <gene>.query.faa, <region>.window.fa, <gene>.tblastn.tsv,
#         window_coords.tsv}
set -euo pipefail

DATA=${1:?usage: 09_tblastn.sh <data_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"

OUT="$DATA/tblastn"
mkdir -p "$OUT"
nixrun() { nix shell "${NIXPKGS_REF}#blast" -c "$@"; }

# 1. triggers: named genes whose best Arabidopsis hit has subject cover < TBLASTN_SCOV_MAX
python3 -I - "$DATA" "$HERE/named_genes_spec.tsv" "$TBLASTN_SCOV_MAX" "$WINDOW_FLANK" > "$OUT/triggers.tsv" <<'EOF'
import csv, sys
data, spec, smax, flank = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
named = sorted({r["gene"] for r in csv.DictReader(open(spec), delimiter="\t")})
genes = {r["transcript_id"]: r for r in csv.DictReader(open(f"{data}/work/genes.tsv"), delimiter="\t")}
best = {}
for r in csv.DictReader(open(f"{data}/work/arabidopsis_loci.tsv"), delimiter="\t"):
    if r["transcript_id"] not in best:          # file is ordered by bitscore within gene
        best[r["transcript_id"]] = r
print("\t".join(["gene", "region_id", "best_locus", "accession", "scov_pct", "trigger", "window"]))
for g in named:
    b = best.get(g)
    if b is None:
        continue
    scov = float(b["scov"])
    reg = genes[g]["region_id"]
    chrom, span = reg.split("_", 1)
    s, e = (int(x) for x in span.split("-"))
    win = f"Chr_{chrom}:{max(1, s - flank)}-{e + flank}"
    print("\t".join([g, reg, b["locus"], b["accession"], f"{scov:.1f}", "yes" if scov < 100 * smax else "no", win]))
EOF
column -t -s $'\t' "$OUT/triggers.tsv"

# 2. per triggered gene: window (once per region), query protein, tblastn
: > "$OUT/window_coords.tsv"
tail -n +2 "$OUT/triggers.tsv" | awk -F'\t' '$6=="yes"' | while IFS=$'\t' read -r gene region locus acc scov trig win; do
  wfa="$OUT/${region}.window.fa"
  if [[ ! -s "$wfa" ]]; then
    ssh -o BatchMode=yes "$LAB_HOST" "$LAB_SAMTOOLS faidx $LAB_FASTA $win" < /dev/null > "$wfa"
    nixrun makeblastdb -in "$wfa" -dbtype nucl -out "$OUT/${region}.window" > "$OUT/${region}.makeblastdb.log"
  fi
  grep -q -P "^${region}\t" "$OUT/window_coords.tsv" || printf '%s\t%s\n' "$region" "$win" >> "$OUT/window_coords.tsv"
  awk -v a="$acc" '/^>/{p=($0 ~ "^>" a "$")} p' "$DATA/blastdb/at_proteome.faa" > "$OUT/${gene}.query.faa"
  [[ -s "$OUT/${gene}.query.faa" ]] || { echo "[09] no sequence for $acc" >&2; exit 1; }
  nixrun tblastn -query "$OUT/${gene}.query.faa" -db "$OUT/${region}.window" -evalue "$EVALUE" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    > "$OUT/${gene}.tblastn.tsv"
  echo "[09] $gene ($acc vs $win): $(wc -l < "$OUT/${gene}.tblastn.tsv") HSPs"
done
echo "[09] done"
