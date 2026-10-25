#!/usr/bin/env bash
# 04_blastp.sh -- BLASTP of the 353 Trifolium proteins against the four UniProt sets.
#
# BLAST+ comes from a pinned nixpkgs revision (config.sh: NIXPKGS_REF). Databases are built from copies of
# the UniProt FASTA files whose headers are cut to the bare accession, so sseqid is the UniProt accession.
# Default blastp scoring (BLOSUM62, gap open 11 / extend 1, composition-based statistics 2, no SEG);
# E-value cutoff and max_target_seqs from config.sh. max_target_seqs is large on purpose: it acts inside the
# search, so the top-N cut is made afterwards by 05_summarise.py.
#
# Usage: bash 04_blastp.sh <data_dir>
# Writes: <data_dir>/blastdb/<set>.*, <data_dir>/blast/<set>.blastp.tsv, <data_dir>/blast/blast_version.txt
set -euo pipefail

DATA=${1:?usage: 04_blastp.sh <data_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"

DB="$DATA/blastdb"; OUT="$DATA/blast"
mkdir -p "$DB" "$OUT"
QUERY="$DATA/work/queries.faa"
[[ -s "$QUERY" ]] || { echo "[04] missing $QUERY (run 02 first)" >&2; exit 1; }

nixrun() { nix shell "${NIXPKGS_REF}#blast" -c "$@"; }
nixrun blastp -version > "$OUT/blast_version.txt"
cat "$OUT/blast_version.txt"

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp"
printf '%s\n' "$OUTFMT" > "$OUT/outfmt.txt"

for set in at_proteome at_swissprot mt_proteome gm_proteome; do
  src="$DATA/uniprot/$set.fasta"
  [[ -s "$src" ]] || { echo "[04] missing $src (run 03 first)" >&2; exit 1; }
  # >sp|P69781|GAT12_ARATH ... -> >P69781
  awk '/^>/ { split(substr($1, 2), a, "|"); print ">" a[2]; next } { print }' "$src" > "$DB/$set.faa"
  nixrun makeblastdb -in "$DB/$set.faa" -dbtype prot -out "$DB/$set" -title "$set" > "$DB/$set.makeblastdb.log"
  echo "[04] blastp vs $set"
  nixrun blastp -query "$QUERY" -db "$DB/$set" -evalue "$EVALUE" -max_target_seqs "$MAX_TARGET_SEQS" \
    -num_threads "$THREADS" -outfmt "$OUTFMT" > "$OUT/$set.blastp.tsv"
  printf '[04] %s: %s HSP rows, %s queries with hits\n' "$set" "$(wc -l < "$OUT/$set.blastp.tsv")" \
    "$(cut -f1 "$OUT/$set.blastp.tsv" | sort -u | wc -l)"
done
echo "[04] done"
