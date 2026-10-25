#!/usr/bin/env bash
# 03_fetch_uniprot.sh -- download the four UniProt protein sets (FASTA) and their annotation tables (TSV)
# through the UniProt REST stream endpoint, recording the UniProt release from each response's headers.
#
# Sets:
#   at_proteome   Arabidopsis thaliana reference proteome     proteome:UP000006548
#   at_swissprot  Arabidopsis thaliana Swiss-Prot             organism_id:3702 AND reviewed:true
#   mt_proteome   Medicago truncatula reference proteome      proteome:UP000002051
#   gm_proteome   Glycine max reference proteome              proteome:UP000008827
# Annotation fields: see UNIPROT_FIELDS in config.sh (accession, gene names, ordered locus name = AGI,
# protein name, length, PANTHER cross-references, ...).
#
# Fails when a download's record count differs from the search endpoint's x-total-results, or when the
# downloads come from different UniProt releases.
#
# Usage: bash 03_fetch_uniprot.sh <data_dir>
# Writes: <data_dir>/uniprot/{<set>.fasta, <set>.annot.tsv, *.headers, downloads.tsv}
set -euo pipefail

DATA=${1:?usage: 03_fetch_uniprot.sh <data_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"

OUT="$DATA/uniprot"          # downloaded data only (untrusted; nothing is executed from here)
mkdir -p "$OUT"

declare -A QUERY=(
  [at_proteome]="proteome:UP000006548"
  [at_swissprot]="organism_id:3702 AND reviewed:true"
  [mt_proteome]="proteome:UP000002051"
  [gm_proteome]="proteome:UP000008827"
)
SETS=(at_proteome at_swissprot mt_proteome gm_proteome)

urlenc() { python3 -I -c 'import sys, urllib.parse; print(urllib.parse.quote(sys.argv[1], safe=""))' "$1"; }
header() { tr -d '\r' < "$1" | awk -F': ' -v k="$2" 'tolower($1)==k {print $2}' | tail -1; }

manifest="$OUT/downloads.tsv"
printf 'set\tquery\tformat\tfile\tuniprot_release\trelease_date\tapi_deployment_date\tx_total_results\tn_records\tmd5\tdownloaded_utc\n' > "$manifest"

for set in "${SETS[@]}"; do
  q=$(urlenc "${QUERY[$set]}")
  curl -sS --fail --retry 3 -D "$OUT/$set.search.headers" -o /dev/null \
    "https://rest.uniprot.org/uniprotkb/search?query=${q}&size=1&format=fasta"
  total=$(header "$OUT/$set.search.headers" x-total-results)
  for fmt in fasta tsv; do
    if [[ $fmt == fasta ]]; then
      file="$OUT/$set.fasta"
      url="${UNIPROT_STREAM}?format=fasta&compressed=true&query=${q}"
    else
      file="$OUT/$set.annot.tsv"
      url="${UNIPROT_STREAM}?format=tsv&compressed=true&fields=${UNIPROT_FIELDS}&query=${q}"
    fi
    echo "[03] $set $fmt"
    curl -sS --fail --retry 3 -D "$file.headers" -o "$file.gz" "$url"
    gunzip -f "$file.gz"
    if [[ $fmt == fasta ]]; then n=$(grep -c '^>' "$file"); else n=$(( $(wc -l < "$file") - 1 )); fi
    if [[ "$n" != "$total" ]]; then
      echo "[03] ERROR: $set $fmt has $n records, search reports $total" >&2
      exit 1
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$set" "${QUERY[$set]}" "$fmt" "$(basename "$file")" \
      "$(header "$file.headers" x-uniprot-release)" "$(header "$file.headers" x-uniprot-release-date)" \
      "$(header "$file.headers" x-api-deployment-date)" "$total" "$n" "$(md5sum "$file" | cut -d' ' -f1)" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$manifest"
  done
done

n_rel=$(tail -n +2 "$manifest" | cut -f5 | sort -u | wc -l)
if [[ "$n_rel" != 1 ]]; then
  echo "[03] ERROR: downloads span more than one UniProt release:" >&2
  tail -n +2 "$manifest" | cut -f1,3,5 >&2
  exit 1
fi
column -t -s $'\t' "$manifest"
echo "[03] done"
