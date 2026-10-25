#!/usr/bin/env bash
# run_all.sh -- run the gene-identity procedure end to end, in order.
#
# Usage: bash run_all.sh <data_dir> <curated_dir> [--cdsearch-reparse]
#   <data_dir>     large inputs/outputs (proteomes, BLAST databases, raw results); created if absent
#   <curated_dir>  holds the curated tables region_gene_evidence.tsv and evidence_rgt12.tsv
#   --cdsearch-reparse  do not submit a new NCBI CD-Search; re-parse <data_dir>/cdsearch/cdsearch_raw.txt
# The three result tables (gene_identity_all_genes.tsv, named_genes_check.tsv, versions.tsv) are written next to
# this script; the supplement's genes table (t_genes.tsv, t_genes.md) to docs/manuscript/figures/supp/.
set -euo pipefail

DATA=${1:?usage: run_all.sh <data_dir> <curated_dir> [--cdsearch-reparse]}
CUR=${2:?usage: run_all.sh <data_dir> <curated_dir> [--cdsearch-reparse]}
REPARSE=${3:-}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"
mkdir -p "$DATA"
EV="$CUR/region_gene_evidence.tsv"; RGT12="$CUR/evidence_rgt12.tsv"
[[ -s "$EV" && -s "$RGT12" ]] || { echo "missing curated tables in $CUR" >&2; exit 1; }

bash "$HERE/01_fetch_lab_inputs.sh" "$DATA"
python3 -I "$HERE/02_prepare_queries.py" "$DATA" "$EV" 353
bash "$HERE/03_fetch_uniprot.sh" "$DATA"
bash "$HERE/04_blastp.sh" "$DATA"
if [[ "$REPARSE" == "--cdsearch-reparse" ]]; then
  python3 -I "$HERE/05_cdsearch.py" "$DATA" --reparse
else
  python3 -I "$HERE/05_cdsearch.py" "$DATA" "$CDSEARCH_POLL_S" "$CDSEARCH_MAX_WAIT_S"
fi
python3 -I "$HERE/06_family_tables.py" "$DATA" "$HERE/li2016_burp_figure_classes.tsv"
python3 -I "$HERE/07_summarise.py" "$DATA" "$HERE" "$TOP_N" "$SYN_BITSCORE_FRAC" "$EV" "$RGT12"
python3 -I "$HERE/08_microsynteny.py" "$DATA" "$SYN_REGIONS" "$SYN_BITSCORE_FRAC" "$SYN_AGI_WINDOW" "$SYN_MIN_SUPPORT"
bash "$HERE/09_tblastn.sh" "$DATA"
python3 -I "$HERE/10_named_genes_check.py" "$DATA" "$HERE" "$HERE/named_genes_spec.tsv" "$EV" "$RGT12" "$SYN_BITSCORE_FRAC"
bash "$HERE/11_versions.sh" "$DATA" "$HERE"
bash "$HERE/12_fetch_associated_snps.sh" "$DATA"
python3 -I "$HERE/13_table_genes.py" "$DATA" "$HERE/gene_identity_all_genes.tsv" "$HERE/../../../figures/supp"
echo "[run_all] done"
