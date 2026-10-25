#!/usr/bin/env bash
# 11_versions.sh -- collect tool versions, database releases, entry counts, web-service records and parameters
# of a run into versions.tsv.
#
# Usage: bash 11_versions.sh <data_dir> <out_dir>
# Writes: <out_dir>/versions.tsv (columns: component, item, value, source)
set -euo pipefail

DATA=${1:?usage: 11_versions.sh <data_dir> <out_dir>}
OUTDIR=${2:?usage: 11_versions.sh <data_dir> <out_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=config.sh
source "$HERE/config.sh"

BLAST_STORE=$(nix eval --raw "${NIXPKGS_REF}#blast.outPath")
PYVER=$(python3 -I -c 'import sys; print(sys.version.split()[0])')

python3 -I - "$DATA" "$OUTDIR/versions.tsv" "$BLAST_STORE" "$NIXPKGS_REF" "$PYVER" "$EVALUE" "$MAX_TARGET_SEQS" \
  "$TOP_N" "$SYN_BITSCORE_FRAC" "$SYN_AGI_WINDOW" "$SYN_MIN_SUPPORT" "$TBLASTN_SCOV_MAX" "$WINDOW_FLANK" \
  "$CDSEARCH_POLL_S" <<'EOF'
import csv, datetime, os, sys
(data, out, blast_store, nixpkgs, pyver, evalue, mts, top_n, frac, win, sup, scov, flank, poll) = sys.argv[1:15]
rows = []
def add(component, item, value, source):
    rows.append([component, item, str(value), source])
def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

today = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")
add("run", "date (UTC) of this versions table", today, "11_versions.sh")
# BLAST+
ver = open(os.path.join(data, "blast", "blast_version.txt")).read().split("\n")
add("BLAST+", "blastp / makeblastdb / tblastn version", ver[0].strip(), "blastp -version")
add("BLAST+", "build", ver[1].strip() if len(ver) > 1 else "", "blastp -version")
add("BLAST+", "nix store path", blast_store, f"nix eval {nixpkgs}#blast.outPath")
add("BLAST+", "nixpkgs revision", nixpkgs, "config.sh NIXPKGS_REF")
add("BLASTP", "E-value cutoff", evalue, "config.sh EVALUE")
add("BLASTP", "max_target_seqs (inside the search)", mts, "config.sh MAX_TARGET_SEQS")
add("BLASTP", "hits kept per gene and database", top_n, "config.sh TOP_N (cut after collapsing to loci)")
add("BLASTP", "scoring", "defaults: BLOSUM62, gap open 11 / extend 1, composition-based statistics 2 "
    "(compositional matrix adjust), SEG off", "blastp -help and the pairwise report footer (checked 2026-10-25)")
add("TBLASTN", "trigger", f"named gene whose best Arabidopsis hit covers < {float(scov) * 100:.0f}% of the Arabidopsis protein",
    "config.sh TBLASTN_SCOV_MAX")
add("TBLASTN", "E-value cutoff / window", f"{evalue}; region +/- {flank} bp", "config.sh EVALUE, WINDOW_FLANK")
add("TBLASTN", "scoring", "defaults: BLOSUM62, composition-based statistics 2, SEG on the query (12 2.2 2.5)",
    "tblastn -help (checked 2026-10-25)")
# lab tools and inputs
for line in open(os.path.join(data, "lab", "lab_tool_versions.txt")):
    k, v = line.rstrip("\n").split("\t", 1)
    add("lab", k, v, "01_fetch_lab_inputs.sh")
for line in open(os.path.join(data, "lab", "lab_md5.txt")):
    md5, path = line.split()
    add("lab input md5", os.path.basename(path), md5, path)
# UniProt
for r in read_tsv(os.path.join(data, "uniprot", "downloads.tsv")):
    if r["format"] != "fasta":
        continue
    add("UniProt", f"{r['set']} ({r['query']})",
        f"release {r['uniprot_release']} ({r['release_date']}); {r['n_records']} entries; fasta md5 {r['md5']}; "
        f"downloaded {r['downloaded_utc']}",
        "X-UniProt-Release / X-UniProt-Release-Date headers; x-total-results; 03_fetch_uniprot.sh")
annot = [r for r in read_tsv(os.path.join(data, "uniprot", "downloads.tsv")) if r["format"] == "tsv"]
add("UniProt", "annotation tables (TSV)", "; ".join(f"{r['set']} {r['n_records']} rows, release {r['uniprot_release']}"
                                                  for r in annot), "03_fetch_uniprot.sh (fields in config.sh)")
# CD-Search
meta = {r["key"]: r["value"] for r in read_tsv(os.path.join(data, "cdsearch", "cdsearch_meta.tsv"))}
add("CD-Search", "service", meta.get("service", ""), "05_cdsearch.py")
add("CD-Search", "search id (cdsid)", meta.get("cdsid", ""), "service response")
add("CD-Search", "submitted / finished (UTC)", f"{meta.get('submitted_utc')} / {meta.get('finished_utc')}", "05_cdsearch.py")
add("CD-Search", "CDD version", meta.get("cdd_version_in_service_output", ""), "service output")
add("CD-Search", "NCBI FTP cdd.info (public release, for context only)", meta.get("ncbi_ftp_cdd_info_first_line", ""),
    "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.info")
add("CD-Search", "parameters", ", ".join(f"{k[6:]}={v}" for k, v in meta.items() if k.startswith("param_")),
    "05_cdsearch.py")
add("CD-Search", "queries / with any hit / with a specific hit",
    f"{meta.get('n_queries')} / {meta.get('n_queries_with_any_hit')} / {meta.get('n_queries_with_specific_hit')}",
    "cdsearch_hits.tsv")
add("CD-Search", "service warnings", meta.get("service_warnings", ""), "service output")
add("CD-Search", "poll interval (s)", poll, "config.sh CDSEARCH_POLL_S")
# family name sources
for r in read_tsv(os.path.join(data, "work", "family_sources.tsv")):
    add("Medicago family names", r["family"], f"{r['citation']}; doi {r['doi']}; PMID {r['pmid']}; {r['pmcid']}"
        if r["doi"] else r["citation"], f"{r['fetched_from']} ({r['fetched_utc']})")
add("Medicago family names", "Li 2016 BURP subfamilies not stated in the text",
    "read by eye from Fig. 2 (tree arcs) and Fig. 1 (label colours) into li2016_burp_figure_classes.tsv; "
    "checked against the 23 members whose class the text states (all agree)",
    "manual figure reading, 2026-10-25")
add("Medicago family names", "Li 2016 figure source",
    "publisher PDF, paper-shelf copy mac-studio:~/Papers/10.3389_fpls.2016.00485.pdf (shelf source: unpaywall, "
    "https://www.frontiersin.org/articles/10.3389/fpls.2016.00485/pdf; md5 e15cad5c90a83bd9781f76ec130533ba); "
    "Fig. 1 and Fig. 2 images extracted from PDF pages 5-6 with pdfimages (poppler-utils 26.10.0, nixpkgs)",
    "manual figure reading, 2026-10-25")
# synteny and ties
add("microsynteny", "candidate loci", f"bitscore >= {frac} x best Arabidopsis locus; nuclear AGI codes only", "config.sh")
add("microsynteny", "link", f"different locus, same Arabidopsis chromosome, |AGI number difference| <= {win}", "config.sh")
add("microsynteny", "syntenic call", f">= {sup} other genes of the region linked; tie-break: most supporting genes, "
    "then bitscore", "config.sh")
add("near-tie", "definition", f"bitscore >= {frac} x the best locus bitscore", "config.sh SYN_BITSCORE_FRAC")
add("python", "version", pyver, "python3 -I (stdlib only)")
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["component", "item", "value", "source"])
    w.writerows(rows)
print(f"[11] {len(rows)} rows -> {out}")
EOF
