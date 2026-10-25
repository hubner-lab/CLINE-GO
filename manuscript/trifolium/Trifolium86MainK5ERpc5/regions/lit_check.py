#!/usr/bin/env python3
"""Verify every literature row against NCBI PubMed.

COPY (2026-10-16, plan ~/.claude/plans/peaceful-kindling-micali.md Phase 7) of ANALYSIS/Trifolium86MainK5/regionsER/
lit_check.py, unchanged except that cached abstracts are also copied from the PC3 run's regionsER/logs/abstracts/.

COPY (2026-10-14, plan ~/.claude/plans/wild-cooking-kite.md Phase 3) of ANALYSIS/Trifolium86Main/deepdive6/lit_check.py:
identical checks; the directory is argv[1] (reads <dir>/tables/literature_raw.tsv, writes <dir>/tables/literature.tsv,
<dir>/logs/lit_dropped.tsv, <dir>/logs/abstracts/), and the output keeps the extra columns direction / role.
Abstracts already cached in ANALYSIS/Trifolium86Main/deepdive6/logs/abstracts/ are copied, not re-fetched.
Usage: python3 -I lit_check.py <dir>

Original header:

Input : tables/literature_raw.tsv, one row per (region_id, gene_id, pmid) as returned by the literature agents
        (columns region_id gene_id role pmid first_author year species legume evidence finding). role = key | btw |
        seed (seed = the 7 Phase 1 R3 PMIDs of genes3of3, verified here for the first time).
Check : E-utilities esummary per PMID. A row is KEPT when the PMID resolves, the claimed first author's surname
        matches PubMed's first author, and the claimed year equals the publication year. Otherwise it is DROPPED and
        listed with the reason. If NCBI cannot be reached the script stops (no fallback).
Output: tables/literature.tsv (kept rows + PubMed title, journal, first author as indexed),
        logs/lit_dropped.tsv, logs/abstracts/<pmid>.txt (efetch abstract, read before a paper is cited).
Run: python3 -I ANALYSIS/Trifolium86Main/deepdive6/lit_check.py   (from the repo root)
"""
import csv
import json
import os
import sys
import time
import unicodedata
import urllib.request

if len(sys.argv) != 2:
    sys.exit("usage: lit_check.py <dir>")
B = os.path.abspath(sys.argv[1])
HERE = os.path.dirname(os.path.abspath(__file__))
OLD_ABS_DIRS = [os.path.join(HERE, "../../Trifolium86MainK5/regionsER/logs/abstracts"),
                os.path.join(HERE, "../../Trifolium86Main/deepdive6/logs/abstracts")]
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read().decode()
        except Exception as e:  # network failure is fatal after retries
            last = e
            time.sleep(2 + 2 * i)
    sys.exit(f"ERROR: NCBI unreachable ({last}); stopping, no fallback")


def norm(s):
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()
    return "".join(ch for ch in s if ch.isalpha())


with open(os.path.join(B, "tables", "literature_raw.tsv")) as fh:
    raw = list(csv.DictReader(fh, delimiter="\t"))
pmids = sorted({r["pmid"] for r in raw if r["pmid"].isdigit()})
summ = {}
for i in range(0, len(pmids), 100):
    d = json.loads(get(EU + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(pmids[i:i + 100])))["result"]
    for u in d.get("uids", []):
        summ[u] = d[u]
    time.sleep(0.4)

kept, dropped = [], []
for r in raw:
    s = summ.get(r["pmid"])
    why = None
    if not r["pmid"].isdigit() or s is None or "error" in s:
        why = "PMID does not resolve"
    else:
        first = s.get("sortfirstauthor", "") or (s.get("authors") or [{"name": ""}])[0]["name"]
        claimed = norm(r["first_author"].split()[0]) if r["first_author"].strip() else ""
        if not claimed or not norm(first).startswith(claimed):
            why = f"first author mismatch: claimed {r['first_author']!r}, PubMed {first!r}"
        elif s.get("pubdate", "")[:4] != r["year"].strip():
            why = f"year mismatch: claimed {r['year']}, PubMed {s.get('pubdate', '')[:4]}"
    if why:
        dropped.append({**r, "reason": why})
    else:
        kept.append({**r, "pubmed_first_author": first, "title": s["title"], "journal": s.get("source", "")})

os.makedirs(os.path.join(B, "logs", "abstracts"), exist_ok=True)
for p in sorted({r["pmid"] for r in kept}):
    f = os.path.join(B, "logs", "abstracts", p + ".txt")
    for d in OLD_ABS_DIRS:
        old = os.path.join(d, p + ".txt")
        if (not os.path.exists(f) or os.path.getsize(f) == 0) and os.path.exists(old) and os.path.getsize(old) > 0:
            with open(old) as a, open(f, "w") as b:
                b.write(a.read())
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        txt = get(EU + f"efetch.fcgi?db=pubmed&id={p}&rettype=abstract&retmode=text")
        with open(f, "w") as fh:
            fh.write(txt)
        time.sleep(0.4)

cols = ["region_id", "gene_id", "role", "pmid", "first_author", "pubmed_first_author", "year", "species", "legume",
        "evidence", "direction", "title", "journal", "finding"]
with open(os.path.join(B, "tables", "literature.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(kept)
with open(os.path.join(B, "logs", "lit_dropped.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, list(raw[0].keys()) + ["reason"], delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(dropped)
print(f"INFO: {len(raw)} literature rows, {len(pmids)} distinct PMIDs; kept {len(kept)} rows "
      f"({len({r['pmid'] for r in kept})} PMIDs), dropped {len(dropped)}")
for r in dropped:
    print(f"  DROPPED {r['region_id']} {r['gene_id']} {r['pmid']}: {r['reason']}")
