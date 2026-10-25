#!/usr/bin/env python3
"""Checks REGIONS_v2.md against the tables (plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7). COPY of
ANALYSIS/Trifolium86MainK5/regionsER/report_check.py (PC3 run) with two checks added, because the carried PC3 text would
otherwise pass with stale numbers:
  - p values: inside a region section every "p <number>" (and every number in a key-table Best-p cell) equals, at its
    printed precision, a v2 value of that region: a significant SNP's EMMAX min p or RDA p, or a region gene's best
    body / promoter p (EMMAX or RDA); the method thresholds are also allowed;
  - stale tokens: the PC3 set descriptors ("3 PCs", "3 conditioning PCs", "6.62e-6", "13 regions", "27 agreed",
    "315 genes", "28 significant SNPs") appear nowhere except on a line that says "PC3".

A region section starts at a `## N. <region_id>` heading and ends at the next `## ` heading. Within a section:
  - every gene id (gNNN) is a gene of that region (region_gene_evidence.tsv) or a +-500 kb neighbourhood gene of it
    (neighbourhood.tsv);
  - gene ids in a key-table row and in the Results paragraph are region genes, never neighbourhood genes;
  - the key-table distance cell: every number in it equals dist_to_sig_snp of one of the row's genes;
  - the Results paragraph has 120-200 words.
Everywhere:
  - every SNP id (chr:pos, pos >= 10,000) is a significant SNP of some region (region_sig_snps.tsv), or a tested SNP
    (EMMAX p-value table), or the line says it is "no longer tested";
  - every PMID (7-8 digits after "PMID", or in a key-table Literature cell) is a non-rejected row of literature.tsv;
  - every "(Author Year, PMID n)" matches that row's first author and year.
Exit 1 on any failure.
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/report_check.py   (from the repo root)
"""
import collections
import csv
import os
import re
import sys

B = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(B, "tables")
EMMAX = os.path.join(os.getcwd(), "Trifolium86MainK5ERpc5_results/GEA/tables/methods/EMMAX/EMMAX_pvalues_K5.tsv")


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


ev = rows(os.path.join(T, "region_gene_evidence.tsv"))
reg_genes = collections.defaultdict(set)
dist = {}
for r in ev:
    g = r["gene_id"].replace(".t1", "")
    reg_genes[r["region_id"]].add(g)
    dist[(r["region_id"], g)] = int(r["dist_to_sig_snp"])
nbh = collections.defaultdict(set)
for r in rows(os.path.join(T, "neighbourhood.tsv")):
    for rid in r["region_ids"].split(","):
        nbh[rid].add(r["gene_id"].replace(".t1", ""))
sig = {r["snp_id"] for r in rows(os.path.join(T, "region_sig_snps.tsv"))}
tested = set()
with open(EMMAX) as fh:
    next(fh)
    for line in fh:
        tested.add(line.split("\t", 1)[0])
# v2 p values per region (sig SNPs + region genes) and the thresholds
pv = collections.defaultdict(set)
for r in rows(os.path.join(T, "region_sig_snps.tsv")):
    pv[r["region_id"]] |= {float(r["EMMAX_minp"]), float(r["RDA_p"])}
for r in ev:
    for k in ("best_p_EMMAX_body", "best_p_RDA_body", "best_p_EMMAX_promoter", "best_p_RDA_promoter"):
        if r[k]:
            pv[r["region_id"]].add(float(r[k]))
THR = set()
for f in ("EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv", "RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv"):
    THR |= {float(r["pval_threshold"]) for r in rows(os.path.join(os.path.dirname(EMMAX), "..", f))}
NUM = r"(\d+(?:\.\d+)?(?:e-\d+)?)"
STALE = ["3 PCs", "3 conditioning PCs", "6.62e-6", "13 regions", "27 agreed", "315 genes", "28 significant SNPs"]


def same(txt, v):
    """txt printed at its precision equals v"""
    if "e" in txt:
        m, e = txt.split("e")
        nd = len(m.replace(".", "").lstrip("0"))
    else:
        nd = len(txt.replace(".", "").lstrip("0"))
    return float(f"{v:.{max(nd, 1)}g}") == float(txt)


lit = {}
for r in rows(os.path.join(T, "literature.tsv")):
    if r["role"] != "rejected":
        lit.setdefault(r["pmid"], set()).add((r["first_author"].split()[0], r["year"]))

fail = []
text = open(os.path.join(B, "REGIONS_v2.md")).read().split("\n")
cur, in_res = None, False
n_res = 0
for i, line in enumerate(text, 1):
    if line.startswith("## "):
        m = re.match(r"## \d+\. (\S+)", line)
        cur = m.group(1) if m and m.group(1) in reg_genes else None
        if m and not cur:
            fail.append(f"{i}: section heading names no region: {line}")
    is_res = line.startswith("**Results paragraph.**")
    is_key = line.startswith("| g")
    # genes
    for g in set(re.findall(r"\bg\d{3,6}\b", line)):
        if cur is None:
            continue
        if g in reg_genes[cur]:
            continue
        if g in nbh[cur] and not (is_res or is_key):
            continue
        fail.append(f"{i}: {g} is not a {'region' if (is_res or is_key) else 'region / neighbourhood'} gene of {cur}")
    # key-table distances
    if is_key and cur:
        cells = [c.strip() for c in line.strip("|").split("|")]
        genes = re.findall(r"\bg\d{3,6}\b", cells[0])
        ds = {dist[(cur, g)] for g in genes if (cur, g) in dist}
        for n in re.findall(r"\d[\d,]*", cells[3]):
            if int(n.replace(",", "")) not in ds:
                fail.append(f"{i}: distance {n} not among {sorted(ds)} for {genes}")
        for p in re.findall(r"\b\d{7,8}\b", cells[-1]):
            if p not in lit:
                fail.append(f"{i}: key-table PMID {p} not a verified, non-rejected literature row")
    # p values against the v2 tables
    if cur:
        nums = re.findall(r"\bp " + NUM, line)
        if is_key:
            cells = [c.strip() for c in line.strip("|").split("|")]
            nums += re.findall(r"(?<![\w.])(\d+\.\d+(?:e-\d+)?|\d(?:\.\d+)?e-\d+)(?![\w.])", cells[5])
        for n in nums:
            if not any(same(n, v) for v in pv[cur] | THR):
                fail.append(f"{i}: p {n} is not a v2 value of {cur}")
    for t in STALE:
        if t in line and "PC3" not in line:
            fail.append(f"{i}: stale PC3 descriptor {t!r}")
    # results paragraph length
    if is_res:
        n_res += 1
        words = len(line.replace("**Results paragraph.**", "").split())
        if not 120 <= words <= 200:
            fail.append(f"{i}: Results paragraph of {cur} has {words} words (120-200)")
    # SNPs
    for s in set(re.findall(r"\b([1-6]):(\d{5,})\b", line)):
        sid = f"{s[0]}:{s[1]}"
        if sid in sig or sid in tested or "no longer tested" in line:
            continue
        fail.append(f"{i}: SNP {sid} is neither significant nor tested")
    # PMIDs and Author-Year pairs
    for p in re.findall(r"PMID (\d{7,8})", line):
        if p not in lit:
            fail.append(f"{i}: PMID {p} not a verified, non-rejected literature row")
    for au, yr, p in re.findall(r"\(([A-Z][\w\-ü óé]+?) (\d{4}), PMID (\d{7,8})\)", line):
        if p in lit and not any(a == au.split()[0] and y == yr for a, y in lit[p]):
            fail.append(f"{i}: ({au} {yr}, PMID {p}) does not match {sorted(lit[p])}")

if n_res != 15:
    fail.append(f"expected 15 Results paragraphs, found {n_res}")
for f in fail:
    print("FAIL", f)
print(f"INFO: {len(text)} lines, {n_res} Results paragraphs, {len(fail)} failures")
sys.exit(1 if fail else 0)
