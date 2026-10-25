#!/usr/bin/env python3
"""Inputs for the region-QC report of the 15 Trifolium86MainK5ERpc5 (5 PCs) regions (plan
~/.claude/plans/peaceful-kindling-micali.md, Phase 7). COPY of ANALYSIS/Trifolium86MainK5/regionsER/prep_regionplots.py
(PC3 run, plan wild-cooking-kite.md): labels of the two dropped chr-6 regions removed, the four new regions and the new
sig-SNP host genes added; key tables parsed from REGIONS_v2.md.

Window = region +-100 kb (user decision 2026-10-15); the region edges are drawn as dotted vertical lines.
Labelled genes (all gene boxes one colour, every other gene an unnamed box):
  - key   : the key-table genes of REGIONS_v2.md, every tandem copy named;
  - host  : genes holding a significant SNP in their body or 5 kb promoter (tables/region_sig_snps.tsv gene_hit);
  - out   : priority-A neighbourhood guesses (tables/neighbourhood_guess.tsv) inside the window.
Priority-A guesses outside the window go to the gene table only.

Writes
  figures/input/topr_labelled.tsv     deepdive6's gff2topr.py output relabelled (the GFF is byte-identical, see plan)
  figures/input/regionplot_calls.tsv  order, region_id, window, traits argument (sig SNPs in the window), edges
  report/region_genes.tsv             one row per gene shown in the report
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/prep_regionplots.py   (from the repo root)
"""
import csv
import os
import re
import sys

ROOT = os.getcwd()
B = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(B, "tables")
FLANK = 100_000
M = os.path.join(ROOT, "Trifolium86MainK5ERpc5_results/GEA/tables/methods")
SIG = {"EMMAX": "EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv", "RDA": "RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv"}
TOPR_SRC = os.path.join(ROOT, "ANALYSIS/Trifolium86Main/deepdive6/figures/input/topr_ids.tsv")

# plot labels: key-table genes (every tandem copy), sig-SNP host genes, priority-A guesses in the window
LABELS = {
    # 1
    "g760": "BURP", "g761": "BURP", "g762": "BURP", "g763": "BURP (RD22 cl.)", "g764": "BURP (RD22 cl.)",
    "g765": "BURP (RD22 cl.)", "g766": "YUC6-like", "g779": "BOBBER1", "g782": "WSD", "g767": "DPH2-like",
    "g755": "RD22-like*", "g757": "RD22-like*",
    # 2
    "g25525": "PP2C clade B", "g25535": "phytol kinase", "g25536": "phytol kinase", "g25537": "phytol kinase",
    "g25527": "GDSL lipase",
    # 3
    "g26553": "caleosin PXG2-like", "g26535": "ETR", "g26545": "CCCH ZnF",
    # 4
    "g26949": "USP", "g26938": "LEA_3", "g26939": "LEA_3", "g26961": "TSPO", "g26944": "MAG lipase",
    # 5
    "g27359": "MYB37-clade", "g27357": "TRL33", "g27362": "GT21", "g27338": "TIP*", "g27342": "beta-CA*",
    # 6
    "g30974": "GSTL1",
    # 7
    "g31323": "B3 REM/VRN1", "g31328": "B3 REM/VRN1", "g31329": "B3 REM/VRN1", "g31337": "RRP12-like",
    "g31338": "MORN",
    # 8
    "g32205": "ClpB3", "g32221": "SnRK1 KIN10", "g32199": "BAG", "g32211": "OMT", "g32225": "OSCA-like*",
    # 9
    "g32874": "CCD8", "g32870": "RCC1", "g32871": "RCC1", "g32865": "GATA", "g32889": "GLO1*",
    # user 2026-10-17 (Fig 3 v3 Phase 4 evidence pass, regions/evidence_rgt12.md): label the resolved orthologues
    "g32863": "SKM1", "g32856": "MLK3", "g32864": "ATHB8",
    # 10
    "g38787": "ACA11-like", "g38779": "TIP", "g38778": "NOA1", "g38772": "mt sHSP",
    # 12
    "g45753": "NHL LEA_2", "g45752": "Yip1", "g45772": "UGT78-like", "g45773": "UGT78-like", "g45777": "CIPK21",
    "g45770": "PB1 kinase", "g45771": "MYB TF", "g45735": "VIN3-like*", "g45784": "DGK*",
    # new at 5 PCs
    "g12454": "GAAP4-like", "g12464": "MSRB5-like", "g12461": "SAMDC4", "g12462": "PP2C cl. L",
    "g12459": "MO25", "g12460": "F-box LRR",
    "g15018": "CNL", "g15021": "CNL", "g15022": "CNL", "g15023": "CNL frag.", "g15024": "CNL", "g15036": "CNL",
    "g15041": "CNL", "g15043": "CNL", "g15020": "DEAD-box",
    "g27356": "FADS",
    "g31011": "HYH-like", "g31005": "unannotated",
    "g33151": "ClpB4", "g33152": "ClpB4",
}


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gid(x):
    return x.split(".")[0]


ev = rows(os.path.join(T, "region_gene_evidence.tsv"))
regions = []
for r in ev:
    if r["region_id"] not in regions:
        regions.append(r["region_id"])
assert len(regions) == 15
ev_by = {(r["region_id"], gid(r["gene_id"])): r for r in ev}
coord = {gid(r["gene_symbol"]): (r["chrom"], int(r["gene_start"]), int(r["gene_end"])) for r in rows(TOPR_SRC)}

# key-table genes per region, parsed from REGIONS_ER.md (the report the user already has)
key = {}
cur = None
for line in open(os.path.join(B, "REGIONS_v2.md")):
    m = re.match(r"## \d+\. (\S+)", line)
    if line.startswith("## "):
        cur = m.group(1) if m and m.group(1) in regions else None
    if cur and line.startswith("| g"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        key.setdefault(cur, []).append((re.findall(r"g\d+", cells[0]), cells))
host = {}
for r in rows(os.path.join(T, "region_sig_snps.tsv")):
    for h in r["gene_hit"].split(","):
        mm = re.match(r"(g\d+)\.t1\((body|promoter)\)", h)
        if mm:
            host.setdefault(r["region_id"], {}).setdefault(mm.group(1), set()).add(r["snp_id"])
guess = [r for r in rows(os.path.join(T, "neighbourhood_guess.tsv")) if r["priority"] == "A"]
EXTRA = {"5_5557888-5883272": ["g32863", "g32856", "g32864"]}
# genes kept in LABELS / the report but drawn without a plot label (user 2026-10-18: no ATHB8 label in Fig 4b; the
# gene box stays)
PLOT_UNLABELLED = {"g32864"}
extra_ev = {r["gene_id"]: r for r in rows(os.path.join(B, "evidence_rgt12.tsv"))}

sig = {m: rows(os.path.join(M, f)) for m, f in SIG.items()}
thr = {m: float(tab[0]["pval_threshold"]) for m, tab in sig.items()}

calls, genes, used = [], [], set()
for i, rid in enumerate(regions, 1):
    c, s, e = re.match(r"(\d+)_(\d+)-(\d+)", rid).groups()
    s, e = int(s), int(e)
    ws, we = max(0, s - FLANK), e + FLANK
    traits = {}
    for m, tab in sig.items():
        for x in tab:
            if x["chr"] == c and ws <= int(x["pos"]) <= we:
                traits.setdefault(x["trait"], set()).add(m)
    arg = ",".join(f"{t}:{'|'.join(sorted(ms))}" for t, ms in sorted(traits.items()))
    calls.append({"order": i, "region_id": rid, "custom_region": f"{c}:{ws}-{we}", "traits_arg": arg,
                  "vlines": f"{s},{e}", "n_series": sum(len(v) for v in traits.values())})

    def add(g, where, label, ident, dist, r=None, pm=""):
        cc = coord[g]
        inwin = cc[0] == c and cc[2] >= ws and cc[1] <= we
        if label and not inwin:
            sys.exit(f"ERROR: labelled gene {g} outside the window of {rid}")
        genes.append({
            "order": i, "region_id": rid, "gene_id": g, "where": where, "label": label, "identity": ident,
            "distance": dist, "start": cc[1], "end": cc[2],
            "snps_body": r["n_snps_body"] if r else "", "snps_promoter": r["n_snps_promoter"] if r else "",
            "best_p_EMMAX_body": r["best_p_EMMAX_body"] if r else "",
            "best_p_EMMAX_promoter": r["best_p_EMMAX_promoter"] if r else "",
            "best_p_RDA_body": r["best_p_RDA_body"] if r else "",
            "best_p_RDA_promoter": r["best_p_RDA_promoter"] if r else "", "pmids": pm})
        if label:
            used.add(g)

    seen = set()
    for gs, cells in key.get(rid, []):
        for g in gs:
            r = ev_by[(rid, g)]
            add(g, "key", LABELS.get(g, ""), cells[1], f"{int(r['dist_to_sig_snp']):,} bp to sig. SNP", r,
                cells[-1])
            seen.add(g)
    for g, snps in sorted(host.get(rid, {}).items()):
        if g in seen:
            continue
        r = ev_by[(rid, g)]
        add(g, "sig-SNP host", LABELS.get(g, ""), r["identity"] or r["desc"][:70],
            "holds " + ", ".join(sorted(snps)), r)
        seen.add(g)
    # genes added by the rgt_12 evidence pass (regions/evidence_rgt12.tsv; user 2026-10-17: label them)
    for g in EXTRA.get(rid, []):
        r, x = ev_by[(rid, g)], extra_ev[g]
        add(g, "evidence pass", LABELS[g], x["identity"], f"{int(r['dist_to_sig_snp']):,} bp to sig. SNP", r,
            ",".join(p.split("(")[0] for p in (x["new_pmids_key"] + "," + x["new_pmids_support"]).split(",") if p))
        seen.add(g)
    for r in guess:
        if rid not in r["region_ids"].split(","):
            continue
        g = gid(r["gene_id"])
        d = dict(zip(r["region_ids"].split(","), r["dist_to_edge"].split(",")))[rid]
        side, kb = ("left" if d[0] == "L" else "right"), int(d[1:]) / 1000
        cc = coord[g]
        inwin = cc[0] == c and cc[2] >= ws and cc[1] <= we
        add(g, "outside (guess)", LABELS.get(g, "") if inwin else "", r["guess"], f"{kb:,.1f} kb {side} of edge")

# every label used exactly where it belongs
missing = set(LABELS) - used
if missing:
    sys.exit("ERROR: labels never placed: " + ", ".join(sorted(missing)))
for rid in regions:
    for gs, _ in key.get(rid, []):
        for g in gs:
            assert g in LABELS, g
    for g in host.get(rid, {}):
        assert g in LABELS, g

# ---- topr annotation ----
src = rows(TOPR_SRC)
n_lab = 0
for r in src:
    g = gid(r["gene_symbol"])
    r["biotype"] = "protein_coding"
    # plot label without brackets or "*" (user review v2, 2026-10-17: "remove brackets from the gene names",
    # "* no needed"); report/region_genes.tsv keeps the full LABELS text
    r["gene_symbol"] = "" if g in PLOT_UNLABELLED else re.sub(r"\s*\([^)]*\)|\*", "", LABELS.get(g, "")).strip()
    n_lab += g in LABELS
assert n_lab == len(LABELS)
os.makedirs(os.path.join(B, "figures", "input"), exist_ok=True)
os.makedirs(os.path.join(B, "report"), exist_ok=True)
for path, data in [(os.path.join(B, "figures", "input", "topr_labelled.tsv"), src),
                   (os.path.join(B, "figures", "input", "regionplot_calls.tsv"), calls),
                   (os.path.join(B, "report", "region_genes.tsv"), genes)]:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(data[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(data)
print(f"INFO: thresholds {thr}; {len(LABELS)} labelled genes; {len(genes)} report gene rows")
for c in calls:
    print(f"  {c['order']:>2} {c['region_id']:<24} {c['custom_region']:<26} n={c['n_series']} {c['traits_arg']}")
