#!/usr/bin/env python3
"""Per-gene evidence inputs for the 15 EMMAX+RDA (2/2) regions of Trifolium86MainK5ERpc5 (5 PCs)
(plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7). COPY of ANALYSIS/Trifolium86MainK5/regionsER/
evidence_inputs.py (the PC3 run): same logic; inputs switched to the PC5 project, and the "old call" columns now carry
the PC3 regionsER curation (tables/curation.tsv there, keyed by gene_id: two PC3 regions changed id at 5 PCs), not the
genes3of3 one. The ``old_map`` against the old 36 is kept as in v1.

Regions / genes: the pipeline's own regions_combined.tsv / genes_combined.tsv of Trifolium86MainK5ERpc5.
Annotation: read from the GFF BY KEY (genes_combined.tsv columns after em_max_annot_lvl are misaligned, defect D7).
The GFF is byte-identical to Trifolium86Main's, so gene ids carry over from the old curation (asserted by the caller
with cmp; asserted here by the old gene_sheet ids all resolving).

Flags: the strict / verbatim keyword lines and the tandem-array ids are taken from the screen this set already went
through (Trifolium86MainK5ERpc5/regions/stress/genes_stress_annotated.tsv, = stress_screen.R, Phase 2). For the +-500 kb
neighbourhood genes the same STRICT keyword regexes are re-applied here (copied from stress_screen.R FAM_STRICT);
the GO line is omitted there (a guess list, not a screen).

Significant SNPs of a region: every SNP of the EMMAX bonf 0.05 or RDA qval 0.05 sig table inside the region span;
`agreed` marks the 30 regions_combined.tsv snp_ids (the 2/2 Overlap members).

Promoter = 5,000 bp, strand-aware as the pipeline does it: [start - L, start] on +, [end, end + L] on -.

Outputs (tables/):
  gene_sheet.tsv        one row per region gene: annotation, KEGG names, array, flags, old curation call + old region
  old_map.tsv           per A region: overlapping old-36 regions, genes kept / new / lost against them
  gene_snps.tsv         per region gene: tested SNPs in body / promoter, best p per method (+ predictor), nearest
                        significant SNP of the region and its distance
  region_sig_snps.tsv   per significant SNP: p per method, agreed flag, gene it sits in (body/promoter) or nearest gene
  neighbourhood.tsv     every transcript within region midpoint +- 500,000 bp that is not a gene of that region
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/evidence_inputs.py   (from the repo root)
"""
import collections
import csv
import os
import re
import sys

ROOT = os.getcwd()
B = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(ROOT, "Trifolium86MainK5ERpc5_results")
T = os.path.join(P, "GEA/tables")
M = os.path.join(T, "methods")
GFF = os.path.join(P, "_intermediate/annotation/normalized.gff3")
STRESS = os.path.join(ROOT, "ANALYSIS/Trifolium86MainK5ERpc5/regions/stress/genes_stress_annotated.tsv")
V1CUR = os.path.join(ROOT, "ANALYSIS/Trifolium86MainK5/regionsER/tables/curation.tsv")
G3 = os.path.join(ROOT, "ANALYSIS/Trifolium86Main/genes3of3")
KO = os.path.join(G3, "input/kegg_ko_list.tsv")
OUT = os.path.join(B, "tables")
PROM = 5000
NBH = 500000
os.makedirs(OUT, exist_ok=True)

FAM_STRICT = {  # copied from ANALYSIS/Trifolium86MainK5/stress_screen.R FAM_STRICT
    "drought_osmotic": r"dehydration|dehydrin|drought|osmot|desicc|LEA_|late embryogenesis|aquaporin|(?<!hydroxy)proline(?!-rich)|betaine|raffinose|trehalose",
    "heat_chaperone": r"heat shock|\bHSP\d|HSP20|HSP70|HSP90|chaperonin|Cpn60|DnaJ|ClpB|thermotoleran",
    "cold": r"cold|freez|chill|\bCBF\d|DREB1|ICE1|antifreeze",
    "salt_ion": r"salt|salin|Na\+/H\+|\bNHX|SOS1|\bHKT|cation/H\+ exchanger|potassium transport",
    "oxidative": r"peroxidase|catalase|superoxide dismutase|glutathione (S-transferase|peroxidase|reductase)|thioredoxin|monodehydroascorbate|dehydroascorbate reductase|ascorbate peroxidase|oxidative stress",
    "aba_signalling": r"abscisic|ABA |PYR/PYL|PP2C|SnRK2|ABF|AREB",
    "stress_tf": r"DREB|WRKY|\bNAC\b|NAC domain|\bMYB\b|bZIP|leucin zipper|AP2/ERF|ethylene-responsive transcription|heat stress transcription|ABSCISIC ACID-INSENSITIVE|G-box-binding|\bWRKY|\bAP2\b|\bNAM\b|HSF_DNA-bind|bZIP_",
    "water_transport": r"aquaporin|plasma membrane intrinsic|tonoplast intrinsic|\bMIP\b|\b(PIP|TIP)[12]\b",
    "flowering_photo": r"flowering|photoperiod|circadian|CONSTANS|FT |vernaliz|phytochrome|cryptochrome",
}
FAM_RE = {k: re.compile(v, re.I) for k, v in FAM_STRICT.items()}


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, header, data):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(data)


def attr(a, k):
    m = re.search(r"(?:^|;)" + re.escape(k) + r"=([^;]*)", a)
    return m.group(1) if m else ""


# ---- GFF transcripts, by key -----------------------------------------------------------------------------
gff = {}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        a = f[8]
        gid = attr(a, "ID")
        assert gid not in gff
        gff[gid] = dict(chr=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], desc=attr(a, "em_desc"),
                        preferred_name=attr(a, "em_Preferred_name"), pfam=attr(a, "em_PFAMs"),
                        kegg_ko=attr(a, "em_KEGG_ko").replace("ko:", ""), ogs=attr(a, "em_OGs"))
print(f"INFO: {len(gff)} GFF transcripts")
ko_name = {r[0].replace("ko:", ""): r[1] for r in csv.reader(open(KO), delimiter="\t")}


def ko_names(kos):
    return " | ".join(f"{k} {ko_name.get(k, '?')}" for k in kos.split(",") if k)


def og_narrow(ogs):
    return ogs.split(",")[-1].split("|")[0] if ogs else ""


# ---- regions and genes of the PC5 set --------------------------------------------------------------------------
reg = rows(os.path.join(T, "regions_combined.tsv"))
assert len(reg) == 15
R = {r["region_id"]: dict(chr=r["chr"], start=int(r["start"]), end=int(r["end"]),
                          agreed=r["snp_ids"].split(","), traits=r["traits"], min_p=r["min_pvalue"]) for r in reg}
assert sum(len(r["agreed"]) for r in R.values()) == 30
gc = rows(os.path.join(T, "genes_combined.tsv"))
assert len(gc) == 353
rg = collections.defaultdict(list)
for r in gc:
    assert r["gene_id"] in gff, r["gene_id"]
    rg[r["region_id"]].append(r["gene_id"])
assert set(rg) == set(R)

st = {(r["region_id"], r["gene_id"]): r for r in rows(STRESS)}
assert len(st) == 353 and set(st) == {(r["region_id"], r["gene_id"]) for r in gc}
KS = [k for k in rows(STRESS)[0] if k.startswith("ks_") and k != "ks_any"]
KW = [k for k in rows(STRESS)[0] if k.startswith("kw_") and k != "kw_any"]
GOC = [k for k in rows(STRESS)[0] if k.startswith("go_GO") or (k.startswith("go_") and k != "go_any")]

# ---- old-36 curation and region membership ---------------------------------------------------------------
old_sheet = rows(os.path.join(G3, "tables/gene_sheet.tsv"))
assert all(r["gene_id"] in gff for r in old_sheet)
old_reg = {r["region_id"]: (r["chr"], int(r["start"]), int(r["end"]), r["rank"], r["tier"])
           for r in rows(os.path.join(G3, "tables/regions_ranked.tsv"))}
assert len(old_reg) == 36
old_genes = collections.defaultdict(set)
for r in old_sheet:
    old_genes[r["region_id"]].add(r["gene_id"])
cur = {}
for r in rows(V1CUR):
    assert r["gene_id"] not in cur
    cur[r["gene_id"]] = r  # PC3 regionsER curation (calls, rule source, lit flag); a gene sits in one PC3 region
assert len(cur) == 315
cur3 = {r["gene_id"]: r for r in rows(os.path.join(G3, "tables/curation.tsv"))}  # old-36 calls, neighbourhood column only
assert len(cur3) == 307

# ---- gene sheet + old map --------------------------------------------------------------------------------
sheet, omap = [], []
for rid in sorted(R, key=lambda x: (R[x]["chr"], R[x]["start"])):
    c, s, e = R[rid]["chr"], R[rid]["start"], R[rid]["end"]
    ov = [o for o, (oc, os_, oe, _, _) in old_reg.items() if oc == c and os_ <= e and oe >= s]
    og = set().union(*[old_genes[o] for o in ov]) if ov else set()
    ag = set(rg[rid])
    omap.append([rid, ",".join(ov), ",".join(f"{old_reg[o][3]}:{old_reg[o][4]}" for o in ov), len(ag),
                 len(ag & og), len(ag - og), len(og - ag), ",".join(sorted(og - ag))])
    for gid in sorted(rg[rid], key=lambda g: gff[g]["start"]):
        g, f = gff[gid], st[(rid, gid)]
        cu = cur.get(gid)
        sheet.append([rid, gid, g["chr"], g["start"], g["end"], g["strand"], g["desc"], g["preferred_name"], g["pfam"],
                      g["kegg_ko"], ko_names(g["kegg_ko"]), og_narrow(g["ogs"]), f["array_id"], f["array_size"],
                      ",".join(k[3:] for k in KS if f[k] == "TRUE"), ",".join(k[3:] for k in KW if f[k] == "TRUE"),
                      ",".join(k[3:] for k in GOC if f[k] == "TRUE"),
                      "TRUE" if "TRUE" in (f["kw_any"], f["ks_any"], f["go_any"]) else "FALSE",
                      cu["call"] if cu else "", cu["category"] if cu else "", cu["reason"] if cu else "",
                      cu["region_id"] if cu else ""])
write("gene_sheet.tsv", ["region_id", "gene_id", "chr", "start", "end", "strand", "desc", "preferred_name", "pfam",
                         "kegg_ko", "ko_names", "og_narrow", "array_id", "array_size", "ks_lines", "kw_lines",
                         "go_terms", "flagged", "old_call", "old_category", "old_reason", "old_region"], sheet)
write("old_map.tsv", ["region_id", "old_regions", "old_rank_tier", "n_genes", "n_kept", "n_new", "n_lost",
                      "lost_genes"], omap)
print(f"INFO: gene_sheet {len(sheet)} rows; {sum(1 for r in sheet if r[17] == 'TRUE')} flagged; "
      f"{sum(1 for r in sheet if r[18])} carry an old call")

# ---- tested SNPs inside the regions and p values ---------------------------------------------------------
def in_region(c, p):
    return [rid for rid, r in R.items() if r["chr"] == c and r["start"] <= p <= r["end"]]


snp = {}  # (rid, sid) -> dict
with open(os.path.join(M, "EMMAX/EMMAX_pvalues_K5.tsv")) as fh:
    rd = csv.reader(fh, delimiter="\t")
    hdr = next(rd)
    pc = hdr[3:]
    n_tested = 0
    for f in rd:
        n_tested += 1
        c, p = f[1], int(f[2])
        for rid in in_region(c, p):
            vals = [float(x) if x not in ("", "NA") else 9.0 for x in f[3:]]
            mp = min(vals)
            snp[(rid, f[0])] = dict(chr=c, pos=p, EMMAX=mp, EMMAX_pred=",".join(pc[i] for i, v in enumerate(vals) if v == mp))
assert n_tested == 177259, n_tested
with open(os.path.join(M, "RDA/RDA_pvalues_K5.tsv")) as fh:
    rd = csv.reader(fh, delimiter="\t")
    next(rd)
    n_rda = 0
    for f in rd:
        n_rda += 1
        for rid in in_region(f[1], int(f[2])):
            snp[(rid, f[0])]["RDA"] = float(f[3]) if f[3] not in ("", "NA") else 9.0
assert n_rda == 177259 and all("RDA" in v for v in snp.values())
print(f"INFO: {len(snp)} tested SNPs inside the 15 regions")

thr, sig = {}, collections.defaultdict(set)
for m, fn in (("EMMAX", "EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv"), ("RDA", "RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv")):
    s = rows(os.path.join(M, fn))
    t = {float(r["pval_threshold"]) for r in s}
    assert len(t) == 1
    thr[m] = t.pop()
    for r in s:
        sig[r["SNPID"]].add(m)
# sig-table membership equals min p <= threshold inside the regions
for (rid, sid), v in snp.items():
    for m in ("EMMAX", "RDA"):
        assert (v[m] <= thr[m]) == (m in sig.get(sid, set())), (sid, m)
print(f"INFO: thresholds EMMAX {thr['EMMAX']:.3g}, RDA {thr['RDA']:.3g}")
rsig = {(rid, sid): v for (rid, sid), v in snp.items() if sid in sig}
agreed = {(rid, s) for rid, r in R.items() for s in r["agreed"]}
assert agreed <= set(rsig), agreed - set(rsig)


def iv(g):
    return (g["start"], g["end"]), ((g["end"], g["end"] + PROM) if g["strand"] == "-" else (max(1, g["start"] - PROM), g["start"]))


gs = []
for rid in sorted(R, key=lambda x: (R[x]["chr"], R[x]["start"])):
    sigpos = [(v["pos"], sid) for (r2, sid), v in rsig.items() if r2 == rid]
    for gid in sorted(rg[rid], key=lambda g: gff[g]["start"]):
        g = gff[gid]
        body, prom = iv(g)
        out = [rid, gid]
        for lo, hi in (body, prom):
            hs = [v for (r2, sid), v in snp.items() if r2 == rid and lo <= v["pos"] <= hi]
            be = min(hs, key=lambda v: v["EMMAX"]) if hs else None
            out += [len(hs), f"{be['EMMAX']:.3g}" if be else "", be["EMMAX_pred"] if be else "",
                    f"{min(v['RDA'] for v in hs):.3g}" if hs else "",
                    ",".join(sorted(sid for (r2, sid), v in rsig.items() if r2 == rid and lo <= v["pos"] <= hi))]
        d = [(0 if g["start"] <= p <= g["end"] else min(abs(g["start"] - p), abs(g["end"] - p)), sid) for p, sid in sigpos]
        dmin = min(d)
        out += [dmin[1], dmin[0], ",".join(sorted(sig[dmin[1]]))]
        gs.append(out)
write("gene_snps.tsv", ["region_id", "gene_id",
                        "n_snps_body", "best_p_EMMAX_body", "pred_EMMAX_body", "best_p_RDA_body", "sig_snps_body",
                        "n_snps_promoter", "best_p_EMMAX_promoter", "pred_EMMAX_promoter", "best_p_RDA_promoter",
                        "sig_snps_promoter", "nearest_sig_snp", "dist_to_sig_snp", "nearest_sig_methods"], gs)

# ---- significant SNP sheet -------------------------------------------------------------------------------
out = []
for (rid, sid), v in sorted(rsig.items(), key=lambda x: (x[1]["chr"], x[1]["pos"])):
    hits = []
    for gid in rg[rid]:
        body, prom = iv(gff[gid])
        if body[0] <= v["pos"] <= body[1]:
            hits.append(f"{gid}(body)")
        elif prom[0] <= v["pos"] <= prom[1]:
            hits.append(f"{gid}(promoter)")
    if not hits:
        dd = min((min(abs(gff[g]["start"] - v["pos"]), abs(gff[g]["end"] - v["pos"])), g) for g in rg[rid])
        hits = [f"{dd[1]}({dd[0]} bp)"]
    out.append([rid, sid, v["pos"], ",".join(sorted(sig[sid])), "TRUE" if (rid, sid) in agreed else "FALSE",
                f"{v['EMMAX']:.3g}", v["EMMAX_pred"], f"{v['RDA']:.3g}", ",".join(hits)])
write("region_sig_snps.tsv", ["region_id", "snp_id", "pos", "sig_methods", "agreed", "EMMAX_minp", "EMMAX_pred",
                              "RDA_p", "gene_hit"], out)
print(f"INFO: {len(out)} significant SNPs inside the 15 regions ({len(agreed)} agreed)")

# ---- +-500 kb neighbourhood ------------------------------------------------------------------------------
nb = collections.defaultdict(lambda: dict(regions=[], dists=[]))
for rid, r in R.items():
    mid = (r["start"] + r["end"]) // 2
    lo, hi = mid - NBH, mid + NBH
    mine = set(rg[rid])
    for gid, g in gff.items():
        if g["chr"] != r["chr"] or g["end"] < lo or g["start"] > hi or gid in mine:
            continue
        d = r["start"] - g["end"] if g["end"] < r["start"] else g["start"] - r["end"]
        nb[gid]["regions"].append(rid)
        nb[gid]["dists"].append(f"{'L' if g['end'] < r['start'] else 'R'}{d}")
in_any_region = {g for v in rg.values() for g in v}
nbo = []
for gid in sorted(nb, key=lambda g: (gff[g]["chr"], gff[g]["start"])):
    g = gff[gid]
    txt = f"{g['desc']} {g['pfam']}"
    ks = ",".join(k for k, rx in FAM_RE.items() if rx.search(txt))
    nbo.append([",".join(nb[gid]["regions"]), gid, g["chr"], g["start"], g["end"], ",".join(nb[gid]["dists"]),
                "TRUE" if gid in in_any_region else "FALSE", g["desc"], g["preferred_name"], g["pfam"],
                ko_names(g["kegg_ko"]), ks, cur3[gid]["call"] if gid in cur3 else ""])
write("neighbourhood.tsv", ["region_ids", "gene_id", "chr", "start", "end", "dist_to_edge", "in_other_region", "desc",
                            "preferred_name", "pfam", "ko_names", "ks_lines", "old_call"], nbo)
print(f"INFO: neighbourhood {len(nbo)} genes ({sum(1 for r in nbo if r[11])} strict-flagged; "
      f"{sum(1 for r in nbo if ',' in r[0])} shared by two regions)")
