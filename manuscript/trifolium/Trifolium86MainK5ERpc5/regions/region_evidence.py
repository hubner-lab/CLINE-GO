#!/usr/bin/env python3
"""Per-gene evidence table of the 15 EMMAX+RDA regions of Trifolium86MainK5ERpc5 (5 PCs)
(plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7). COPY of ANALYSIS/Trifolium86MainK5/regionsER/
region_evidence.py (PC3 run). IDENTITY / EVIDENCE / FG_GENE keep the PC3 entries (filtered to genes still in a region:
the entries of the two dropped chr-6 regions are left out and counted) and add NEW_* entries for the Phase 7
literature genes.

One row per region gene (all 353). Columns:
  annotation      desc, preferred_name, pfam, ko_names, array_id / array_size
  curation        call / category / reason / source / lit (tables/curation.tsv)
  identity        literature-pass identity of the gene (IDENTITY below, from the literature agents' NCBI protein / UniProt
                  look-ups of the eggNOG best hit; "" = not looked up). It can contradict the curation call: that is
                  reported, the call is not rewritten.
  function_group  declared from category + identity before any text was written (FG below)
  evidence        the agents' verdict for the gene class (strong / moderate / weak / none; "" = not sent)
  SNPs            n SNPs and best EMMAX / RDA p in body and promoter, sig SNPs there, nearest significant SNP, distance,
                  its methods, carrier signature, pattern_n, locus_unique (tables/gene_snps.tsv + sig_snp_patterns.tsv)
  literature      PMIDs by role (key / support / conflict) from tables/literature.tsv (verified); rejected rows excluded;
                  n_refs, n_legume_refs, key finding
Output: tables/region_gene_evidence.tsv
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/region_evidence.py
"""
import collections
import csv
import os

B = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(B, "tables")

IDENTITY = {
    "g760.t1": "BURP; best hit MtBURP9 (RD22 class not stated)", "g761.t1": "BURP; best hit MtBURP9",
    "g762.t1": "BURP; best hit MtBURP6, not in the RD22-like class of Li 2016 (contests the RD22-like call)",
    "g763.t1": "BURP; best hit MtBURP5, RD22-like class", "g764.t1": "BURP; best hit MtBURP5, RD22-like class",
    "g765.t1": "BURP; best hit MtBURP5, RD22-like class",
    "g766.t1": "YUCCA; best hit MtYUC1, clade I-2, close to AtYUC6",
    "g778.t1": "CCT-only protein (no B-box): CO-like membership uncertain",
    "g779.t1": "BOBBER1 (NudC-domain non-canonical sHSP)",
    "g782.t1": "WSD-family wax ester synthase (member unresolved)",
    "g25525.t1": "clade B PP2C, AP2C1-like MAPK phosphatase (not an ABA clade-A PP2C)",
    "g25535.t1": "phytol / farnesol kinase (FOLK vs VTE5 unresolved)", "g25536.t1": "phytol / farnesol kinase",
    "g25537.t1": "phytol / farnesol kinase",
    "g26535.t1": "ethylene receptor (subfamily unresolved)",
    "g26553.t1": "caleosin / peroxygenase; best hit PXG2-like (CLO2), not RD20 / CLO3",
    "g26938.t1": "ARG2-like atypical LEA (Pfam LEA_3)", "g26939.t1": "ARG2-like atypical LEA (Pfam LEA_3)",
    "g26948.t1": "EID1-like F-box (EID1 / EDL clade)",
    "g26949.t1": "universal stress protein, At3g01520-like clade",
    "g26961.t1": "TSPO (TspO/MBR) membrane protein",
    "g27357.t1": "thioredoxin-like 3-3 (TRL33 clade), not an h-type thioredoxin",
    "g27359.t1": "R2R3 MYB subgroup 14 (RAX / MYB37 clade)",
    "g30974.t1": "lambda-class GST (GSTL1-like)",
    "g31323.t1": "B3, REM / VRN1 clade (two B3 domains)", "g31328.t1": "B3, REM / VRN1 clade",
    "g31329.t1": "B3, REM / VRN1 clade",
    "g32199.t1": "BAG co-chaperone (member unresolved)",
    "g32205.t1": "ClpB3 / APG6, chloroplast HSP100",
    "g32206.t1": "EID1-like F-box (EID1 / EDL unresolved)",
    "g32221.t1": "SnRK1 catalytic alpha, KIN10-like",
    "g32870.t1": "RCC1-family; UVR8 identity unresolved (best hit could not be checked)",
    "g32871.t1": "RCC1-family, NOT UVR8 (best hit lacks the Trp-triad GWRHT motifs)",
    "g32874.t1": "CCD8 / MAX4 (strigolactone biosynthesis)",
    "g38772.t1": "mitochondrial small HSP (HSP22 / MT-sHSP type)",
    "g38778.t1": "NOA1 / RIF1 plastid GTPase",
    "g38779.t1": "TIP aquaporin (TIP subfamily unresolved)",
    "g38787.t1": "P2B autoinhibited Ca2+-ATPase (ACA11-like), not an ATP2C pump",
    "g44931.t1": "CAS chloroplast calcium-sensing receptor",
    "g45753.t1": "NDR1/HIN1-like (NHL) subgroup of LEA_2",
    "g45772.t1": "flavonoid glucosyltransferase (UGT78-like; 3-O vs 7-O unresolved)",
    "g45773.t1": "flavonoid glucosyltransferase (UGT78-like)",
    "g45777.t1": "CIPK (CIPK21 by eggNOG name)",
    "g45780.t1": "CONSTANS-like, two B-boxes + CCT (member unresolved)",
    "g45794.t1": "SMAX1 / SMXL2 clade (not ClpC / ClpD)",
    "g45795.t1": "trehalose-6-phosphate synthase with TPP domain (class unresolved)",
    "g45808.t1": "CER3 / WAX2 (alkane-forming wax biosynthesis)",
    "g45817.t1": "FLX-like 4 (FLL4)",
}
EVIDENCE = {  # agents' verdict for the gene class (decision tables, 2026-10-15)
    "g760.t1": "weak", "g761.t1": "weak", "g762.t1": "weak", "g763.t1": "moderate", "g764.t1": "moderate",
    "g765.t1": "moderate", "g766.t1": "moderate", "g778.t1": "weak", "g779.t1": "moderate", "g782.t1": "moderate",
    "g25525.t1": "weak", "g25535.t1": "weak", "g25536.t1": "weak", "g25537.t1": "weak", "g26535.t1": "weak",
    "g26553.t1": "weak", "g26938.t1": "moderate", "g26939.t1": "moderate", "g26948.t1": "weak",
    "g26949.t1": "weak", "g26961.t1": "weak", "g27357.t1": "none", "g27359.t1": "moderate", "g30974.t1": "moderate",
    "g31323.t1": "weak", "g31328.t1": "weak", "g31329.t1": "weak", "g32199.t1": "weak", "g32205.t1": "moderate",
    "g32206.t1": "weak", "g32221.t1": "moderate", "g32870.t1": "none", "g32871.t1": "none", "g32874.t1": "strong",
    "g38772.t1": "moderate", "g38778.t1": "moderate", "g38779.t1": "moderate", "g38787.t1": "moderate",
    "g44931.t1": "moderate", "g45753.t1": "moderate", "g45772.t1": "weak", "g45773.t1": "weak",
    "g45777.t1": "moderate", "g45780.t1": "weak", "g45794.t1": "strong", "g45795.t1": "moderate",
    "g45808.t1": "moderate", "g45817.t1": "weak",
}
FG = {  # function group, declared per gene before writing
    "drought_osmotic": "dehydration / osmotic protection", "heat": "heat protection", "cold": "cold",
    "salt_ion": "salt / ion and Ca2+ signalling", "oxidative": "ROS homeostasis", "aba": "ABA / stomatal water use",
    "water_transport": "water transport", "uv_light": "UV / light",
}
FG_GENE = {
    "g766.t1": "drought (YUC6 thiol reductase)", "g782.t1": "cuticle / wax", "g45808.t1": "cuticle / wax",
    "g26553.t1": "ABA / stomatal water use", "g44931.t1": "ABA / stomatal water use",
    "g38787.t1": "salt / ion and Ca2+ signalling", "g45777.t1": "salt / ion and Ca2+ signalling",
    "g32874.t1": "strigolactone (drought / temperature)", "g45794.t1": "strigolactone / karrikin (drought)",
    "g32221.t1": "energy-stress signalling", "g38778.t1": "chilling / NO", "g45753.t1": "salt / ABA (NHL)",
    "g27359.t1": "ABA / drought (MYB37 clade)", "g25525.t1": "MAPK phosphatase (ABA, weak)",
    "g32870.t1": "unresolved (RCC1)", "g32871.t1": "unresolved (RCC1)", "g27357.t1": "unresolved (TRL33)",
}
NEW_FG_GENE = {  # hand-set after reading the Phase 7 agents' decision tables
    "g12454.t1": "drought / salt (GAAP4)", "g12461.t1": "polyamine synthesis (development only)",
    "g12462.t1": "light signalling (clade L PP2C, not ABA)", "g12464.t1": "ROS homeostasis",
    "g33151.t1": "heat protection (expression only)", "g33152.t1": "heat protection (expression only)",
}


def _ident():
    """NEW_IDENTITY / NEW_EVIDENCE from the Phase 7 agents' identity tables (new lit genes + defence-class genes)"""
    idn, evd = {}, {}
    for f in ("agent_batch_v2_identity.tsv", "agent_batch_biotic_A_identity.tsv", "agent_batch_biotic_B_identity.tsv"):
        with open(os.path.join(B, "logs", f)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = r["gene_id"] if r["gene_id"].endswith(".t1") else r["gene_id"] + ".t1"
                assert g not in idn, g
                assert r["evidence_verdict"] in ("strong", "moderate", "weak", "none"), r
                idn[g], evd[g] = r["identity"], r["evidence_verdict"]
    return idn, evd


NEW_IDENTITY, NEW_EVIDENCE = _ident()


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


sh = rows(os.path.join(T, "gene_sheet.tsv"))
cu = {(r["region_id"], r["gene_id"]): r for r in rows(os.path.join(T, "curation.tsv"))}
gs = {(r["region_id"], r["gene_id"]): r for r in rows(os.path.join(T, "gene_snps.tsv"))}
pt = {(r["region_id"], r["snp_id"]): r for r in rows(os.path.join(T, "sig_snp_patterns.tsv"))}
lit = rows(os.path.join(T, "literature.tsv"))
here = {r["gene_id"] for r in sh}
n_out = len(set(IDENTITY) - here)
IDENTITY = {k: v for k, v in IDENTITY.items() if k in here}
EVIDENCE = {k: v for k, v in EVIDENCE.items() if k in here}
FG_GENE = {k: v for k, v in FG_GENE.items() if k in here}
assert not (set(NEW_IDENTITY) & set(IDENTITY))
IDENTITY.update(NEW_IDENTITY)
EVIDENCE.update(NEW_EVIDENCE)
FG_GENE.update(NEW_FG_GENE)
print(f"INFO: {n_out} PC3 identity entries left out with the dropped regions; {len(NEW_IDENTITY)} new")
assert set(IDENTITY) <= here and set(EVIDENCE) <= set(IDENTITY)
L = collections.defaultdict(list)
for r in lit:
    if r["role"] != "rejected":
        L[r["gene_id"] + ".t1"].append(r)
assert set(L) <= set(IDENTITY), set(L) - set(IDENTITY)

out = []
for s in sh:
    k = (s["region_id"], s["gene_id"])
    c, g = cu[k], gs[k]
    p = pt[(s["region_id"], g["nearest_sig_snp"])]
    if c["call"] == "phenology":
        fg = "phenology"
    elif c["call"] == "reject" or not c["call"]:
        fg = ""
    else:
        fg = FG_GENE.get(s["gene_id"], FG.get(c["category"], "candidate (class not resolved)"))
    ls = L.get(s["gene_id"], [])
    byrole = lambda role: ",".join(f"{x['pmid']}({x['first_author']} {x['year']})" for x in ls if x["role"] == role)
    key = next((x for x in ls if x["role"] == "key"), None)
    out.append({
        "region_id": s["region_id"], "gene_id": s["gene_id"], "start": s["start"], "end": s["end"],
        "strand": s["strand"], "desc": s["desc"], "preferred_name": s["preferred_name"], "pfam": s["pfam"],
        "ko_names": s["ko_names"], "array_id": s["array_id"], "array_size": s["array_size"],
        "call": c["call"], "category": c["category"], "source": c["source"], "lit_sent": c["lit"],
        "curation_reason": c["reason"], "identity": IDENTITY.get(s["gene_id"], ""), "function_group": fg,
        "evidence": EVIDENCE.get(s["gene_id"], ""),
        "n_snps_body": g["n_snps_body"], "best_p_EMMAX_body": g["best_p_EMMAX_body"],
        "best_p_RDA_body": g["best_p_RDA_body"], "sig_snps_body": g["sig_snps_body"],
        "n_snps_promoter": g["n_snps_promoter"], "best_p_EMMAX_promoter": g["best_p_EMMAX_promoter"],
        "best_p_RDA_promoter": g["best_p_RDA_promoter"], "sig_snps_promoter": g["sig_snps_promoter"],
        "nearest_sig_snp": g["nearest_sig_snp"], "dist_to_sig_snp": g["dist_to_sig_snp"],
        "nearest_sig_methods": g["nearest_sig_methods"], "nearest_sig_signature": p["carrier_signature"],
        "nearest_sig_pattern_n": p["pattern_n"], "nearest_sig_locus_unique": p["locus_unique"],
        "pmids_key": byrole("key"), "pmids_support": byrole("support"), "pmids_conflict": byrole("conflict"),
        "n_refs": len(ls), "n_legume_refs": sum(1 for x in ls if x["legume"] == "yes"),
        "key_finding": key["finding"] if key else "",
    })
with open(os.path.join(T, "region_gene_evidence.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, list(out[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
print(f"INFO: {len(out)} rows; {sum(1 for r in out if r['n_refs'])} genes with verified literature; "
      f"{sum(1 for r in out if r['function_group'])} with a function group")
