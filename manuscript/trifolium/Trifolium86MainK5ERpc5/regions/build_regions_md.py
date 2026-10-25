#!/usr/bin/env python3
"""Builds REGIONS_v2.md from the hand-written REGIONS_v2.src.md (plan ~/.claude/plans/peaceful-kindling-micali.md,
Phase 7). The prose is written by hand; every number cell that the PC3 report typed by hand is filled here from the
v2 tables, so a carried PC3 sentence cannot bring a PC3 number with it:
  {{OVERVIEW}}           the region table (kb, genes, significant SNPs, best p per method, EMMAX predictors, support)
  {{KEY:<region_id>}}    the key-gene table of that region: one row per KEY[region] entry; the distance, SNP-count and
                         best-p cells come from tables/region_gene_evidence.tsv, the identity / function-group /
                         literature cells from KEY.
Then report_check.py checks the result.
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/build_regions_md.py
"""
import collections
import csv
import os
import re

B = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(B, "tables")
R = os.path.join(B, "../../../Trifolium86MainK5ERpc5_results/GEA/tables/regions_combined.tsv")

# region -> [(genes, identity, function group, literature cell)]
KEY = {
    "1_8649433-8911055": [
        ("g763,g764,g765", "BURP, best hit MtBURP5 (RD22-like class)", "dehydration / osmotic protection", "27148311, 22548236; conflict 25333723"),
        ("g760,g761,g762", "BURP, best hits MtBURP9 / MtBURP6 (not RD22-like)", "dehydration / osmotic protection", "family level only"),
        ("g766", "YUCCA (MtYUC1, close to AtYUC6)", "drought (YUC6 thiol reductase)", "26314500, 22986790"),
        ("g779", "BOBBER1 (NudC non-canonical sHSP)", "heat protection", "19571304"),
        ("g782", "WSD-family wax ester synthase", "cuticle / wax", "30729606, 34068347"),
    ],
    "3_147045758-147307380": [
        ("g25525", "clade B PP2C, AP2C1-like MAPK phosphatase", "MAPK phosphatase (ABA, weak)", "20488890"),
        ("g25535,g25536,g25537", "phytol / farnesol kinase (FOLK vs VTE5 unresolved)", "ABA / stomatal water use", "21395888"),
    ],
    "4_8145960-8407582": [
        ("g26553", "caleosin / peroxygenase, PXG2-like", "ABA / stomatal water use", "20952421; conflict 21471120"),
        ("g26535", "ethylene receptor (subfamily unresolved)", "candidate (class not resolved)", "24820022, 22706288"),
    ],
    "4_11490368-11751990": [
        ("g26949", "universal stress protein, At3g01520-like", "candidate", "22915741; conflict 35080785"),
        ("g26938,g26939", "ARG2-like atypical LEA (Pfam LEA_3)", "dehydration / osmotic protection", "41168775"),
        ("g26961", "TSPO", "dehydration / osmotic protection", "19548979; conflict 38623073"),
    ],
    "4_15212239-15556144": [
        ("g27359", "R2R3 MYB subgroup 14 (RAX / MYB37 clade)", "ABA / drought (MYB37 clade)", "26646286, 40743810, 35903240"),
        ("g27357", "thioredoxin-like 3-3 (TRL33)", "unresolved", "no abiotic evidence found"),
    ],
    "4_100931231-101192853": [
        ("g30974", "lambda-class GST (GSTL1-like)", "ROS homeostasis", "24399237, 23380449"),
    ],
    "4_107336728-107600612": [
        ("g31323,g31328,g31329", "B3, REM / VRN1 clade", "candidate (phenology context)", "12114624 (phenology only)"),
    ],
    "4_121924771-122186393": [
        ("g32205", "ClpB3 / APG6, chloroplast HSP100", "heat protection", "16995899, 16572840; conflict 17144892"),
        ("g32221", "SnRK1 catalytic alpha, KIN10-like", "energy-stress signalling", "17671505, 28740502"),
        ("g32199", "BAG co-chaperone (member unresolved)", "heat protection", "class only"),
    ],
    "5_5557888-5883272": [
        ("g32874", "CCD8 / MAX4, strigolactone biosynthesis", "strigolactone (drought / temperature)", "24379380, 25716094, 34719688; conflict 29538660"),
        ("g32870,g32871", "RCC1 family, not UVR8 (g32871 lacks the Trp triad)", "unresolved", "no claim"),
        ("g32865", "GATA transcription factor (class not resolved)", "not called", "not searched"),
    ],
    "5_111174564-111436186": [
        ("g38787", "P2B Ca2+-ATPase, ACA11-like", "salt / ion and Ca2+ signalling", "38923303; conflict 38799098"),
        ("g38779", "TIP aquaporin", "water transport", "28646139, 32544243; conflict 21397356"),
        ("g38778", "NOA1 / RIF1 plastid GTPase", "chilling / NO", "27891134, 17351048"),
        ("g38772", "mitochondrial small HSP", "heat protection", "14741379, 31877784"),
    ],
    "6_116698620-117089945": [
        ("g45753", "NDR1/HIN1-like (NHL) LEA_2", "salt / ABA (NHL)", "38005766, 41171317; conflict 26849212"),
        ("g45772,g45773", "flavonoid glucosyltransferase (UGT78-like)", "ROS homeostasis", "40617078, 11396938"),
        ("g45777", "CIPK (CIPK21 by name)", "salt / ion and Ca2+ signalling", "26198257, 35642834"),
    ],
}
SUPPORT = {  # overview "candidate gene support" column, one short phrase per region
    "1_8649433-8911055": "moderate (BURP RD22-like, YUC6, BOBBER1, WSD)",
    "3_147045758-147307380": "weak",
    "4_8145960-8407582": "weak",
    "4_11490368-11751990": "weak to moderate",
    "4_15212239-15556144": "moderate (MYB37-clade)",
    "4_100931231-101192853": "moderate (lambda GST)",
    "4_107336728-107600612": "none",
    "4_121924771-122186393": "moderate (ClpB3, SnRK1)",
    "5_5557888-5883272": "strong class (CCD8)",
    "5_111174564-111436186": "moderate (TIP, NOA1, ACA, mt-sHSP)",
    "6_116698620-117089945": "moderate (NHL LEA_2 at 1.7 kb, CIPK21)",
}
NEW_KEY = {       # the four new regions (Phase 7 literature pass)
    "2_145944878-146206511": [
        ("g12454", "LFG/GAAP-family TMBIM, GAAP4-like (not BI-1)", "drought / salt (GAAP4)", "40369998, 35034107"),
        ("g12464", "cytosolic MSRB5-type methionine sulfoxide reductase", "ROS homeostasis", "37035039"),
        ("g12461", "SAMDC4 / BUD2-type", "polyamine synthesis (development only)", "16699540; family 33138071"),
        ("g12462", "clade L PP2C (PP2C35 / PP2C19-like), not clade A", "light signalling (not ABA)", "41719245"),
    ],
    "2_176097944-176448788": [
        ("g15041,g15043", "Rx-type CNL (no characterised orthologue)", "biotic; abiotic class-level only", "37399697; conflict 26617631"),
        ("g15021,g15022,g15023,g15024", "Rx-type CNL array (g15023 a fragment)", "biotic; abiotic class-level only", "37399697, 38588338; conflict 33747005"),
        ("g15018", "Rx-type CNL, full length", "biotic; abiotic class-level only", "37399697, 38588338; conflict 26617631"),
        ("g15036", "CNL (R13L1 / DRL21-like LRRs)", "biotic; abiotic class-level only", "37399697; conflict 33747005, 37074757"),
    ],
    "4_101307505-101569127": [
        ("g31011", "HYH-like bZIP (GmSTF3/4 lineage)", "UV / light", "18055587, 37881032"),
    ],
    "5_8217606-8479228": [
        ("g33151,g33152", "mitochondrial ClpB4 (ClpB-m), not HSP101", "heat (expression only)", "17144892, 20141629, 42353767"),
    ],
}
NEW_SUPPORT = {
    "2_145944878-146206511": "moderate (GAAP4, MSRB5)",
    "2_176097944-176448788": "none (CNL cluster, class-level only)",
    "4_101307505-101569127": "weak (HYH-like, UV-B); same haplotype as region 8",
    "5_8217606-8479228": "weak (mitochondrial ClpB4)",
}
KEY.update(NEW_KEY)
SUPPORT.update(NEW_SUPPORT)


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def g2(x):
    return f"{x:.2g}".replace("e-0", "e-")


ev = {(r["region_id"], r["gene_id"].replace(".t1", "")): r for r in rows(os.path.join(T, "region_gene_evidence.tsv"))}
sig = collections.defaultdict(list)
for r in rows(os.path.join(T, "region_sig_snps.tsv")):
    sig[r["region_id"]].append(r)
ngenes = collections.Counter(r["region_id"] for r in rows(os.path.join(T, "region_gene_evidence.tsv")))
reg = rows(R)
assert len(reg) == 15
THR = {"EMMAX": 0.05 / 177259}


def rng(vals, fmt=lambda v: f"{v:,}"):
    lo, hi = min(vals), max(vals)
    return fmt(lo) if lo == hi else f"{fmt(lo)}-{fmt(hi)}"


def key_table(rid):
    out = ["| Gene | Identity | Function group | Nearest sig. SNP (bp) | SNPs body / promoter | Best p body / promoter (EMMAX) | Literature |",
           "|---|---|---|---|---|---|---|"]
    for genes, ident, fg, litc in KEY[rid]:
        gs = genes.split(",")
        rs = [ev[(rid, g)] for g in gs]
        dist = rng([int(r["dist_to_sig_snp"]) for r in rs])
        nb = rng([int(r["n_snps_body"]) for r in rs], str)
        npr = rng([int(r["n_snps_promoter"]) for r in rs], str)
        best = [(float(r[k]), w) for r in rs for k, w in (("best_p_EMMAX_body", "body"), ("best_p_EMMAX_promoter", "promoter")) if r[k]]
        if best:
            v, w = min(best)
            bp = g2(v) + (f" ({w}{', significant' if v <= THR['EMMAX'] else ''})" if v < 1e-3 else "")
        else:
            bp = "—"
        out.append(f"| {', '.join(gs)} | {ident} | {fg} | {dist} | {nb} / {npr} | {bp} | {litc} |")
    return "\n".join(out)


def overview():
    out = ["| # | Region | kb | Genes | Sig. SNPs | Best p (EMMAX / RDA) | Predictors | Candidate gene support |",
           "|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(reg, 1):
        rid = r["region_id"]
        xs = sig[rid]
        preds = [t for t in r["traits"].split(",") if t != "climate_multivariate"]
        pr = ", ".join(sorted(preds, key=lambda t: int(t.split("_")[1]))) if len(preds) <= 4 else f"{len(preds)} bio"
        kb = round((int(r["end"]) - int(r["start"])) / 1000)
        out.append(f"| {i} | {rid} | {kb} | {ngenes[rid]} | {len(xs)} | "
                   f"{g2(min(float(x['EMMAX_minp']) for x in xs))} / {g2(min(float(x['RDA_p']) for x in xs))} | {pr} | "
                   f"{SUPPORT.get(rid, 'TODO')} |")
    return "\n".join(out)


cur = rows(os.path.join(T, "curation.tsv"))
newc = [r for r in cur if r["source"] != "carried"]
nr = {r["gene_id"] for r in rows(os.path.join(T, "neighbour_read.tsv"))}
lit = rows(os.path.join(T, "literature.tsv"))
raw = rows(os.path.join(T, "literature_raw.tsv"))
dropped = rows(os.path.join(B, "logs/lit_dropped.tsv"))
guess = rows(os.path.join(T, "neighbourhood_guess.tsv"))
COUNTS = {
    "N_FLAGGED": sum(r["source"] == "flagged" for r in newc),
    "N_NEIGHBOUR": sum(r["source"] == "neighbour" for r in newc),
    "N_ARRAY": sum(r["source"] == "array" for r in newc),
    "N_READ_NEW": sum(r["gene_id"] in nr for r in newc),
    "N_LIT_NEW": len(rows(os.path.join(B, "logs/agent_batch_v2_identity.tsv"))),  # the first agent's genes
    "N_GUESS": len(guess), "N_GUESS_A": sum(r["priority"] == "A" for r in guess),
    "N_BIO2": sum(1 for x in rows(os.path.join(T, "region_sig_snps.tsv")) if x["EMMAX_pred"] == "bio_2"),
}
assert COUNTS["N_FLAGGED"] + COUNTS["N_NEIGHBOUR"] + COUNTS["N_ARRAY"] + COUNTS["N_READ_NEW"] == 83, COUNTS
v1raw = {(r["gene_id"], r["pmid"], r["role"]) for r in rows(os.path.join(B, "../../Trifolium86MainK5/regionsER/tables/literature_raw.tsv"))}
n_new_rows = sum(1 for r in raw if (r["gene_id"], r["pmid"], r["role"]) not in v1raw)
COUNTS["LIT_SUMMARY"] = (
    f"{len(raw)} rows and {len({r['pmid'] for r in raw})} PMIDs ({len(raw) - n_new_rows} rows carried from the PC3 pass for "
    f"genes still in a region, {n_new_rows} new), all re-verified against NCBI esummary (first author and year): "
    f"{len(lit)} kept, {len(dropped)} dropped; {sum(r['role'] == 'rejected' for r in lit)} kept rows are `rejected`.")
src = open(os.path.join(B, "REGIONS_v2.src.md")).read()
for k, v in COUNTS.items():
    src = src.replace("{{" + k + "}}", str(v))
src = src.replace("{{OVERVIEW}}", overview())
src = re.sub(r"\{\{KEY:(\S+?)\}\}", lambda m: key_table(m.group(1)), src)
assert "{{" not in src, re.findall(r"\{\{[^}]*\}\}", src)
with open(os.path.join(B, "REGIONS_v2.md"), "w") as fh:
    fh.write(src)
print(f"INFO: REGIONS_v2.md written ({src.count(chr(10))} lines)")
