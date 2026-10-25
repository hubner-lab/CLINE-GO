#!/usr/bin/env python3
"""Hand calls for the 353 genes of the 15 EMMAX+RDA regions of Trifolium86MainK5ERpc5 (5 PCs)
(plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7). Adapted from ANALYSIS/Trifolium86MainK5/regionsER/
write_curation.py (the PC3 run, 315 genes).

Rules: exactly the curation rules in ANALYSIS/Trifolium86Main/genes3of3/gene_sheet.R header (calls stress / phenology /
reject / unclear; seven ranked categories; G1-G3, G2b; R1-R5; tandem copies take the array's call).
  source = carried    the gene is one of the 315 genes of the PC3 run (same GFF, same gene id): its PC3 call, category,
                      reason and lit flag are kept unchanged, and so is its PC3 "read, not added" reason when it had no
                      call. v1_source keeps the PC3 source (carried from genes3of3 / flagged / neighbour / array).
           flagged    flagged by the screen, never curated: called here (83 new genes: the 82 of the four new regions
                      plus g27344, which entered 4_15.2M when its left edge moved 5.2 kb)
           neighbour  unflagged, read here, added because the annotation names a stress / phenology class the keyword
                      and GO lines missed. Every unflagged new gene read and NOT added is in tables/neighbour_read.tsv.
           array      an unflagged tandem copy of a called gene in the same region
  lit = TRUE          the gene goes to the literature (PC3 genes: their PC3 literature rows are carried; new genes: the
                      Phase 7 literature agents): every stress and phenology call, and the `unclear` genes whose
                      annotation names a stress-linked class.

Outputs: tables/curation.tsv (one row per region gene), tables/neighbour_read.tsv
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/write_curation.py   (from anywhere)
"""
import csv
import os

B = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(B, "tables")
V1 = os.path.join(B, "../../Trifolium86MainK5/regionsER/tables")

# ---- flagged, never curated ------------------------------------------------------------------------------
FLAGGED = {
    # 2_145944878-146206511
    "g12446.t1": ("reject", "", "G2b: plastid glutamine synthetase GS2 (N assimilation / photorespiration); GO salt term transferred from ortholog"),
    "g12450.t1": ("reject", "", "G1: broad GO only; FAR1-RELATED SEQUENCE (MULE transposase-derived family), member not FHY3 / FAR1"),
    "g12462.t1": ("unclear", "aba", "R2: PP2C with no ko, clade unresolved"),
    "g12464.t1": ("unclear", "oxidative", "candidate: methionine sulfoxide reductase B (msrB), repair of oxidised methionine; GO cold / oxidative transferred; class not in list"),
    "g12470.t1": ("reject", "", "lexical: 'proline' hit the proline-rich nuclear receptor coactivator (PNRC) motif"),
    # 2_176097944-176448788
    "g15018.t1": ("reject", "", "biotic: NB-LRR (RPM1 / RPS3 ko)"),
    "g15021.t1": ("reject", "", "biotic: NB-LRR (RPM1 / RPS3 ko), 4-copy array"),
    "g15022.t1": ("reject", "", "biotic: NB-LRR (RPM1 / RPS3 ko), 4-copy array"),
    "g15023.t1": ("reject", "", "biotic: NB-LRR (RPM1 / RPS3 ko), 4-copy array"),
    "g15029.t1": ("reject", "", "G1: broad GO only (ethylene); XLG3 extra-large G protein (root growth / immunity)"),
    "g15036.t1": ("reject", "", "biotic: NB-LRR (SHOC2 ko)"),
    "g15048.t1": ("reject", "", "G1: broad GO only; ABCB auxin / MDR transporter (same call as g27361 in the PC3 run)"),
    # 4_101307505-101569127
    "g30998.t1": ("reject", "", "G1: broad GO only; TPR protein (same call as g30978 in the PC3 run)"),
    # user 2026-10-16: UV-B photoreception is called `unclear` / uv_light as in the PC3 run (was phenology here)
    "g31011.t1": ("unclear", "uv_light", "candidate: HY5-like bZIP (ko HY5; HYH-like by blastp), UVR8 UV-B signalling; UV-B is abiotic but not a ranked category"),
    # 5_8217606-8479228
    "g33148.t1": ("reject", "", "ARR-B type-B response regulator (K14491): cytokinin signalling; keyword hit is lexical (same call as g17901 in genes3of3); GO water deprivation transferred"),
    "g33151.t1": ("stress", "heat", "ClpB (K03695), HSP100 class (same call as g32205)"),
    "g33152.t1": ("stress", "heat", "ClpB (K03695), HSP100 class, tandem copy of g33151"),
    "g33153.t1": ("reject", "", "G1: broad GO only; U-box / ARM-repeat protein (VAC8 ko)"),
    "g33154.t1": ("reject", "", "G1: broad GO only; SDF2-like ER protein (ER stress / immunity)"),
    "g33160.t1": ("unclear", "", "bare WRKY (TF rule)"),
}

# ---- neighbour read: unflagged new genes added -----------------------------------------------------------
NEIGHBOUR = {
    "g12454.t1": ("unclear", "", "candidate: BI1 / TMBIM family; Arabidopsis BI-1 attenuates heat- and ER-stress cell death; member unresolved"),
    "g12461.t1": ("unclear", "drought_osmotic", "candidate: S-adenosylmethionine decarboxylase (SAMDC), polyamine synthesis; SAMDC overexpression reported to raise salt / drought tolerance"),
}

# ---- unflagged new genes read with a specific reason for not adding them (rest get the generic reason) ---
NOT_ADDED = {
    "g12445.t1": "elongation factor Tu (TUFM ko): translation (housekeeping)",
    "g12452.t1": "FGAM synthase (purL): purine synthesis (housekeeping)",
    "g12453.t1": "FGAM synthase (purL): purine synthesis (housekeeping)",
    "g12459.t1": "MO25 / CAB39 scaffold of STE20 kinases: holds the region's 2 significant promoter SNPs; no named abiotic role",
    "g12460.t1": "F-box / LRR protein: holds the region's 2 significant promoter SNPs; member unknown",
    "g12463.t1": "SAC3 / GANP mRNA export factor (housekeeping)",
    "g12467.t1": "RING E3 ligase (zf-C3HC4): large family, member unknown",
    "g12469.t1": "DUF679 / DMP membrane protein: development / reproduction",
    "g15020.t1": "DEAD-box helicase DBP3 ko: ribosome biogenesis; holds 2 of the region's significant promoter SNPs",
    "g15039.t1": "ABCC1-like MRP transporter: vacuolar glutathione-conjugate / metal transport, not a ranked category",
    "g15041.t1": "biotic: NB-LRR (SHOC2 ko)", "g15043.t1": "biotic: NB-LRR (SHOC2 ko)",
    "g15046.t1": "biotic: NB-LRR (SHOC2 ko)",
    "g15019.t1": "PGG ankyrin-repeat protein: member unknown", "g15037.t1": "PGG ankyrin-repeat protein: member unknown",
    "g15038.t1": "PGG ankyrin-repeat protein: member unknown", "g15040.t1": "PGG ankyrin-repeat protein: member unknown",
    "g15042.t1": "PGG ankyrin-repeat protein: member unknown", "g15044.t1": "PGG ankyrin-repeat protein: member unknown",
    "g15045.t1": "PGG ankyrin-repeat protein: member unknown",
    "g15047.t1": "copia LTR retrotransposon", "g31014.t1": "LTR retrotransposon gag / protease",
    "g27344.t1": "pectate lyase: cell-wall remodelling (new at 5 PCs: the 4_15.2M left edge moved 5.2 kb)",
    "g30997.t1": "MULE / SWIM transposase-derived protein", "g31000.t1": "FAR1-RELATED SEQUENCE 5-like (MULE transposase-derived family)",
    "g31002.t1": "annexin D5-like (ANN5 is pollen-specific); annexin stress links family-level only",
    "g31005.t1": "no annotation; holds the region's significant promoter SNP 4:101438316",
    "g31007.t1": "BET bromodomain protein (GTE family): chromatin", "g31008.t1": "BET bromodomain protein (GTE family): chromatin",
    "g31009.t1": "BET bromodomain protein (GTE family): chromatin",
    "g33149.t1": "CHD3 chromatin remodeller (PICKLE-like): development; carries a promoter SNP at EMMAX p 6.0e-5",
    "g33150.t1": "CCR4-NOT deadenylase CAF1: mRNA turnover; stress links family-level only",
    "g33161.t1": "POLAR stomatal-lineage polarity protein: stomatal development; water-use link indirect",
    "g33162.t1": "chorismate mutase CM2: aromatic amino-acid synthesis (primary metabolism)",
}
GENERIC = "housekeeping / development / metabolism; the annotation names no abiotic-stress or phenology class"
# defence-class genes literature-tested for abiotic roles (user request 2026-10-16); calls unchanged, lit set TRUE,
# including 6 carried genes (g38791, g38792, g32849-g32852, g26559) whose PC3 lit flag was FALSE
LIT_DEFENCE = {g + ".t1" for g in ("g15018", "g15021", "g15022", "g15023", "g15024", "g15036", "g15041", "g15043",
                                   "g15046", "g38791", "g38792", "g32849", "g32850", "g32851", "g32852", "g26559")}
LIT_NEW = {"g12462.t1", "g12464.t1", "g12454.t1", "g12461.t1", "g31011.t1", "g33151.t1", "g33152.t1"}


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


sh = rows(os.path.join(OUT, "gene_sheet.tsv"))
assert len(sh) == 353
by_id = {r["gene_id"]: r for r in sh}
v1 = {r["gene_id"]: r for r in rows(os.path.join(V1, "curation.tsv"))}
v1_nr = {r["gene_id"]: r["not_added_because"] for r in rows(os.path.join(V1, "neighbour_read.tsv"))}
assert len(v1) == 315
new = [g for g in by_id if g not in v1]
assert len(new) == 83, len(new)
for d in (FLAGGED, NEIGHBOUR, NOT_ADDED):
    assert set(d) <= set(new), set(d) - set(new)
assert not set(FLAGGED) & set(NEIGHBOUR) and not set(NEIGHBOUR) & set(NOT_ADDED) and not set(FLAGGED) & set(NOT_ADDED)
assert LIT_NEW <= set(FLAGGED) | set(NEIGHBOUR)

cur, nbr = [], []
for r in sh:
    g = r["gene_id"]
    if g in v1:
        o = v1[g]
        call, cat, reason, src, lit, v1src = o["call"], o["category"], o["reason"], "carried", o["lit"] == "TRUE", o["source"]
        if not call:
            assert g in v1_nr, g
            nbr.append([r["region_id"], g, r["desc"], v1_nr[g]])
    elif r["flagged"] == "TRUE":
        assert g in FLAGGED, f"flagged gene without a call: {g}"
        (call, cat, reason), src, v1src = FLAGGED[g], "flagged", ""
        lit = g in LIT_NEW or call in ("stress", "phenology")
    elif g in NEIGHBOUR:
        (call, cat, reason), src, v1src = NEIGHBOUR[g], "neighbour", ""
        lit = g in LIT_NEW or call in ("stress", "phenology")
    else:
        call, cat, reason, src, lit, v1src = "", "", "", "", False, ""
        nbr.append([r["region_id"], g, r["desc"], NOT_ADDED.get(g, GENERIC)])
    lit = lit or g in LIT_DEFENCE
    cur.append([r["region_id"], g, call, cat, reason, src, "TRUE" if lit else "FALSE", v1src])

# an unflagged tandem copy of a called gene takes the array's call (source = array) and leaves the neighbour list
called = {}
for (rid, g, call, cat, reason, *_), r in zip(cur, sh):
    if call and int(r["array_size"]) > 1:
        called.setdefault((rid, r["array_id"]), (call, cat, f"tandem copy of {g}: {reason}"))
inherited = set()
for x, r in zip(cur, sh):
    k = (x[0], r["array_id"])
    if not x[2] and k in called:
        assert x[5] != "carried", x  # the PC3 run already applied this rule to its genes
        x[2], x[3], x[4] = called[k]
        x[5] = "array"
        inherited.add(x[1])
nbr = [n for n in nbr if n[1] not in inherited]
arr = {}
for (rid, g, call, cat, *_), r in zip(cur, sh):
    if int(r["array_size"]) > 1:
        arr.setdefault((rid, r["array_id"]), set()).add((call, cat))
bad = {k: v for k, v in arr.items() if len(v) > 1}
assert not bad, bad

with open(os.path.join(OUT, "curation.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["region_id", "gene_id", "call", "category", "reason", "source", "lit", "v1_source"])
    w.writerows(cur)
with open(os.path.join(OUT, "neighbour_read.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["region_id", "gene_id", "desc", "not_added_because"])
    w.writerows(nbr)
n = lambda f: sum(1 for x in cur if f(x))
nn = lambda f: sum(1 for x in cur if x[5] != "carried" and f(x))
print(f"INFO: {len(cur)} genes: carried {n(lambda x: x[5] == 'carried')}, flagged {n(lambda x: x[5] == 'flagged')}, "
      f"neighbour {n(lambda x: x[5] == 'neighbour')}, array {n(lambda x: x[5] == 'array')}, read-not-added "
      f"{len(nbr)} ({sum(1 for x in nbr if x[1] not in v1)} new); stress {n(lambda x: x[2] == 'stress')}, "
      f"phenology {n(lambda x: x[2] == 'phenology')}, unclear {n(lambda x: x[2] == 'unclear')}; to literature "
      f"{n(lambda x: x[6] == 'TRUE')} (new genes {nn(lambda x: x[6] == 'TRUE')})")
