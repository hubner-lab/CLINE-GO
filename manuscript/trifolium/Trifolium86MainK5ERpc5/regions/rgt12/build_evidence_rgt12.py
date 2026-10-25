#!/usr/bin/env python3
"""Per-gene evidence table for Fig 3 region rgt_12 (5_5557888-5883272).

Plan ~/.claude/plans/idempotent-conjuring-noodle.md, Phase 4. Read-only on every v2 table: region_gene_evidence.tsv,
literature.tsv and REGIONS_v2.md are what report_check.py validates, so the new literature pass lives in this
directory (rgt12/tables/literature_raw.tsv -> lit_check.py rgt12 -> rgt12/tables/literature.tsv).

Stage `base`  : rgt12/tables/evidence_base.tsv (31 genes) + rgt12/logs/lit_batch_{A,B}.txt (gene lines for the agents)
Stage `merge` : base + rgt12/tables/literature.tsv (verified new rows) + rgt12/logs/agent_*_identity.tsv
                -> regions/evidence_rgt12.tsv
Gene set      = the pipeline's genes_per_region.tsv for the region (asserted equal as a set to the v2 evidence rows).
Predictors    = per significant SNP from the pipeline's selected_snps_per_trait.tsv (the significant (SNP, trait) pairs),
                ordered by p. The `pred_EMMAX_*` columns of gene_snps.tsv are the best EMMAX predictor of ANY tested SNP
                in the gene, significant or not, and are labelled as such.
SNP counts    = the tested set of gene_snps.tsv (imputed-MAF-0.1 matrix used by EMMAX/RDA), body and 5 kb promoter,
                strand-aware (evidence_inputs.py:224).
Usage (from the repo root): python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/rgt12/build_evidence_rgt12.py base|merge
"""
import csv
import os
import sys

ROOT = os.getcwd()
B = os.path.dirname(os.path.abspath(__file__))
REG = os.path.dirname(B)
T2 = os.path.join(REG, "tables")
GEA = os.path.join(ROOT, "Trifolium86MainK5ERpc5_results/GEA/tables")
R = "5_5557888-5883272"
PROM = 5000
# plot labels drawn in rgt_12 (prep_regionplots.py LABELS, region 9 block)
LABELS = {"g32874": "CCD8", "g32870": "RCC1", "g32871": "RCC1", "g32865": "GATA",
          "g32863": "SKM1", "g32856": "MLK3", "g32864": "ATHB8"}  # last three added 2026-10-17 (user)
# half LD window: config_Trifolium86MainK5ERpc5.yaml snp_clumping_distance (genome-wide r2 0.2 distance); regions are the
# significant-SNP cluster +- this distance. User 2026-10-17: rank genes only by inside / outside it, not by distance.
HALF_LD = 130811
# genes not sent to the literature pass, with the reason (stays in the table)
NOT_SEARCHED = {
    "g32849": "class V GH18 chitinase array, literature carried from the v2 defence-gene pass",
    "g32850": "class V GH18 chitinase array, literature carried from the v2 defence-gene pass",
    "g32851": "class V GH18 chitinase array, literature carried from the v2 defence-gene pass",
    "g32852": "class V GH18 chitinase array, literature carried from the v2 defence-gene pass",
    "g32871": "RCC1, not UVR8 (Trp triad absent), identity carried from the PC3 pass",
    "g32874": "CCD8 / MAX4, literature carried from the v2 pass",
    "g32877": "plastid psaA/psaB copy in the nuclear assembly (NUPT-like insertion); no nuclear gene function",
    "g32878": "plastid psaA copy in the nuclear assembly (NUPT-like insertion); no nuclear gene function",
    "g32879": "plastid psbD copy in the nuclear assembly (NUPT-like insertion); no nuclear gene function",
}
BATCH_A = ["g32865", "g32872", "g32869", "g32864", "g32873", "g32863", "g32870", "g32856", "g32855", "g32854",
           "g32867"]


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(p, fields, data):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(data)


def short(g):
    return g.split(".")[0]


def base():
    ev = [r for r in rows(os.path.join(T2, "region_gene_evidence.tsv")) if r["region_id"] == R]
    gs = {r["gene_id"]: r for r in rows(os.path.join(T2, "gene_snps.tsv")) if r["region_id"] == R}
    pipe = {r["gene_id"]: r for r in rows(os.path.join(GEA, "genes_per_region.tsv"))
            if r["region_id"].startswith(R + "_")}
    if {r["gene_id"] for r in ev} != set(pipe) or set(gs) != set(pipe):
        sys.exit(f"ERROR: gene sets differ: evidence {len(ev)}, gene_snps {len(gs)}, pipeline {len(pipe)}")
    sig = {r["snp_id"]: r for r in rows(os.path.join(T2, "region_sig_snps.tsv")) if r["region_id"] == R}
    preds = {}
    for r in rows(os.path.join(GEA, "selected_snps_per_trait.tsv")):
        if r["SNPID"] in sig:
            preds.setdefault(r["SNPID"], []).append((float(r["min_pvalue"]), r["trait"]))
    if set(preds) != set(sig):
        sys.exit(f"ERROR: significant SNPs without a trait row: {sorted(set(sig) - set(preds))}")
    print(f"INFO: {len(ev)} genes (set = genes_per_region.tsv), {len(sig)} significant SNPs")

    out = []
    for e in sorted(ev, key=lambda r: int(r["start"])):
        g, p, s = e["gene_id"], pipe[e["gene_id"]], gs[e["gene_id"]]
        lo, hi, strand = int(e["start"]), int(e["end"]), e["strand"]
        prom = (hi, hi + PROM) if strand == "-" else (lo - PROM, lo)
        dist = {k: (0 if lo <= int(w["pos"]) <= hi else min(abs(int(w["pos"]) - lo), abs(int(w["pos"]) - hi)))
                for k, w in sig.items()}
        d = min(dist.values())
        sid, *ties = sorted(k for k in sig if dist[k] == d)
        pos = int(sig[sid]["pos"])
        if lo <= pos <= hi:
            where = "body"
        elif prom[0] <= pos <= prom[1]:
            where = "promoter (5 kb)"
        else:
            up = (pos < lo) if strand == "+" else (pos > hi)
            where = ("upstream" if up else "3' flank (< 5 kb)" if d <= PROM else "downstream")
        # the v2 table's nearest SNP is one of possibly several tied SNPs; the distance must agree
        if int(e["dist_to_sig_snp"]) != d:
            sys.exit(f"ERROR: {g}: distance {d} disagrees with region_gene_evidence {e['dist_to_sig_snp']}")
        v = sig[sid]
        pr = sorted(preds[sid])
        both = min((k for k in sig if sig[k]["sig_methods"] == "EMMAX,RDA"), key=lambda k: (dist[k], k))
        out.append({
            "gene_id": short(g), "transcript": g, "start": lo, "end": hi, "strand": strand,
            "annotation": e["desc"], "pfam": e["pfam"], "kegg": e["ko_names"], "eggnog_best_hit": p["em_target"],
            "array_id": e["array_id"], "array_size": e["array_size"], "plot_label": LABELS.get(short(g), ""),
            "v2_call": e["call"], "v2_category": e["category"], "v2_function_group": e["function_group"],
            "v2_identity": e["identity"],
            "n_snps_body": s["n_snps_body"], "best_p_EMMAX_body_any": s["best_p_EMMAX_body"],
            "pred_EMMAX_body_any": s["pred_EMMAX_body"], "best_p_RDA_body_any": s["best_p_RDA_body"],
            "n_snps_promoter": s["n_snps_promoter"], "best_p_EMMAX_promoter_any": s["best_p_EMMAX_promoter"],
            "pred_EMMAX_promoter_any": s["pred_EMMAX_promoter"], "best_p_RDA_promoter_any": s["best_p_RDA_promoter"],
            "sig_snps_in_gene": ";".join(x for x in (s["sig_snps_body"], s["sig_snps_promoter"]) if x),
            "nearest_sig_snp": sid + ("," + ",".join(ties) if ties else ""), "dist_to_sig_snp": d,
            "nearest_sig_position": where, "nearest_sig_methods": v["sig_methods"],
            "nearest_sig_p_EMMAX": v["EMMAX_minp"] if "EMMAX" in v["sig_methods"] else "",
            "nearest_sig_p_RDA": v["RDA_p"],
            "nearest_sig_predictors": ",".join(t for _, t in pr), "nearest_sig_n_bio": sum(t.startswith("bio_") for _, t in pr),
            "nearest_both_methods_snp": both, "dist_to_both_methods_snp": dist[both],
            "both_methods_snp_predictors": ",".join(t for _, t in sorted(preds[both])),
            "v2_pmids_key": e["pmids_key"], "v2_pmids_support": e["pmids_support"],
            "v2_pmids_conflict": e["pmids_conflict"], "not_searched_reason": NOT_SEARCHED.get(short(g), ""),
        })
    fields = list(out[0])
    write(os.path.join(B, "tables/evidence_base.tsv"), fields, out)

    for name, keep in (("A", lambda r: r["gene_id"] in BATCH_A),
                       ("B", lambda r: r["gene_id"] not in BATCH_A and not r["not_searched_reason"])):
        lines = ["GENES:"]
        for r in (r for r in out if keep(r)):
            snp = (f"significant SNPs IN this gene: {r['sig_snps_in_gene']}" if r["sig_snps_in_gene"] else
                   f"nearest significant SNP {r['nearest_sig_snp']} at {r['dist_to_sig_snp']} bp ({r['nearest_sig_position']})")
            lines.append(f"- {r['transcript']} | region {R} | {r['annotation']} | PFAM={r['pfam'] or '-'} | "
                         f"KEGG={r['kegg'] or '-'} | eggNOG best hit={r['eggnog_best_hit']} | "
                         f"v2 call={r['v2_call'] or '-'}/{r['v2_category'] or '-'} | {snp}; "
                         f"predictors of that SNP: {r['nearest_sig_predictors']}")
        with open(os.path.join(B, f"logs/lit_batch_{name}.txt"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"INFO: batch {name}: {len(lines) - 1} genes")


LIT_COLS = ["gene_id", "role", "pmid", "first_author", "year", "species", "legume", "evidence", "direction", "finding"]
ID_COLS = ["gene_id", "identity", "evidence_verdict", "note"]
VOCAB = {"role": {"key", "support", "conflict", "rejected"}, "legume": {"yes", "no"},
         "evidence": {"mutant", "transgenic", "natural_variation", "expression", "review"},
         "direction": {"improves_tolerance", "reduces_tolerance", "flowering", "expression_only", "none"}}


def raw():
    """agent_batch_{A,B}.tsv -> tables/literature_raw.tsv (the lit_check.py input), after a schema check."""
    base_ids = {r["gene_id"] for r in rows(os.path.join(B, "tables/evidence_base.tsv"))}
    out, ids = [], []
    for name in ("A", "B"):
        lit = rows(os.path.join(B, f"logs/agent_batch_{name}.tsv"))
        idt = rows(os.path.join(B, f"logs/agent_batch_{name}_identity.tsv"))
        want = {line.split(" | ")[0][2:].split(".")[0] for line in open(os.path.join(B, f"logs/lit_batch_{name}.txt"))
                if line.startswith("- ")}
        if {r["gene_id"] for r in idt} != want:
            sys.exit(f"ERROR: batch {name}: identity rows {sorted(r['gene_id'] for r in idt)} != batch {sorted(want)}")
        for r in lit:
            if list(r) != LIT_COLS or r["gene_id"] not in want:
                sys.exit(f"ERROR: batch {name}: bad literature row {r}")
            bad = [k for k, v in VOCAB.items() if r[k] not in v]
            if bad:
                sys.exit(f"ERROR: batch {name}: {r['gene_id']} {r['pmid']}: values outside the vocabulary in {bad}")
            out.append({"region_id": R, **r})
        if any(list(r) != ID_COLS or r["evidence_verdict"] not in ("strong", "moderate", "weak", "none") for r in idt):
            sys.exit(f"ERROR: batch {name}: bad identity row")
        ids += idt
    assert {r["gene_id"] for r in ids} <= base_ids
    write(os.path.join(B, "tables/literature_raw.tsv"), ["region_id"] + LIT_COLS, out)
    write(os.path.join(B, "tables/identity.tsv"), ID_COLS, ids)
    print(f"INFO: literature_raw.tsv: {len(out)} rows, {len({r['pmid'] for r in out})} PMIDs; identity.tsv: {len(ids)} genes")


# function group per gene, assigned from the resolved identity (tables/identity.tsv) or, for genes not searched here,
# from the v2 row; one place, so every group is reviewable
GROUP = {
    "g32849": "defence (class V chitinase; ABA / osmotic induction is expression only)",
    "g32850": "defence (class V chitinase; ABA / osmotic induction is expression only)",
    "g32851": "defence (class V chitinase; ABA / osmotic induction is expression only)",
    "g32852": "defence (class V chitinase; ABA / osmotic induction is expression only)",
    "g32853": "housekeeping (GPI-anchor synthesis)",
    "g32854": "chromatin (RdDM, SHH1)",
    "g32855": "transcription factor (IDD, member unresolved)",
    "g32856": "phenology (flowering repressor, FLC)",
    "g32858": "biotic response (ACRE-like)",
    "g32859": "unknown (LFS-like Bet v1 fold)",
    "g32860": "unknown (LFS-like Bet v1 fold)",
    "g32861": "lipid metabolism (MAG lipase)",
    "g32862": "unknown (DUF1677)",
    "g32863": "heat tolerance of reproduction (pollen tube growth at 30 C)",
    "g32864": "vascular development (xylem differentiation)",
    "g32865": "transcription factor; seed dormancy (paralogue evidence)",
    "g32866": "mitochondrial / biotic (NOXY2)",
    "g32867": "housekeeping (peptidyl-tRNA hydrolase)",
    "g32869": "unknown (RING-H2 + SWIM)",
    "g32870": "unresolved (RCC1, At5g11580-like)",
    "g32871": "unresolved (RCC1, At5g11580-like)",
    "g32872": "unknown (DUF810)",
    "g32873": "housekeeping (YebC / TACO1)",
    "g32874": "strigolactone (drought / temperature)",
    "g32875": "unknown",
    "g32877": "plastid insertion (NUPT)",
    "g32878": "plastid insertion (NUPT)",
    "g32879": "plastid insertion (NUPT)",
    "g32881": "unknown (CCHC zinc knuckle; possibly TE-derived)",
    "g32883": "cell-wall glycosyltransferase (GT47, uncharacterized clade)",
    "g32884": "cell-wall glycosyltransferase (GT47, uncharacterized clade)",
}
# genes not searched here: verdict from the v2 report (REGIONS_v2.md section 12 and the defence side line)
VERDICT_V2 = {"g32874": "strong", "g32849": "weak", "g32850": "weak", "g32851": "weak", "g32852": "weak",
              "g32871": "none", "g32877": "none", "g32878": "none", "g32879": "none"}
# identity amendments from this pass for genes that were not re-searched (batch A, g32870 note)
ID_AMEND = {"g32871": "RCC1, At5g11580-like (PTHR22870:SF155), NOT UVR8; At5g11580 aligns to g32870 at residues 12-89 "
                      "and to g32871 at 143-549, same strand, 84 bp apart: probably one gene split into two models"}
TIER_EV = {"strong": 3, "moderate": 2, "weak": 1, "none": 0}


def proximity(r):
    return 1 if int(r["dist_to_sig_snp"]) <= HALF_LD else 0


def merge():
    base_rows = rows(os.path.join(B, "tables/evidence_base.tsv"))
    base_keys = list(base_rows[0])  # before the rows are updated in place
    lit = rows(os.path.join(B, "tables/literature.tsv"))
    raw_n = len(rows(os.path.join(B, "tables/literature_raw.tsv")))
    if len(lit) != raw_n:
        sys.exit(f"ERROR: lit_check dropped {raw_n - len(lit)} rows; resolve logs/lit_dropped.tsv first")
    idt = {r["gene_id"]: r for r in rows(os.path.join(B, "tables/identity.tsv"))}
    if set(GROUP) != {r["gene_id"] for r in base_rows} or set(idt) | set(VERDICT_V2) != set(GROUP):
        sys.exit("ERROR: GROUP / identity / VERDICT_V2 do not cover the 31 genes exactly")
    out = []
    for r in base_rows:
        g = r["gene_id"]
        mine = [x for x in lit if x["gene_id"] == g]
        cite = lambda role: ",".join(f"{x['pmid']}({x['first_author']} {x['year']})" for x in mine if x["role"] == role)
        verdict = idt[g]["evidence_verdict"] if g in idt else VERDICT_V2[g]
        lead = next((x for x in mine if x["role"] == "key"), next((x for x in mine if x["role"] == "support"), None))
        r.update({
            "identity": idt[g]["identity"] if g in idt else ID_AMEND.get(g, r["v2_identity"]),
            "identity_check": idt[g]["note"] if g in idt else "carried from v2",
            "function_group": GROUP[g], "evidence_verdict": verdict,
            "new_pmids_key": cite("key"), "new_pmids_support": cite("support"),
            "new_pmids_conflict": cite("conflict"), "new_pmids_rejected": cite("rejected"),
            "n_new_refs_cited": sum(x["role"] != "rejected" for x in mine),
            "n_new_legume_refs_cited": sum(x["role"] != "rejected" and x["legume"] == "yes" for x in mine),
            "new_lead_finding": lead["finding"] if lead else "",
            "tier_evidence": TIER_EV[verdict], "in_half_ld_window": proximity(r),
        })
        out.append(r)
    out.sort(key=lambda r: (-r["tier_evidence"], -r["in_half_ld_window"], int(r["start"])))
    for i, r in enumerate(out, 1):
        r["interest_rank"] = i
        r["unlabelled_candidate"] = "TRUE" if (r["tier_evidence"] >= 1 and not r["plot_label"] and
                                               not r["function_group"].startswith(("defence", "biotic"))) else ""
    out.sort(key=lambda r: int(r["start"]))
    fields = ["interest_rank", "gene_id", "plot_label", "unlabelled_candidate", "identity", "function_group",
              "evidence_verdict", "tier_evidence", "in_half_ld_window"] + \
             [k for k in base_keys if k not in ("gene_id", "plot_label")] + \
             ["identity_check", "new_pmids_key", "new_pmids_support", "new_pmids_conflict", "new_pmids_rejected",
              "n_new_refs_cited", "n_new_legume_refs_cited", "new_lead_finding"]
    assert len(fields) == len(set(fields)), "duplicate columns"
    write(os.path.join(REG, "evidence_rgt12.tsv"), fields, out)
    print(f"INFO: evidence_rgt12.tsv: {len(out)} genes")
    for r in sorted(out, key=lambda r: r["interest_rank"])[:10]:
        print(f"  {r['interest_rank']:>2} {r['gene_id']} [{r['plot_label'] or '-'}] {r['evidence_verdict']}/"
              f"{r['in_half_ld_window']}: {r['identity'][:60]}")


if __name__ == "__main__":
    {"base": base, "raw": raw, "merge": merge}[sys.argv[1]]()
