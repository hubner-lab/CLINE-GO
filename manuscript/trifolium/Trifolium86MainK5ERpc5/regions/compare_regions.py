#!/usr/bin/env python3
"""PC5 regions vs the PC3 regions already described (plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 2).

PC3 = Trifolium86MainK5ER (EMMAX n_pcs 3, RDA condition_pcs 3; 13 regions described in
ANALYSIS/Trifolium86MainK5/regionsER/REGIONS_ER.md). PC5 = Trifolium86MainK5ERpc5 (5 / 5), same imputation, same rule
(pipeline Overlap, EMMAX bonf 0.05 x RDA qval 0.05 within 130,811 bp). A PC5 region is `old` when its span overlaps a
PC3 region span on the same chromosome (spans already carry the +-130,811 bp clumping pad), else `new`.

Stress: stress_screen.R run on the PC5 regions (stress/), lines ks_any (keyword strict) and go_any (GO_STRESS) in
tandem-array units; set level as in ANALYSIS/Trifolium86MainK5/compare_sets.py (observed arrays, expected = sum of
n_arrays x genome array rate, one-sided hypergeometric p of the pooled arrays).

Outputs (this dir): regions_v2.tsv (one row per PC5 region), regions_lost.tsv (PC3 regions with no PC5 match),
set_stress.tsv (PC3 vs PC5 set-level enrichment).
Run from the repo root: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/compare_regions.py
"""
import csv
import math
import os
import re

ROOT = os.getcwd()
B = os.path.join(ROOT, "ANALYSIS/Trifolium86MainK5ERpc5/regions")
PC3 = os.path.join(ROOT, "Trifolium86MainK5ER_results/GEA/tables")
PC5 = os.path.join(ROOT, "Trifolium86MainK5ERpc5_results/GEA/tables")
STRESS = {"PC3": os.path.join(ROOT, "ANALYSIS/Trifolium86MainK5/stress/A_ER_2of2"), "PC5": os.path.join(B, "stress")}
REPORT = os.path.join(ROOT, "ANALYSIS/Trifolium86MainK5/regionsER/REGIONS_ER.md")


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, table):
    with open(os.path.join(B, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]), delimiter="\t")
        w.writeheader()
        w.writerows(table)
    print(f"wrote {name}: {len(table)} rows")


def hyper_sf(k, K, N, n):
    """P(X >= k), X ~ Hypergeometric(N, K, n) (copied from compare_sets.py)."""
    if k <= 0:
        return 1.0
    lc = lambda a, b: math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    den = lc(N, n)
    return min(1.0, sum(math.exp(lc(K, i) + lc(N - K, n - i) - den) for i in range(k, min(K, n) + 1)))


ov = lambda a, b: a["chr"] == b["chr"] and int(a["start"]) <= int(b["end"]) and int(b["start"]) <= int(a["end"])
# report section number + title per PC3 region id
sect = {m.group(2): (m.group(1), m.group(3)) for m in
        re.finditer(r"^## (\d+)\. (\S+) — (.+)$", open(REPORT).read(), re.M)}
r3, r5 = rows(os.path.join(PC3, "regions_combined.tsv")), rows(os.path.join(PC5, "regions_combined.tsv"))
g3, g5 = ({} for _ in range(2))
for d, f in ((g3, PC3), (g5, PC5)):
    for g in rows(os.path.join(f, "genes_combined.tsv")):
        d.setdefault(g["region_id"], set()).add(g["gene_id"])
s3 = {r["SNPID"] for r in rows(os.path.join(PC3, "selected_snps.tsv"))}
s5 = {r["SNPID"] for r in rows(os.path.join(PC5, "selected_snps.tsv"))}
ks5 = {r["region_id"]: r for r in rows(os.path.join(STRESS["PC5"], "regions_stress_long.tsv")) if r["line"] == "ks_any"}
go5 = {r["region_id"]: r for r in rows(os.path.join(STRESS["PC5"], "regions_stress_long.tsv")) if r["line"] == "go_any"}

out = []
for r in r5:
    m = [o for o in r3 if ov(r, o)]
    snps = r["snp_ids"].split(",")
    genes = g5.get(r["region_id"], set())
    old_genes = set().union(*(g3.get(o["region_id"], set()) for o in m)) if m else set()
    out.append(dict(
        region_id=r["region_id"], status="old" if m else "new",
        pc3_region=",".join(o["region_id"] for o in m),
        report_section=",".join(f"{sect[o['region_id']][0]}. {sect[o['region_id']][1]}" for o in m if o["region_id"] in sect),
        investigated="yes" if m else "no",
        snp_count=r["snp_count"], snps_shared_with_pc3=sum(x in s3 for x in snps), snp_ids=r["snp_ids"],
        traits=r["traits"], min_pvalue=r["min_pvalue"],
        n_genes=len(genes), genes_in_pc3_region=len(genes & old_genes),
        strict_stress_arrays=ks5.get(r["region_id"], {}).get("hit_arrays", 0),
        strict_exp=ks5.get(r["region_id"], {}).get("exp_arrays", ""),
        go_stress_arrays=go5.get(r["region_id"], {}).get("hit_arrays", 0),
        go_exp=go5.get(r["region_id"], {}).get("exp_arrays", ""),
        stress_genes=ks5.get(r["region_id"], {}).get("hit_desc", "")[:300]))
write("regions_v2.tsv", out)

lost = [dict(region_id=o["region_id"], report_section=" ".join(sect.get(o["region_id"], ("", "")))
             , snp_count=o["snp_count"], snp_ids=o["snp_ids"], n_genes=len(g3.get(o["region_id"], set())),
             snps_still_selected=sum(x in s5 for x in o["snp_ids"].split(",")))
        for o in r3 if not any(ov(o, r) for r in r5)]
if lost:
    write("regions_lost.tsv", lost)

st = []
for tag, regs, d in (("PC3", r3, STRESS["PC3"]), ("PC5", r5, STRESS["PC5"])):
    long = rows(os.path.join(d, "regions_stress_long.tsv"))
    rates = {r["line"]: r for r in rows(os.path.join(d, "genome_stress_rates.tsv"))}
    row = dict(set=tag, regions=len(regs), snps=len(s3 if tag == "PC3" else s5),
               genes=len(set().union(*(g3 if tag == "PC3" else g5).values())))
    for line, lab in (("ks_any", "kw_strict"), ("go_any", "go")):
        L = [r for r in long if r["line"] == line]
        obs, exp = sum(int(r["hit_arrays"]) for r in L), sum(float(r["exp_arrays"]) for r in L)
        n_a, rt = sum(int(r["n_arrays"]) for r in L), rates[line]
        Na = round(int(rt["genome_arrays"]) / float(rt["array_rate"]))
        row.update({f"{lab}_obs": obs, f"{lab}_exp": f"{exp:.1f}", f"{lab}_ratio": f"{obs / exp:.2f}" if exp else "",
                    f"{lab}_p": f"{hyper_sf(obs, int(rt['genome_arrays']), Na, n_a):.2g}"})
    st.append(row)
write("set_stress.tsv", st)
n_old = sum(r["status"] == "old" for r in out)
print(f"PC5 regions {len(out)}: old {n_old}, new {len(out) - n_old}; PC3 regions lost {len(lost)} of {len(r3)}; "
      f"SNPs PC3 {len(s3)} / PC5 {len(s5)} / shared {len(s3 & s5)}")
