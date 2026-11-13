"""gwas_check.py: post-run checks for mode=gwas and mode=gea_x_gwas on one Arabidopsis project (Phase 10).

usage: python3 -I manuscript/arabidopsis/gwas_check.py <project> <K> <expected clumping distance bp>
"""
import os
import sys

import pandas as pd

project, K, dist = sys.argv[1], sys.argv[2], int(sys.argv[3])
ROOT = "/mnt/data/eugene/ADAPTOGENE"
R = f"{ROOT}/{project}_results"
L = f"{ROOT}/{project}_logs/GWAS"
fails = []
# Flowering genes the first project's notes name (TAIR ids).
KNOWN = {"FT": "AT1G65480", "DOG1": "AT5G45830", "FLC": "AT5G10140", "FRI": "AT4G00650",
         "SVP": "AT2G22540"}


def check(ok, what):
    print(("  OK   " if ok else "  FAIL ") + what)
    if not ok:
        fails.append(what)


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


print(f"== {project} GWAS, K {K}")
print(read(f"{R}/GWAS/tables/phenotype_missing_summary.tsv").iloc[:, :4].to_string(index=False))
sig = {}
for m in ["EMMAX", "BLINK"]:
    th = read(f"{R}/GWAS/tables/methods/{m}/{m}_thresholds_K{K}.tsv")
    la = read(f"{R}/GWAS/tables/methods/{m}/{m}_lambda_K{K}.tsv")
    s = read(f"{R}/GWAS/tables/methods/{m}/{m}_pvalues_K{K}_sig_snps_bonf_0.05.tsv")
    sig[m] = s
    print(f"{m}: " + "; ".join(f"{r.trait} n_tested {r.n_tested} thr {float(r.threshold):.3g} "
                               f"calls {r.n_significant}" for r in th.itertuples()))
    for r in la.itertuples():
        lam = float(r.lambda_gc)
        print(f"  lambda {r.trait} {lam:.3f} ({r.correction})"
              + ("" if 0.90 <= lam <= 1.15 else "  <- outside the plan's 0.90-1.15 band"))
    od = sorted(set(s["overlap_distance"])) if len(s) else []
    check(od in ([str(dist)], []), f"{m} sig table overlap_distance {od} == {dist}")
    if len(s):
        top = s.assign(p=s["pvalue"].astype(float)).sort_values("p").head(8)
        print("  top: " + ", ".join(f"{r.trait} {r.SNPID} {r.p:.2g}" for r in top.itertuples()))
ids = set(sig["EMMAX"]["SNPID"]) | set(sig["BLINK"]["SNPID"])
sel = read(f"{R}/GWAS/tables/selected_snps.tsv")
check(set(sel["SNPID"]) == ids, f"selected_snps ({len(sel)}) == union of the two sig tables ({len(ids)})")
rpt = read(f"{R}/GWAS/tables/regions_per_trait.tsv")
rc = read(f"{R}/GWAS/tables/regions_combined.tsv")
gc = read(f"{R}/GWAS/tables/genes_combined.tsv")
print(f"distinct sig SNPs {len(ids)}; regions per trait {len(rpt)}, combined {len(rc)}; "
      f"genes {gc['gene_id'].nunique()}")
for f in ["create_regions.log"]:
    p = f"{L}/{f}"
    if os.path.exists(p):
        lines = [ln.strip() for ln in open(p) if "distance" in ln.lower()]
        print(f"{f}: {lines[:2]}")
        check(any(str(dist) in ln for ln in lines), f"GWAS {f} names {dist} bp")
for name, tair in KNOWN.items():
    hit = gc[gc["gene_id"].str.contains(tair)]
    print(f"  {name} ({tair}): " + (f"in combined region(s) {sorted(set(hit['region_id']))}"
                                    if len(hit) else "no combined region"))
for m in ["EMMAX", "BLINK"]:
    w = read(f"{R}/GWAS/tables/methods/{m}/{m}_wza_K{K}_sig_windows_bonf_0.05.tsv")
    print(f"WZA {m}: {len(w)} significant window rows")

print(f"\n== {project} GEA x GWAS")
pw_path = f"{R}/GEAxGWAS/tables/pairwise_overlap_table.tsv"
if not os.path.exists(pw_path):
    print("pairwise table absent: mode=gea_x_gwas has not run on this project")
else:
    pw = read(pw_path)
    x = pw[pw["comparison_type"].str.contains("GWAS") & pw["comparison_type"].str.contains("GEA")]
    print(f"pairwise rows {len(pw)}; GEA-GWAS pairs {len(x)}: sharing a window "
          f"{(x['window_overlaps'].astype(int) > 0).sum()}, sharing a SNP {(x['exact_matches'].astype(int) > 0).sum()}")
    for r in x[(x["window_overlaps"].astype(int) > 0) | (x["exact_matches"].astype(int) > 0)].itertuples():
        print(f"  {r.trait_a} ({r.source_a}) x {r.trait_b} ({r.source_b}): exact {r.exact_matches}, "
              f"windows {r.window_overlaps}")
gea_rc = read(f"{R}/GEA/tables/regions_combined.tsv")
ov = []
for g in rc.itertuples():
    same = gea_rc[gea_rc["chr"] == g.chr]
    o = same[(same["start"].astype(int) <= int(g.end)) & (same["end"].astype(int) >= int(g.start))]
    if len(o):
        ov.append((g.region_id, ",".join(o["region_id"]), g.methods, ",".join(o["methods"])))
print(f"GWAS combined regions overlapping a GEA combined region: {len(ov)} of {len(rc)}")
for o in ov:
    print(f"  GWAS {o[0]} ({o[2]}) x GEA {o[1]} ({o[3]})")
print(f"\n{len(fails)} check(s) failed" + (": " + "; ".join(fails) if fails else ""))
sys.exit(1 if fails else 0)
