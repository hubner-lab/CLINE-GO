"""gea_check.py: post-run checks for mode=gea on one Arabidopsis project (plan shiny-seeking-scott Phase 10).

usage: python3 -I manuscript/arabidopsis/gea_check.py <project> <K> <expected clumping distance bp>

Reads only; prints counts, lambdas, set sizes and every consistency check, and exits 1 if a check fails.
"""
import glob
import os
import re
import sys

import pandas as pd
import yaml

project, K, dist = sys.argv[1], sys.argv[2], int(sys.argv[3])
ROOT = "/mnt/data/eugene/ADAPTOGENE"
R = f"{ROOT}/{project}_results"
L = f"{ROOT}/{project}_logs/GEA"
fails = []


def check(ok, what):
    print(("  OK   " if ok else "  FAIL ") + what)
    if not ok:
        fails.append(what)


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


methods = ["EMMAX", "LFMM", "RDA"]
sig = {}
print(f"== {project}, K {K}, expected distance {dist} bp")
print("\n== thresholds and calls")
for m in methods:
    th = read(f"{R}/GEA/tables/methods/{m}/{m}_thresholds_K{K}.tsv")
    s = read(f"{R}/GEA/tables/methods/{m}/{m}_pvalues_K{K}_sig_snps_bonf_0.05.tsv")
    sig[m] = set(s["SNPID"])
    n_tested = sorted(set(th["n_tested"]))
    thr = sorted(set(th["threshold"]))
    per = ", ".join(f"{t} {n}" for t, n in zip(th["trait"], th["n_significant"]))
    print(f"{m}: n_tested {n_tested}, threshold {thr}, distinct sig SNPs {len(sig[m])}")
    print(f"  per predictor: {per}")
    none = [t for t, n in zip(th["trait"], th["n_significant"]) if n == "0"]
    print(f"  predictors with no call: {', '.join(none) if none else 'none'}")
    od = sorted(set(s["overlap_distance"])) if len(s) else []
    check(od in ([str(dist)], []), f"{m} sig table overlap_distance {od} == {dist}")

print("\n== lambda")
em = read(f"{R}/GEA/tables/methods/EMMAX/EMMAX_lambda_K{K}.tsv")
lam = em["lambda_gc"].astype(float)
print(f"EMMAX lambda_gc min {lam.min():.3f} ({em.loc[lam.idxmin(), 'trait']}), "
      f"max {lam.max():.3f} ({em.loc[lam.idxmax(), 'trait']}), median {lam.median():.3f}; "
      f"correction {sorted(set(em['correction']))}")
print("  " + ", ".join(f"{t} {float(v):.3f}" for t, v in zip(em["trait"], em["lambda_gc"])))
gifs = {}
for f in sorted(glob.glob(f"{L}/assoc_lfmm_bio_*.log")):
    t = re.search(r"assoc_lfmm_(bio_\d+)\.log", f).group(1)
    hits = re.findall(r"\$ gif\s*: num ([0-9.eE+-]+)", open(f).read())
    gifs[t] = float(hits[-1]) if hits else float("nan")
g = pd.Series(gifs)
print(f"LFMM gif (job logs, the factor applied) min {g.min():.2f} ({g.idxmin()}), max {g.max():.2f} "
      f"({g.idxmax()}), median {g.median():.2f}, n {g.notna().sum()}")
print("  " + ", ".join(f"{t} {v:.2f}" for t, v in sorted(gifs.items(), key=lambda kv: int(kv[0][4:]))))
lf = read(f"{R}/GEA/tables/methods/LFMM/LFMM_lambda_K{K}.tsv")
print(f"  LFMM_lambda table values {sorted(set(lf['lambda_gc']))} (1 by construction)")
d = dict(read(f"{R}/GEA/tables/methods/RDA/RDA_diagnostics_K{K}.tsv").values.tolist())
keys = ["fit_mode", "n_markers_fitted", "rda_axes", "adj_r_squared", "max_vif", "gif_lambda", "cpu",
        "marker_envelope_status", "n_candidates_bonf_0.01", "n_candidates_bonf_0.05",
        "n_candidates_q_0.05", "n_candidates_q_0.10", "candidate_threshold", "n_candidates_total"]
print("RDA: " + "; ".join(f"{k} {d.get(k, 'NA')}" for k in keys))
an = read(f"{R}/GEA/tables/methods/RDA/RDA_anova_K{K}.tsv")
print("  anova: " + "; ".join(f"{r.term} F {r.F} p {r.p}" for r in an.itertuples() if r.term == "Model"))

print("\n== imputed MAF filter")
rep = glob.glob(f"{R}/_work/*/*K{K}imp_maf*_report.tsv")
for f in rep:
    print(os.path.basename(f))
    print(read(f).to_string(index=False))

print("\n== six sets")
sel = read(f"{R}/GEA/tables/selected_snps.tsv")
sel_ids = set(sel["SNPID"])
sets = {}
for name in ["emmax", "lfmm", "rda", "union", "two_of_three", "three_of_three"]:
    d_ = f"{R}/_intermediate/snp_sets/gea.{name}"
    y = yaml.safe_load(open(f"{d_}/snp_set.yaml"))
    ids = set(read(f"{d_}/selected_snps.tsv")["SNPID"])
    sets[name] = ids
    print(f"gea.{name}: {len(ids)} SNPs, strategy {y.get('strategy')}, clumping_distance "
          f"{y.get('clumping_distance')}")
    check(int(y.get("clumping_distance")) == dist, f"gea.{name} clumping_distance == {dist}")
    if len(ids) < 3:
        print(f"  NOTE gea.{name} has under 3 SNPs: Maladaptation.snp_sets all skips it (attempts row)")
check(sets["three_of_three"] <= sets["two_of_three"], "three_of_three inside two_of_three")
check(sets["two_of_three"] <= sets["union"], "two_of_three inside union")
check(sets["union"] == sel_ids, "union == GEA/tables/selected_snps.tsv")
for m in methods:
    check(sets[m.lower()] == sig[m], f"gea.{m.lower()} == {m} sig table")
check(sel_ids == sig["EMMAX"] | sig["LFMM"] | sig["RDA"], "selected_snps == union of the three sig tables")

print("\n== regions and genes")
for f in [f"{L}/create_regions.log", f"{L}/find_sig_snps_EMMAX_bonf_0.05.log"]:
    lines = [ln.strip() for ln in open(f) if "distance" in ln.lower()]
    print(f"{os.path.basename(f)}: {lines[:2]}")
    check(any(str(dist) in ln for ln in lines), f"{os.path.basename(f)} names {dist} bp")
rpt = read(f"{R}/GEA/tables/regions_per_trait.tsv")
rc = read(f"{R}/GEA/tables/regions_combined.tsv")
w = rc["length"].astype(int)
print(f"regions: per trait {len(rpt)}, combined {len(rc)}; combined span {w.sum() / 1e6:.1f} Mb, "
      f"median width {w.median() / 1e3:.0f} kb, largest {w.max() / 1e3:.0f} kb")
gc = read(f"{R}/GEA/tables/genes_combined.tsv")
print(f"genes in combined regions: {gc['gene_id'].nunique()} distinct ({len(gc)} rows)")

print("\n== WZA significant windows")
for m in methods:
    f = f"{R}/GEA/tables/methods/{m}/{m}_wza_K{K}_sig_windows_bonf_0.05.tsv"
    s = read(f)
    print(f"{m}: {len(s)} rows, {s['SNPID'].nunique() if len(s) else 0} distinct windows")
    if len(s):
        wid = s["SNPID"].str.extract(r":(\d+)-(\d+)$").astype(int)
        span = sorted(set((wid[1] - wid[0] + 1).tolist()))[:3]
        print(f"  window widths (first 3 distinct): {span}")

print(f"\n{len(fails)} check(s) failed" + (": " + "; ".join(fails) if fails else ""))
sys.exit(1 if fails else 0)
