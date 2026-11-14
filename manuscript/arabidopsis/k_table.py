#!/usr/bin/env python3
"""k_table.py: the evidence for K, one row per K, and for EMMAX's n_pcs, one row per PC count.

Joins tables that earlier steps wrote; it computes nothing new beyond sums, maxima and the length of the
broken-stick leading run. Plan valiant-mixing-sparrow Phase 4: the user picks one K (sNMF.k_best = LFMM K =
BLINK n_pcs = the geometric offset's K) and GEA EMMAX n_pcs from these two tables.

Reads, for <project>:
  ANALYSIS/<project>/k_evidence/{curve,solutions,reproducibility,composition}.tsv   k_evidence.py
  ANALYSIS/<project>/pregea/rung_hit_overlap.tsv                                   pregea_hit_overlap.R
  ANALYSIS/<project>/pregea/pc_climate_r2.tsv                                      pc_climate_r2.R
  <project>_results/PreGEA/tables/{lfmm,emmax}/{lfmm,emmax}_ladder.tsv              mode=pregea, ladders
  <project>_results/PreGEA/tables/structure/pruned_pca_screeplot.tsv               mode=pregea, screeplot
  <project>_results/PreStructure/tables/tracy_widom.tsv                            mode=prestructure

Writes:
  ANALYSIS/<project>/k_evidence/k_table.tsv, one row per K of the LFMM ladder:
    K, cross-entropy (mean, sd, step from K-1, step over the pooled run-to-run sd, step ratio), sNMF solutions
    (count, runs in the best run's), clusters replicated, median worst r, share of samples over 60 %, smallest
    main cluster; LFMM (pooled lambda, largest per-predictor lambda, hits summed over predictors as production
    calls them, predictors with a hit, Jaccard with the reference K and with K + 1); and PC K - 1, the axis K
    clusters need beside the first K - 2 (variance share, above the broken stick, inside the leading run,
    Tracy-Widom p, R2 with the predictors over samples)
  ANALYSIS/<project>/pregea/emmax_table.tsv, one row per n_pcs of the EMMAX ladder:
    n_pcs, pooled lambda, smallest and largest per-predictor lambda, hits summed over predictors, Jaccard with
    n_pcs 0 and with one PC more, and for n_pcs >= 1 the last PC added (variance share, R2 with the predictors)

The broken-stick leading run is the number of axes from PC1 on that all exceed the broken stick; the
pipeline's own count includes tail axes above it (filed 2026-11-06), so the run is the number to read.

Usage: python3 -I manuscript/arabidopsis/k_table.py <project>   (from the repo root; stdlib only)
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    if not path.is_file():
        sys.exit(f"k_table.py: missing {path.relative_to(ROOT)}")
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path.relative_to(ROOT)} ({len(rows)} rows)")


def check(ok, what):
    if not ok:
        sys.exit(f"k_table.py: check failed: {what}")


def rnd(x, d=4):
    return "" if x is None or x == "" else f"{float(x):.{d}g}"


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: k_table.py <project>")
    project = sys.argv[1]
    an = ROOT / "ANALYSIS" / project
    res = ROOT / f"{project}_results"
    kev = an / "k_evidence"

    curve = {r["K"]: r for r in read(kev / "curve.tsv")}
    sol = read(kev / "solutions.tsv")
    rep = {r["K"]: r for r in read(kev / "reproducibility.tsv")}
    comp = {r["K"]: r for r in read(kev / "composition.tsv")}
    overlap = read(an / "pregea" / "rung_hit_overlap.tsv")
    pcr2 = {r["pc"]: r for r in read(an / "pregea" / "pc_climate_r2.tsv")}
    lfmm = read(res / "PreGEA/tables/lfmm/lfmm_ladder.tsv")
    emmax = read(res / "PreGEA/tables/emmax/emmax_ladder.tsv")
    scree = {r["pc"]: r for r in read(res / "PreGEA/tables/structure/pruned_pca_screeplot.tsv")}
    tw = {r["axis"]: r for r in read(res / "PreStructure/tables/tracy_widom.tsv")}

    leading = 0
    while str(leading + 1) in scree and scree[str(leading + 1)]["above_bs"] == "TRUE":
        leading += 1
    n_above = sum(r["above_bs"] == "TRUE" for r in scree.values())
    tw_sig = sum(float(r["tw_pvalue"]) < 0.05 for r in tw.values())

    def ladder_by_rung(rows):
        out = {}
        for r in rows:
            out.setdefault(r["rung_value"], {})[r["trait"]] = r
        return out

    lf = ladder_by_rung(lfmm)
    em = ladder_by_rung(emmax)
    ov = {(r["engine"], r["rung"]): r for r in overlap}
    lfmm_ks = sorted(lf, key=int)
    check(all(k in curve for k in lfmm_ks), "every LFMM rung has a cross-entropy row")
    check(all(("lfmm", f"K{k}") in ov for k in lfmm_ks), "every LFMM rung has an overlap row")

    rows = []
    for k in lfmm_ks:
        c, rp, cm = curve[k], rep[k], comp[k]
        ks = [s for s in sol if s["K"] == k]
        best = [s for s in ks if s["holds_best_run"] == "True"]
        check(len(best) == 1, f"K {k}: one solution holds the best run")
        preds = {t: r for t, r in lf[k].items() if t != "__pooled__"}
        hits = {t: int(r["hits_bonf"]) for t, r in preds.items()}
        o = ov[("lfmm", f"K{k}")]
        check(int(o["n_hits"]) == sum(hits.values()), f"K {k}: overlap n_hits = summed ladder hits")
        pc = str(int(k) - 1)
        s = scree[pc]
        rows.append({
            "K": k, "ce_mean": rnd(c["mean"], 6), "ce_sd": rnd(c["sd"], 3), "ce_step": rnd(c["step"], 3),
            "ce_step_over_sd": rnd(c["step_over_pooled_sd"], 3), "ce_step_ratio": rnd(c["step_ratio"], 3),
            "snmf_solutions": len(ks), "runs_in_best_solution": best[0]["n_runs"],
            "clusters_replicated": rp["replicated"], "median_worst_r": rnd(rp["median_worst_r"], 3),
            "share_over_60": rnd(cm["share_over_60"], 3), "smallest_main_cluster": cm["smallest_main"],
            "lfmm_lambda_pooled": rnd(lf[k]["__pooled__"]["lambda_gc"], 3),
            "lfmm_lambda_max": rnd(max(float(r["lambda_gc"]) for r in preds.values()), 3),
            "lfmm_lambda_max_predictor": max(preds, key=lambda t: float(preds[t]["lambda_gc"])),
            "lfmm_hits": sum(hits.values()), "lfmm_predictors_with_hits": sum(v > 0 for v in hits.values()),
            "lfmm_hits_by_predictor": ",".join(f"{t}:{hits[t]}" for t in sorted(hits, key=lambda t: int(t.split("_")[1]))),
            "jaccard_reference": o["reference"], "jaccard_with_reference": o["jaccard_with_reference"],
            "jaccard_with_next": o["jaccard_with_next"],
            "pc": pc, "pc_var_pct": rnd(100 * float(s["var_explained"]), 3), "pc_above_broken_stick": s["above_bs"],
            "pc_in_leading_run": str(int(pc) <= leading).upper(),
            "pc_tw_p": rnd(tw[pc]["tw_pvalue"], 3) if pc in tw else "",
            "pc_climate_r2": rnd(pcr2[pc]["r2_samples"], 3) if pc in pcr2 else "",
        })
    write(kev / "k_table.tsv", rows)

    erows = []
    for n in sorted(em, key=int):
        preds = {t: r for t, r in em[n].items() if t != "__pooled__"}
        lams = [float(r["lambda_gc"]) for r in preds.values()]
        hits = sum(int(r["hits_bonf"]) for r in preds.values())
        o = ov[("emmax", f"npc{n}")]
        check(int(o["n_hits"]) == hits, f"n_pcs {n}: overlap n_hits = summed ladder hits")
        erows.append({
            "n_pcs": n, "emmax_lambda_pooled": rnd(em[n]["__pooled__"]["lambda_gc"], 3),
            "emmax_lambda_min": rnd(min(lams), 3), "emmax_lambda_max": rnd(max(lams), 3), "emmax_hits": hits,
            "jaccard_with_npc0": o["jaccard_with_reference"], "jaccard_with_next": o["jaccard_with_next"],
            "last_pc_var_pct": rnd(100 * float(scree[n]["var_explained"]), 3) if n != "0" else "",
            "last_pc_climate_r2": rnd(pcr2[n]["r2_samples"], 3) if n != "0" else "",
        })
    write(an / "pregea" / "emmax_table.tsv", erows)
    print(f"PCA: {leading} leading axes above the broken stick ({n_above} above it in all); "
          f"Tracy-Widom p < 0.05 on {tw_sig} of {len(tw)} axes")


if __name__ == "__main__":
    main()
