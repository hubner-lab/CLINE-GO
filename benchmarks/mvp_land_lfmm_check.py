#!/usr/bin/env python3
"""Job A of the landscapes pack: is the LFMM2 offset result an error?

LFMM2 ranks the causal loci best on all three landscapes, yet its offset is the least accurate in
82 % of SS-Mtn and 69 % of Est-Clines replicates. This script checks that result by a code path of
its own (Python + scipy; none of the R scorers), reading what the runs wrote:

  sets      every replicate x LFMM, EMMAX, RDA: the set file holds exactly the SNPs the p-value table
            gives at the stated rule. The rule is compute_pval_threshold()'s toppct
            (scripts/R/utils/pval_threshold.R): N = ceiling(share * rows / 100 - 1e-9), NA rows
            counted in rows, cutoff = the N-th lowest p, selection p <= cutoff (so a tie at the
            cutoff can give more than N), the predictors' lists pooled. Share 0.25 per predictor,
            RDA 0.5 on its one multivariate p. Also: the K in the file name is the manifest's k_best.
  fitted    every replicate x six sets: the set size in the file, in snp_set.yaml, and the SNPs each
            offset method says it fitted (Gradient Forest and geometric logs, RDA diagnostics, whose
            snp_set_name must also name the folder). A set under 3 SNPs has no offset run.
  accuracy  20 replicates per landscape (4 per demography, fixed draw, written down) x the LFMM, EMMAX
            and RDA sets x the three offset methods: read the tree's own
            Maladaptation/tables/<method>/gea.<set>_nospatial/<garden>/genetic_offset_site.tsv, keep
            the 100 landscape gardens of garden_fitness_<seed>.tsv, take the mean offset per source
            site, Kendall tau-b against fitness per garden, the median over gardens, accuracy =
            -median tau, then the median over the three methods; compare with rank_long.tsv to 1e-9.
            SS-Clines is the control. Every tree file read is also compared byte for byte with its
            harvested copy under gardens/.
  labels    every replicate: rank_long.tsv's accuracy equals -tau of phase1_seed_medians_solo.tsv
            under this script's own mapping of set and offset-method names.

Writes <out>/lfmm_check.tsv (check, landscape, n, n_agree, max_abs_diff, note), the per-replicate
tables under <out>/lfmm_check_details/, and the drawn replicates. Exit 0 whatever it finds: the
summary is the result.
Usage: python3 -I benchmarks/mvp_land_lfmm_check.py [--procs 8] [--part sets,fitted,accuracy,labels]
"""
import argparse
import glob
import math
import multiprocessing as mp
import os
import random
import re

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = "/mnt/data/eugene/ADAPTOGENE"
EVAL = f"{ROOT}/benchmarks/mvp_eval"
OUT = f"{EVAL}/landscape_ms"
DET = f"{OUT}/lfmm_check_details"
LANDS = {
    "SS-Clines": dict(seeds=f"{ROOT}/work/ld_phase3_2026-10-26/scripts/seeds600.txt", off=f"{EVAL}/offset14",
                      rank=f"{EVAL}/figures_ld_offset14/rank/rank_long.tsv",
                      tree=f"{ROOT}/_archive/mvp_ssclines/results/MVP{{s}}hd_results",
                      logs=f"{ROOT}/_archive/mvp_ssclines/logs/MVP{{s}}hd_logs"),
    "SS-Mtn": dict(seeds=f"{ROOT}/work/mvp_land/seeds_ssmtn.txt", off=f"{EVAL}/offset15_ssmtn",
                   rank=f"{EVAL}/figures_land_ssmtn_offset15_ssmtn/rank/rank_long.tsv",
                   tree=f"{ROOT}/MVP{{s}}hd_results", logs=f"{ROOT}/MVP{{s}}hd_logs"),
    "Est-Clines": dict(seeds=f"{ROOT}/work/mvp_land/seeds_estclines.txt", off=f"{EVAL}/offset15_estclines",
                       rank=f"{EVAL}/figures_land_estclines_offset15_estclines/rank/rank_long.tsv",
                       tree=f"{ROOT}/MVP{{s}}hd_results", logs=f"{ROOT}/MVP{{s}}hd_logs"),
}
# set folder -> rank_long `panel`; offset folder -> rank_long `engine`; offset folder -> method_label
PANEL = {"solo_lfmm": "LFMM", "solo_emmax": "EMMAX", "solo_rda": "RDA", "union": "1/3 methods",
         "best": "2/3 methods", "intersect3": "3/3 methods"}
ENGINE = {"gradient_forest": "gradientForest", "geometric_offset": "LFMM2", "rda_offset": "RDA"}
METHOD_LABEL = {"GFoffset": "gradientForest", "LFMM2offset": "LFMM2", "RDA-uncorrected": "RDA"}
MARKER_SET = {"gea_best": "2/3 methods", "gea_union": "1/3 methods", "gea_strict": "3/3 methods",
              "gea_lfmm_only": "LFMM", "gea_rda_only": "RDA", "gea_emmax_only": "EMMAX"}
RULE = {"LFMM": ("solo_lfmm", ["bio_1", "bio_2"], 0.25), "EMMAX": ("solo_emmax", ["bio_1", "bio_2"], 0.25),
        "RDA": ("solo_rda", ["climate_multivariate"], 0.5)}


def seeds_of(land):
    with open(LANDS[land]["seeds"]) as fh:
        return [x.strip() for x in fh if x.strip()]


def toppct(p, share):
    """compute_pval_threshold(p, 'toppct', share) + the inclusive selection; p with NaN for NA."""
    n_total = len(p)
    n = math.ceil(share * n_total / 100 - 1e-9)
    ok = np.sort(p[~np.isnan(p)])
    if len(ok) < n:
        return None, n, 0
    cutoff = ok[n - 1]
    keep = p <= cutoff
    return keep, n, int(keep.sum()) - n


def check_sets(job):
    land, s, k_best = job
    tree = LANDS[land]["tree"].format(s=s)
    rows = []
    for meth, (setname, cols, share) in RULE.items():
        fs = [f for f in glob.glob(f"{tree}/GEA/tables/methods/{meth}/{meth}_pvalues_K*.tsv") if "_sig_" not in f]
        if len(fs) != 1:
            rows.append(dict(landscape=land, seed=s, method=meth, status=f"{len(fs)} p tables"))
            continue
        k_file = int(re.search(r"_K(\d+)\.tsv$", fs[0]).group(1))
        pv = pd.read_csv(fs[0], sep="\t", dtype={"SNPID": str, "chr": str})
        want = np.zeros(len(pv), dtype=bool)
        ns, extra = [], 0
        for c in cols:
            keep, n, ex = toppct(pv[c].to_numpy(dtype=float), share)
            ns.append(n); extra += ex
            want |= keep
        want_ids = set(pv.loc[want, "SNPID"])
        have = set(pd.read_csv(f"{tree}/_intermediate/snp_sets/gea.{setname}/selected_snps.tsv", sep="\t",
                               dtype={"SNPID": str}, usecols=["SNPID"])["SNPID"])
        rows.append(dict(landscape=land, seed=s, method=meth, status="ok", k_file=k_file, k_best=k_best,
                         n_tested=len(pv), n_rule=ns[0], n_want=len(want_ids), n_have=len(have),
                         n_missing=len(want_ids - have), n_extra=len(have - want_ids), tie_extra=extra))
    return rows


def first_int(path, pattern):
    """All integer groups of the first line of `path` matching `pattern`; None when no such line."""
    if not os.path.exists(path):
        return None
    rx = re.compile(pattern)
    with open(path, errors="replace") as fh:
        for line in fh:
            m = rx.search(line)
            if m:
                return [int(x) for x in m.groups()]
    return None


def diag_values(path):
    if not os.path.exists(path):
        return None
    out = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def check_fitted(job):
    land, s = job
    tree = LANDS[land]["tree"].format(s=s)
    logs = LANDS[land]["logs"].format(s=s) + "/maladaptation"
    rows = []
    for setname in PANEL:
        d = f"{tree}/_intermediate/snp_sets/gea.{setname}"
        f = f"{d}/selected_snps.tsv"
        n_file = None
        if os.path.exists(f):
            with open(f) as fh:
                n_file = sum(1 for _ in fh) - 1
        ny = first_int(f"{d}/snp_set.yaml", r"^n_snps:\s*(\d+)")
        n_yaml = ny[0] if ny else None
        gf = first_int(f"{logs}/gradient_forest_adaptive_gea.{setname}_nospatial.log", r"Response variables \(SNPs\):\s*(\d+)")
        geo = first_int(f"{logs}/geometric_offset_gea.{setname}_nospatial.log", r"Candidate loci:\s*(\d+) of (\d+)")
        dg = diag_values(f"{tree}/Maladaptation/tables/rda_offset/gea.{setname}_nospatial/rda_offset_diagnostics.tsv")
        rows.append(dict(landscape=land, seed=s, set=setname, n_file=n_file, n_yaml=n_yaml,
                         n_gf=gf[0] if gf else None, n_geometric=geo[0] if geo else None,
                         n_rda_matched=int(dg["n_candidates_matched"]) if dg else None,
                         n_rda_zero_var=int(dg["n_zero_variance_dropped"]) if dg else None,
                         rda_set_name=dg.get("snp_set_name") if dg else None))
    return rows


def check_accuracy(job):
    land, s = job
    tree = LANDS[land]["tree"].format(s=s)
    off = LANDS[land]["off"]
    # round_trip: pandas' default parser is not exact and merges offsets one floating-point step apart, which R keeps
    # apart; Kendall's tau-b sees that as a tie (measured: tau -0.85444 vs -0.85327 in one garden).
    fit = pd.read_csv(f"{off}/garden_fitness_{s}.tsv", sep="\t", dtype={"seed": str, "source_site": str},
                      float_precision="round_trip")
    fit = fit[fit["type"] == "landscape"]
    gardens = sorted(fit["garden_id"].unique())
    rows, n_files, n_diff = [], 0, 0
    for setname in ("solo_lfmm", "solo_emmax", "solo_rda"):
        accs, accs_tied = {}, {}
        for meth in ENGINE:
            taus, taus_tied = [], []
            for g in gardens:
                f = f"{tree}/Maladaptation/tables/{meth}/gea.{setname}_nospatial/{g}/genetic_offset_site.tsv"
                h = f"{off}/gardens/{s}/{g}/{meth}__{setname}.tsv"
                with open(f, "rb") as a, open(h, "rb") as b:
                    n_files += 1
                    if a.read() != b.read():
                        n_diff += 1
                go = pd.read_csv(f, sep="\t", dtype={"site": str, "sample": str}, float_precision="round_trip")
                os_ = go.groupby("site")["genetic_offset"].mean()
                ft = fit[fit["garden_id"] == g].set_index("source_site")["fitness"]
                common = os_.index.intersection(ft.index)
                x, y = os_.loc[common].to_numpy(), ft.loc[common].to_numpy()
                ok = np.isfinite(x) & np.isfinite(y)
                t = kendalltau(x[ok], y[ok], variant="b").statistic
                if np.isfinite(t):
                    taus.append(t)
                # the same tau when offsets within 1e-12 count as tied: how much the floating-point steps move it
                t2 = kendalltau(np.round(x[ok], 12), y[ok], variant="b").statistic
                if np.isfinite(t2):
                    taus_tied.append(t2)
            accs[meth] = -float(np.median(taus))
            accs_tied[meth] = -float(np.median(taus_tied))
            rows.append(dict(landscape=land, seed=s, set=setname, engine=ENGINE[meth], accuracy=accs[meth],
                             accuracy_tied_1e12=accs_tied[meth], n_gardens=len(taus)))
        rows.append(dict(landscape=land, seed=s, set=setname, engine="median of 3",
                         accuracy=float(np.median(list(accs.values()))),
                         accuracy_tied_1e12=float(np.median(list(accs_tied.values()))), n_gardens=len(gardens)))
    return rows, n_files, n_diff


def summary_row(check, land, n, n_agree, max_abs_diff, note):
    return dict(check=check, landscape=land, n=n, n_agree=n_agree, max_abs_diff=max_abs_diff, note=note)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--part", default="sets,fitted,accuracy,labels")
    a = ap.parse_args()
    parts = a.part.split(",")
    os.makedirs(DET, exist_ok=True)
    man = pd.read_csv(f"{ROOT}/benchmarks/mvp_seeds.tsv", sep="\t", dtype={"seed": str}).set_index("seed")
    summ = []
    with mp.Pool(a.procs) as pool:
        if "sets" in parts:
            jobs = [(L, s, int(man.loc[s, "k_best"])) for L in LANDS for s in seeds_of(L)]
            R = pd.DataFrame([r for rows in pool.imap_unordered(check_sets, jobs, 8) for r in rows])
            R.to_csv(f"{DET}/sets.tsv", sep="\t", index=False)
            for (L, m), D in R.groupby(["landscape", "method"], sort=False):
                good = (D.status == "ok") & (D.n_missing == 0) & (D.n_extra == 0)
                note = (f"{int((D.tie_extra > 0).sum())} replicates with a tie at the cutoff "
                        f"(the inclusive rule adds {int(D.tie_extra.sum())} SNPs)")
                summ.append(summary_row(f"set = p-table rule ({m})", L, len(D), int(good.sum()),
                                        float((D.n_missing + D.n_extra).max()), note))
                good_k = D.k_file == D.k_best
                summ.append(summary_row(f"K in the file = k_best ({m})", L, len(D), int(good_k.sum()),
                                        float((D.k_file - D.k_best).abs().max()), ""))
            pd.DataFrame(summ).to_csv(f"{OUT}/lfmm_check.tsv", sep="\t", index=False)
            print("sets done", flush=True)
        if "fitted" in parts:
            jobs = [(L, s) for L in LANDS for s in seeds_of(L)]
            R = pd.DataFrame([r for rows in pool.imap_unordered(check_fitted, jobs, 8) for r in rows])
            R.to_csv(f"{DET}/fitted.tsv", sep="\t", index=False)
            for (L, st), D in R.groupby(["landscape", "set"], sort=False):
                has = D.n_file.fillna(0) >= 3
                U = D[has]
                note = f"{int((~has).sum())} replicates under 3 SNPs: no offset run expected"
                summ.append(summary_row(f"set file = snp_set.yaml ({st})", L, len(D),
                                        int((D.n_file.fillna(0) == D.n_yaml.fillna(0)).sum()),
                                        float((D.n_file.fillna(0) - D.n_yaml.fillna(0)).abs().max()), ""))
                for lab, col in (("Gradient Forest", U.n_gf), ("Geometric offset", U.n_geometric),
                                 ("RDA-based offset", U.n_rda_matched)):
                    miss = int(col.isna().sum())
                    diff = (col - U.n_file).abs()
                    summ.append(summary_row(f"SNPs fitted = set size, {lab} ({st})", L, len(U), int((diff == 0).sum()),
                                            float(diff.max()) if len(U) else 0.0,
                                            f"{note}; {miss} without a record"))
                zv = U[U.n_rda_zero_var.fillna(0) > 0]
                name_ok = (U.rda_set_name == f"gea.{st}")
                summ.append(summary_row(f"RDA diagnostics name the set ({st})", L, len(U), int(name_ok.sum()), 0.0,
                                        f"{len(zv)} replicates with zero-variance SNPs dropped from the RDA fit"))
            N = R[R.n_file.fillna(0) < 3]
            ran = N[N[["n_gf", "n_geometric", "n_rda_matched"]].notna().any(axis=1)]
            summ.append(summary_row("no offset run for a set under 3 SNPs", "all", len(N), len(N) - len(ran), float(len(ran)), ""))
            pd.DataFrame(summ).to_csv(f"{OUT}/lfmm_check.tsv", sep="\t", index=False)
            print("fitted done", flush=True)
        if "accuracy" in parts:
            rng = random.Random(20261106)
            drawn = []
            for L in LANDS:
                s = sorted(seeds_of(L))
                dem = man.loc[s, "demog_level_sub"]
                for d in sorted(dem.unique()):
                    pool_d = sorted(dem.index[dem == d])
                    drawn += [(L, x) for x in rng.sample(pool_d, 4)]
            pd.DataFrame(drawn, columns=["landscape", "seed"]).to_csv(f"{DET}/accuracy_seeds.tsv", sep="\t", index=False)
            RR, nf, nd = [], 0, 0
            for rows, a_nf, a_nd in pool.imap_unordered(check_accuracy, drawn, 1):
                RR += rows; nf += a_nf; nd += a_nd
            R = pd.DataFrame(RR)
            R["panel"] = R["set"].map(PANEL)
            rl = pd.concat([pd.read_csv(LANDS[L]["rank"], sep="\t", dtype={"seed": str}).assign(landscape=L) for L in LANDS])
            M = R.merge(rl[["landscape", "seed", "engine", "panel", "accuracy"]], on=["landscape", "seed", "engine", "panel"],
                        suffixes=("", "_rank_long"), how="left")
            M["abs_diff"] = (M.accuracy - M.accuracy_rank_long).abs()
            M.to_csv(f"{DET}/accuracy.tsv", sep="\t", index=False)
            for (L, st, e), D in M.groupby(["landscape", "set", "engine"], sort=False):
                good = D.abs_diff.notna() & (D.abs_diff < 1e-9)
                summ.append(summary_row(f"accuracy = rank_long ({PANEL[st]}, {e})", L, len(D), int(good.sum()),
                                        float(D.abs_diff.max()), "own Kendall tau-b from the tree's site offsets; 100 landscape gardens"
                                        + ("; control" if L == "SS-Clines" else "")))
            for (L, e), D in M.groupby(["landscape", "engine"], sort=False):
                d = (D.accuracy - D.accuracy_tied_1e12).abs()
                summ.append(summary_row(f"floating-point near ties ({e})", L, len(D), int((d < 1e-9).sum()), float(d.max()),
                                        "information, not a check: accuracy when site offsets within 1e-12 count as tied"))
            # the swap the check is for: would another set's accuracy have matched as well?
            for L in LANDS:
                wrong = 0
                W = rl[(rl.landscape == L) & (rl.engine == "median of 3")]
                D = M[(M.landscape == L) & (M.engine == "median of 3")]
                for _, r in D.iterrows():
                    others = W[(W.seed == r.seed) & (W.panel != r.panel)]
                    wrong += int(((others.accuracy - r.accuracy).abs() < 1e-9).any())
                summ.append(summary_row("accuracy matches another set's too (median of 3)", L, len(D), len(D) - wrong, float(wrong),
                                        "counts ties in rank_long; a swap would show as a mismatch above"))
            summ.append(summary_row("tree offsets = harvested copies", "all", nf, nf - nd, float(nd), "byte comparison of every site table read"))
            pd.DataFrame(summ).to_csv(f"{OUT}/lfmm_check.tsv", sep="\t", index=False)
            print("accuracy done", flush=True)
    if "labels" in parts:
        for L in LANDS:
            med = pd.read_csv(f"{LANDS[L]['off']}/phase1_seed_medians_solo.tsv", sep="\t", dtype={"seed": str})
            seeds = set(seeds_of(L))
            med = med[med.seed.isin(seeds) & med.marker_set.isin(MARKER_SET) & med.method_label.isin(METHOD_LABEL)].copy()
            med["panel"] = med.marker_set.map(MARKER_SET); med["engine"] = med.method_label.map(METHOD_LABEL)
            rl = pd.read_csv(LANDS[L]["rank"], sep="\t", dtype={"seed": str})
            rl = rl[rl.engine != "median of 3"]
            M = rl.merge(med[["seed", "panel", "engine", "tau"]], on=["seed", "panel", "engine"], how="left")
            E = M[M["empty"]]
            U = M[~M["empty"]]
            d = (U.accuracy + U.tau).abs()
            summ.append(summary_row("rank_long accuracy = -tau of the scorer's table", L, len(U), int((d < 1e-9).sum()), float(d.max()),
                                    f"{len(E)} rows of sets without an offset; {int(E.tau.notna().sum())} of them have a tau"))
    pd.DataFrame(summ).to_csv(f"{OUT}/lfmm_check.tsv", sep="\t", index=False)
    S = pd.DataFrame(summ)
    bad = S[S.n != S.n_agree]
    print(S.to_string(index=False, max_colwidth=70), flush=True)
    print(f"\n{len(S)} checks, {len(bad)} with n != n_agree", flush=True)


if __name__ == "__main__":
    main()
