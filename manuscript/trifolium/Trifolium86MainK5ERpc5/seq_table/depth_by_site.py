#!/usr/bin/env python3
"""Does read depth, or the missing calls it causes, differ among the sites of Trifolium israeliticum?

Input: t_sequencing.tsv (build_t_sequencing.py). Over the plants the framework retained:
  - Kruskal-Wallis tests of mean depth and of the proportion of missing calls across the 11 sites,
    each with a permutation check (9,999 permutations of the site labels, seed 1);
  - Spearman's rho between the site medians of depth and the site order (1-11, south to north).
Over the retained and over all called plants: Spearman's rho between mean depth and missing calls.

Usage: python3 -I depth_by_site.py <t_sequencing.tsv>
Needs numpy and scipy (scipy 1.17.1 and numpy 2.4.2 in the run of 2026-10-25).
"""
import csv
import statistics
import sys

import numpy as np
from scipy import stats


def kruskal_with_permutation(values, sites, n_perm=9999, seed=1):
    groups = sorted(set(sites))

    def h(lab):
        return stats.kruskal(*[values[lab == g] for g in groups])[0]

    res = stats.kruskal(*[values[sites == g] for g in groups])
    rng = np.random.default_rng(seed)
    k = sum(h(rng.permutation(sites)) >= res.statistic for _ in range(n_perm))
    return res.statistic, len(groups) - 1, res.pvalue, (k + 1) / (n_perm + 1)


def main(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    kept = [r for r in rows if r["outcome"] == "retained"]
    sites = np.array([int(r["site"]) for r in kept])
    depth = np.array([float(r["mean_depth"]) for r in kept])
    miss = np.array([float(r["missing_calls_pct"]) for r in kept])
    print(f"retained plants: {len(kept)} at {len(set(sites))} sites")
    for label, v in (("mean depth", depth), ("missing calls", miss)):
        H, df, p, p_perm = kruskal_with_permutation(v, sites)
        print(f"Kruskal-Wallis, {label} across sites: H = {H:.2f}, df = {df}, P = {p:.3f} "
              f"(permutation P = {p_perm:.3f})")
    order = sorted(set(sites))
    medians = [statistics.median(depth[sites == s]) for s in order]
    rho = stats.spearmanr(order, medians)
    print("site median depth, south to north:", [round(float(m), 1) for m in medians])
    print(f"Spearman, site median depth vs site order: rho = {rho.statistic:.2f}, P = {rho.pvalue:.2f}")
    for label, sel in (("retained", kept), ("all called", rows)):
        d = [float(r["mean_depth"]) for r in sel]
        m = [float(r["missing_calls_pct"]) for r in sel]
        rho = stats.spearmanr(d, m)
        print(f"Spearman, mean depth vs missing calls, {label} (n = {len(sel)}): rho = {rho.statistic:.3f}, "
              f"P = {rho.pvalue:.1e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
