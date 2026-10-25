#!/usr/bin/env python3
"""Spearman vs Kendall for every rank statistic in the offset pack (user 2026-10-15: "Kendall, not Spearman?").

Recomputes, for the 27-SNP set (bio_1,3,5, mean of 3 GCMs), with Spearman rho, Kendall tau-b and weighted Kendall tau
(scipy weightedtau, hyperbolic weights on the highest-offset sites first):
  - random-consensus test: adaptive vs the mean-rank consensus of 199 random 27-SNP sets, leave-one-out null;
  - climate-null test: adaptive vs the climate-only site distance, against the 199 random sets;
  - method agreement (GF / geometric / RDA);
  - each predictor's projected change at the site vs the GF and linear-family ranks.
Also Kendall's W of the 199 random sets (agreement of many rankings at once).
Reads site_maps/ outputs and paper_pack/consensus_sites.tsv; writes site_maps/kendall_check.tsv. Run from the repo root:
  python3 -I ANALYSIS/Trifolium86MainK5/maladaptationER/kendall_check.py
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

B = os.path.join(os.getcwd(), os.environ.get("AN", "ANALYSIS/Trifolium86MainK5/maladaptationER"))  # env AN: v2 = PC5 rebuild
SET = "gea.ER_overlap"
STATS = {
    "spearman": lambda a, b: stats.spearmanr(a, b)[0],
    "kendall": lambda a, b: stats.kendalltau(a, b)[0],
    # high offset first: weightedtau weighs the top of the ranking by the larger value -> pass offsets as they are
    "wtau": lambda a, b: stats.weightedtau(a, b)[0],
}

so = pd.read_csv(os.path.join(B, "site_maps/site_offsets.tsv"), sep="\t", dtype={"site": str},
                 keep_default_na=False)  # method "null" would read as NaN
adapt = {m: so[(so.method == m) & (so.set == SET)].set_index("site").offset for m in ("gf", "geometric", "rda")}
null = so[so.method == "null"].set_index("site").offset
sites = sorted(null.index)

gfd = pd.read_csv(os.path.join(B, "site_maps/gf_random_null/site_offsets_draws.tsv"), sep="\t", dtype={"site": str})
gfd = gfd[gfd.kind == "random"].rename(columns={"seed": "draw"})
lrd = pd.read_csv(os.path.join(B, "site_maps/random_null/site_offsets_draws.tsv"), sep="\t", dtype={"site": str})
draws = {"gf": gfd.pivot(index="draw", columns="site", values="offset")[sites]}
for m in ("geometric", "rda"):
    draws[m] = lrd[lrd.engine == m].pivot(index="draw", columns="site", values="offset")[sites]

rows = []
for m, D in draws.items():
    a = adapt[m][sites].values
    R = D.rank(axis=1).values                      # rank 1 = lowest offset, per draw
    n, k = R.shape
    tot = R.sum(0)
    # Kendall's W of the random sets (no ties)
    S = ((tot - n * (k + 1) / 2) ** 2).sum()
    W = 12 * S / (n ** 2 * (k ** 3 - k))
    for s, f in STATS.items():
        loo = np.array([f(R[i], (tot - R[i]) / (n - 1)) for i in range(n)])
        obs = f(a, tot / n)
        nl = np.array([f(D.values[i], null[sites].values) for i in range(n)])
        obs_n = f(a, null[sites].values)
        rows += [
            dict(test="random_consensus", method=m, stat=s, adaptive=obs, random_median=np.median(loo),
                 q025=np.quantile(loo, .025), q975=np.quantile(loo, .975), p=(1 + (loo <= obs).sum()) / (1 + n)),
            dict(test="climate_null", method=m, stat=s, adaptive=obs_n, random_median=np.median(nl),
                 q025=np.quantile(nl, .025), q975=np.quantile(nl, .975), p=(1 + (nl <= obs_n).sum()) / (1 + n)),
        ]
    rows.append(dict(test="random_kendall_W", method=m, stat="W", adaptive=np.nan, random_median=W))

for x, y in (("gf", "geometric"), ("gf", "rda"), ("geometric", "rda")):
    for s, f in STATS.items():
        rows.append(dict(test="method_agreement", method=f"{x}~{y}", stat=s, adaptive=f(adapt[x][sites], adapt[y][sites])))

cs = pd.read_csv(os.path.join(B, "paper_pack/consensus_sites.tsv"), sep="\t", dtype={"site": str})
for p in ("bio_1", "bio_3", "bio_5"):
    for fam in ("gf_rank", "linear_rank"):
        # rank 1 = highest risk; correlate with -rank so a positive value = larger change, higher risk
        for s in ("spearman", "kendall"):
            rows.append(dict(test="delta_vs_rank", method=f"d_{p}~{fam}", stat=s,
                             adaptive=STATS[s](cs[f"d_{p}"].abs(), -cs[fam])))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(B, "site_maps/kendall_check.tsv"), sep="\t", index=False, float_format="%.3f")
pd.set_option("display.width", 160)
print(out.round(3).to_string(index=False))
