#!/usr/bin/env python3
"""Job B of the landscapes pack: the SS-Mtn and Est-Clines rows of the new tables against each landscape's own scored tables.

The gate in mvp_land_tables.R covers SS-Clines only, because only SS-Clines has published tables. This covers the two new
landscapes, by a second code path, against what the per-landscape runs wrote before this pack existed:
  t_ranks and t_ranks_strata   rank_summary.tsv (mvp_panel_rank.R): top 1-3 and bottom 1-3, all offset methods, pooled + strata
  t_detection_arch, t_detection  detection_stats.tsv: median true and false positives
  t_offset_accuracy            C_delta_stats.tsv: median difference from the causal loci
  fig_s1_numbers               B4_aucpr_stats.tsv: median AUC-PR of RDA, LFMM2, EMMAX
Writes tables/gate_new_landscapes.tsv (check, landscape, n, n_agree, max_abs_diff).
"""
import pandas as pd

ROOT = "/mnt/data/eugene/ADAPTOGENE/benchmarks/mvp_eval"
T = f"{ROOT}/landscape_ms/tables"
L = {"SS-Clines": ("figures_ld_offset14/rank", "figures_ld_gallery14"),
     "SS-Mtn": ("figures_land_ssmtn_offset15_ssmtn/rank", "figures_land_ssmtn_gallery"),
     "Est-Clines": ("figures_land_estclines_offset15_estclines/rank", "figures_land_estclines_gallery")}
ENG = {"gradientForest": "Gradient Forest", "LFMM2": "Geometric offset", "RDA": "RDA-based offset", "median of 3": "Median of the three"}
SET = {"2/3 methods": "2/3 set", "1/3 methods": "1/3 set", "3/3 methods": "3/3 set", "LFMM": "LFMM2", "RDA": "RDA", "EMMAX": "EMMAX"}
STAT = {"top1": "rank 1", "top2": "ranks 1–2", "top3": "ranks 1–3", "bottom1": "rank 6", "bottom2": "ranks 5–6", "bottom3": "ranks 4–6"}
BLOCK = {"ssclines_nvar_mvar": "N variable, m variable", "ssclines_ncline_ns": "N cline north-south", "ssclines_nequal_mconst": "N equal, m constant",
         "ssclines_ncline_ctredge": "N cline center-edge", "ssclines_nequal_mbreaks": "N equal, m breaks"}
ARCH = ["oligogenic", "moderately polygenic", "highly polygenic"]
TR = pd.read_csv(f"{T}/t_ranks.tsv", sep="\t")
TS = pd.read_csv(f"{T}/t_ranks_strata.tsv", sep="\t")
TD = pd.read_csv(f"{T}/t_detection.tsv", sep="\t")
TA = pd.read_csv(f"{T}/t_detection_arch.tsv", sep="\t")
TO = pd.read_csv(f"{T}/t_offset_accuracy.tsv", sep="\t")
S1 = pd.read_csv(f"{T}/fig_s1_numbers.tsv", sep="\t")
rows = []


def add(check, land, a, b):
    """a, b: Series aligned on the same index; a missing pair counts as a disagreement."""
    x = pd.concat([a.rename("a"), b.rename("b")], axis=1)
    d = (x.a - x.b).abs()
    ok = d < 1e-9
    rows.append(dict(check=check, landscape=land, n=len(x), n_agree=int(ok.sum()), max_abs_diff=float(d.max())))


for land, (rk, gal) in L.items():
    RS = pd.read_csv(f"{ROOT}/{rk}/rank_summary.tsv", sep="\t")
    RS = RS[RS.metric.isin(STAT)].assign(offset_method=lambda d: d.engine.map(ENG), set=lambda d: d.panel.map(SET), statistic=lambda d: d.metric.map(STAT))
    P = RS[RS.stratum == "pooled"].set_index(["offset_method", "set", "statistic"]).value
    N = TR[(TR.landscape == land) & TR.statistic.isin(STAT.values())].set_index(["offset_method", "set", "statistic"]).value
    add("t_ranks = rank_summary (6 places, 4 rankings)", land, N, P)
    M = RS[(RS.engine == "median of 3") & (RS.stratum != "pooled")].copy()
    M["label"] = M.stratum.map(lambda s: BLOCK.get(s, s))
    Q = M.set_index(["label", "set", "statistic"]).value
    S = TS[(TS.landscape == land) & (TS.stratum_kind != "design cell") & TS.statistic.isin(STAT.values())].set_index(["stratum", "set", "statistic"]).value
    add("t_ranks_strata = rank_summary (architecture, scenario)", land, S, Q)
    DS = pd.read_csv(f"{ROOT}/{gal}/detection_stats.tsv", sep="\t")
    DS = DS[DS.label.isin(SET)].assign(set=lambda d: d.label.map(SET))
    for lab, col in (("True positives (causal or linked SNPs)", "TP_median"), ("False positives (neutral SNPs)", "FP_median")):
        a = TA[(TA.landscape == land) & (TA.label == lab)].set_index(["architecture", "set"])["median"]
        b = DS[DS.arch.isin(ARCH)].set_index(["arch", "set"])[col]
        add(f"t_detection_arch = detection_stats ({col})", land, a, b)
        a = TD[(TD.landscape == land) & (TD.label == lab)].set_index("set")["median"]
        b = DS[DS.arch == "pooled"].set_index("set")[col]
        add(f"t_detection = detection_stats ({col})", land, a, b)
    CD = pd.read_csv(f"{ROOT}/{gal}/C_delta_stats.tsv", sep="\t")
    CD = CD[CD.label.isin(SET) & (CD.arch == "pooled")].assign(set=lambda d: d.label.map(SET)).set_index("set")["median"]
    a = TO[(TO.landscape == land) & (TO.offset_method == "Median of the three") & (TO.set != "Causal loci")].set_index("set").median_difference
    add("t_offset_accuracy = C_delta_stats (median difference)", land, a, CD)
    B4 = pd.read_csv(f"{ROOT}/{gal}/B4_aucpr_stats.tsv", sep="\t")
    B4 = B4.set_index(["label", "arch"])["median"]
    X = S1[S1.landscape == land]
    for key, col in (("LFMM", "median_aucpr_lfmm2"), ("EMMAX", "median_aucpr_emmax"), ("RDA", "median_aucpr_rda_uncorrected")):
        a = X.set_index("architecture")[col]
        b = B4.xs(key, level="label")
        add(f"fig_s1_numbers = B4_aucpr_stats ({key})", land, a, b)
R = pd.DataFrame(rows)
R.to_csv(f"{T}/gate_new_landscapes.tsv", sep="\t", index=False)
print(R.to_string(index=False))
print(f"{len(R)} checks, {int((R.n != R.n_agree).sum())} with a difference")
