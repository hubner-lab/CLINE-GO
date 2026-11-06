#!/usr/bin/env python3
"""k_evidence.py: what the sNMF runs say about K beyond the cross-entropy curve.

The pipeline's plot draws the best run per K, and its table lists every run but says nothing about
the ancestry matrices. Where cross-entropy has no minimum, K has to be argued from the runs
themselves, so this reads all Q matrices of mode=prestructure and asks three things per K:

  the curve         mean and spread of cross-entropy, the step from K-1, how far the step falls
  reproducibility   do independent runs find the same clusters, column-matched by correlation
  composition       how many accessions a cluster holds, where they come from, what K adds to K-1

Columns are matched with the Hungarian algorithm on Pearson r. A cluster of the best run counts as
replicated when its match reaches r >= 0.9 in at least 8 of the other 9 runs. Two runs share a
solution when every cluster matches at r >= 0.9, and solutions are the connected groups of runs.

One repetition of LEA::snmf() uses one seed for every K, and that seed also draws its masked
genotypes, so a run's cross-entropy carries a mask effect across K: the same run can be lowest at
most K with no better fit behind it. Solutions are compared by mean cross-entropy, with mask_effect
removed: the run's mean deviation over K 2 to 6, where all runs reach one solution.

A sample belongs to a cluster at more than 60 % of its ancestry, the 1001 Genomes Consortium's rule
(Alonso-Blanco et al. 2016); the rest are admixed. Its main cluster is the largest column, which is
how the pipeline's per-cluster LD decay assigns it (scripts/ld_decay_prepare.R).

Stops if the sample order of metadata.tsv differs from the pruned .fam, or if a best-run Q matrix
differs from the pipeline's clusters_K{k}.tsv.

Writes under ANALYSIS/Arabidopsis1135/k_evidence/:
  curve.tsv              K, runs, mean, sd, min, max, step, step_over_sd, step_ratio
  reproducibility.tsv    K, median over run pairs of mean and worst r and of samples keeping their
                         main cluster; clusters of the best run replicated
  composition.tsv        K, share over 60 %, smallest main cluster, clusters with >= 10 samples
  solutions.tsv          K, solution, runs in it, their mean cross-entropy raw and without the
                         mask effect, whether it holds the best run
  clusters.tsv           K, cluster, n_main, n_over_60, runs replicated, centre, countries
  added_by_K.tsv         K, the cluster with no partner at K-1, and how well the others persisted
  assignments.tsv        sample, country, latitude, longitude, one main-cluster column per K

Run from the repo root, in the image:
  docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -v $PWD:/data -w /data \
    cline-go:6cdedcc python3 -B manuscript/arabidopsis/k_evidence.py
"""
import itertools
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

PROJECT = "Arabidopsis1135"
RESULTS = pathlib.Path(f"{PROJECT}_results")
WORK = RESULTS / "_work/maf0.05_miss0.1_smiss0.5_rel0.998/ld0.2_win100_step20"
SNMF = WORK / "Arabidopsis_1001G.snmf"
FAM = WORK / "Arabidopsis_1001G.fam"
TABLES = RESULTS / "PreStructure/tables"
METADATA = RESULTS / "Processing/tables/metadata.tsv"
SOURCE_METADATA = pathlib.Path("data/Arabidopsis_1001G_metadata.tsv.bak")
OUT = pathlib.Path("ANALYSIS/Arabidopsis1135/k_evidence")

K_VALUES = range(2, 13)
STABLE_K = range(2, 7)
RUNS = range(1, 11)
REPLICATED_R = 0.9
REPLICATED_RUNS = 8
MEMBER_SHARE = 0.6
LD_MIN_SAMPLES = 10


def read_q(k, run):
    return np.loadtxt(SNMF / f"K{k}" / f"run{run}" / f"Arabidopsis_1001G_r{run}.{k}.Q")


def cross_correlation(a, b):
    """Pearson r between each column of a and each column of b."""
    n_a = a.shape[1]
    r = np.corrcoef(a, b, rowvar=False)[:n_a, n_a:]
    return np.nan_to_num(r)


def match(a, b):
    """Columns of b matched to columns of a: returns the b column per a column and their r."""
    r = cross_correlation(a, b)
    rows, cols = linear_sum_assignment(-r)
    partner = np.full(a.shape[1], -1)
    partner[rows] = cols
    return partner, r[rows, cols]


def same_main(a, b, partner):
    """Share of samples whose main cluster in a is the partner of their main cluster in b."""
    return np.mean(partner[a.argmax(axis=1)] == b.argmax(axis=1))


def solutions(q):
    """Connected groups of runs in which every cluster matches at r >= REPLICATED_R."""
    group = {run: run for run in q}

    def find(run):
        while group[run] != run:
            run = group[run]
        return run

    for a, b in itertools.combinations(q, 2):
        _, r = match(q[a], q[b])
        if r.min() >= REPLICATED_R:
            group[find(a)] = find(b)
    found = {}
    for run in q:
        found.setdefault(find(run), []).append(run)
    return sorted(found.values(), key=len, reverse=True)


def top_countries(countries, n=3):
    counts = countries.value_counts()
    return ", ".join(f"{name} {count}" for name, count in counts.head(n).items())


def check_sample_order(meta):
    fam = pd.read_csv(FAM, sep=r"\s+", header=None, dtype=str)
    if list(fam[1]) != list(meta["sample"]):
        sys.exit(f"k_evidence: sample order differs between {FAM} and {METADATA}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(METADATA, sep="\t", dtype={"site": str, "sample": str})
    check_sample_order(meta)
    source = pd.read_csv(SOURCE_METADATA, sep="\t", dtype={"site": str, "sample": str})
    country = dict(zip("S" + source["sample"], source["site"]))
    meta["country"] = meta["sample"].map(country)
    if meta["country"].isna().any():
        sys.exit("k_evidence: a sample has no country in the source metadata")

    ce = pd.read_csv(TABLES / f"cross_entropy_K{K_VALUES[0]}-{K_VALUES[-1]}.tsv", sep="\t")
    curve = ce.groupby("K")["cross_entropy"].agg(
        runs="count", mean="mean", sd="std", min="min", max="max").reset_index()
    pooled_sd = np.sqrt((curve["sd"] ** 2).mean())
    curve["step"] = -curve["mean"].diff()
    curve["step_over_pooled_sd"] = curve["step"] / pooled_sd
    curve["step_ratio"] = curve["step"] / curve["step"].shift(-1)
    curve.to_csv(OUT / "curve.tsv", sep="\t", index=False)

    wide = ce.pivot(index="run", columns="K", values="cross_entropy")
    mask_effect = (wide - wide.mean()).loc[:, list(STABLE_K)].mean(axis=1)

    reproducibility, composition, clusters, added, solution_rows = [], [], [], [], []
    assignments = meta[["sample", "country", "latitude", "longitude"]].copy()
    previous_best = None
    for k in K_VALUES:
        q = {run: read_q(k, run) for run in RUNS}
        best_run = int(ce[(ce["K"] == k) & ce["best_run"]]["run"].iloc[0])
        best = q[best_run]

        table = pd.read_csv(TABLES / f"K{k}" / f"clusters_K{k}.tsv", sep="\t",
                            dtype={"site": str, "sample": str})
        if list(table["sample"]) != list(meta["sample"]):
            sys.exit(f"k_evidence: clusters_K{k}.tsv is not in metadata order")
        pipeline_q = table[[f"C{i}" for i in range(1, k + 1)]].to_numpy()
        if np.abs(pipeline_q - best).max() > 1e-4:
            sys.exit(f"k_evidence: run {best_run} is not the Q matrix of clusters_K{k}.tsv")

        mean_r, worst_r, kept = [], [], []
        for a, b in itertools.combinations(RUNS, 2):
            partner, r = match(q[a], q[b])
            mean_r.append(r.mean())
            worst_r.append(r.min())
            kept.append(same_main(q[a], q[b], partner))

        ce_k = ce[ce["K"] == k].set_index("run")["cross_entropy"]
        for number, runs in enumerate(solutions(q), start=1):
            solution_rows.append({
                "K": k, "solution": number, "n_runs": len(runs),
                "runs": ",".join(map(str, runs)), "cross_entropy": ce_k[runs].mean(),
                "cross_entropy_no_mask": (ce_k[runs] - mask_effect[runs]).mean(),
                "holds_best_run": best_run in runs})

        hits = np.zeros(k, dtype=int)
        for run in RUNS:
            if run != best_run:
                _, r = match(best, q[run])
                hits += r >= REPLICATED_R
        replicated = hits >= REPLICATED_RUNS

        main_cluster = best.argmax(axis=1)
        over = best.max(axis=1) > MEMBER_SHARE
        n_main = np.bincount(main_cluster, minlength=k)
        reproducibility.append({
            "K": k, "median_mean_r": np.median(mean_r), "median_worst_r": np.median(worst_r),
            "median_same_main": np.median(kept), "replicated": int(replicated.sum()),
            "best_run": best_run})
        composition.append({
            "K": k, "share_over_60": over.mean(), "smallest_main": int(n_main.min()),
            "clusters_ge_10": int((n_main >= LD_MIN_SAMPLES).sum())})
        for i in range(k):
            members = meta[main_cluster == i]
            clusters.append({
                "K": k, "cluster": f"C{i + 1}", "n_main": int(n_main[i]),
                "n_over_60": int((over & (main_cluster == i)).sum()),
                "runs_replicated": int(hits[i]),
                "latitude": members["latitude"].mean(), "longitude": members["longitude"].mean(),
                "countries": top_countries(members["country"])})
        assignments[f"K{k}"] = [f"C{i + 1}" for i in main_cluster]

        if previous_best is not None:
            partner, r = match(previous_best, best)
            new = sorted(set(range(k)) - set(partner.tolist()))[0]
            added.append({
                "K": k, "new_cluster": f"C{new + 1}", "n_main": int(n_main[new]),
                "n_over_60": int((over & (main_cluster == new)).sum()),
                "countries": top_countries(meta[main_cluster == new]["country"], 4),
                "previous_mean_r": r.mean(), "previous_worst_r": r.min()})
        previous_best = best

    pd.DataFrame(reproducibility).to_csv(OUT / "reproducibility.tsv", sep="\t", index=False)
    pd.DataFrame(composition).to_csv(OUT / "composition.tsv", sep="\t", index=False)
    pd.DataFrame(solution_rows).to_csv(OUT / "solutions.tsv", sep="\t", index=False)
    pd.DataFrame(clusters).to_csv(OUT / "clusters.tsv", sep="\t", index=False)
    pd.DataFrame(added).to_csv(OUT / "added_by_K.tsv", sep="\t", index=False)
    assignments.to_csv(OUT / "assignments.tsv", sep="\t", index=False)
    print(f"k_evidence: pooled run-to-run sd {pooled_sd:.5f}; tables in {OUT}")


if __name__ == "__main__":
    main()
