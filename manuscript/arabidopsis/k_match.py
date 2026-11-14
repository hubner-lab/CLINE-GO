#!/usr/bin/env python3
"""k_match.py: how one project's sNMF clusters at K map onto another project's, sample by sample.

Iterations of the Arabidopsis section that keep the same accessions (Arabidopsis1135maf01 and Arabidopsis1135it3
share their 902) can be compared cluster by cluster: each sample's main cluster (the largest Q column, as
k_evidence.py assigns it) in one project against its main cluster in the other. Clusters are matched one to one by
the Hungarian algorithm on the counts of shared samples; the agreement is the share of samples whose two main
clusters are matched, and the adjusted Rand index of the two partitions (Hubert & Arabie 1985), which needs no
matching.

Reads  ANALYSIS/<project>/k_evidence/{assignments,clusters}.tsv and the same for <other> (k_evidence.py).
Writes ANALYSIS/<project>/k_evidence/k_match_<other>_K<K>.tsv, one row per cluster of <project>:
  cluster, n_main, matched (the other project's cluster), n_main_other, shared, share_of_cluster, share_of_matched,
  runs_replicated and runs_replicated_other (clusters.tsv), countries and countries_other (clusters.tsv)
and prints the agreement and the adjusted Rand index.

Stops when the two projects' sample sets differ: matching then has no common ground.

Usage, from the repo root, in the image (scipy):
  docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -v $PWD:/data -w /data \
    cline-go:6cdedcc python3 -B manuscript/arabidopsis/k_match.py <project> <other> <K>
"""
import csv
import sys
from collections import Counter
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


def read(path):
    if not path.is_file():
        sys.exit(f"k_match.py: missing {path}")
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ari(pairs):
    """Adjusted Rand index of two labelings given as (a, b) pairs."""
    n = len(pairs)
    sum_ab = sum(comb(c, 2) for c in Counter(pairs).values())
    sum_a = sum(comb(c, 2) for c in Counter(a for a, _ in pairs).values())
    sum_b = sum(comb(c, 2) for c in Counter(b for _, b in pairs).values())
    expected = sum_a * sum_b / comb(n, 2)
    return (sum_ab - expected) / ((sum_a + sum_b) / 2 - expected)


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: k_match.py <project> <other> <K>")
    project, other, k = sys.argv[1], sys.argv[2], sys.argv[3]
    col = f"K{k}"
    kev = Path("ANALYSIS") / project / "k_evidence"
    kev_o = Path("ANALYSIS") / other / "k_evidence"
    a = {r["sample"]: r[col] for r in read(kev / "assignments.tsv")}
    b = {r["sample"]: r[col] for r in read(kev_o / "assignments.tsv")}
    if set(a) != set(b):
        sys.exit(f"k_match.py: sample sets differ ({len(set(a) - set(b))} only in {project}, "
                 f"{len(set(b) - set(a))} only in {other})")
    info = {r["cluster"]: r for r in read(kev / "clusters.tsv") if r["K"] == k}
    info_o = {r["cluster"]: r for r in read(kev_o / "clusters.tsv") if r["K"] == k}

    pairs = [(a[s], b[s]) for s in sorted(a)]
    ca = sorted({x for x, _ in pairs}, key=lambda c: int(c[1:]))
    cb = sorted({y for _, y in pairs}, key=lambda c: int(c[1:]))
    counts = Counter(pairs)
    m = np.array([[counts[(x, y)] for y in cb] for x in ca])
    rows, cols = linear_sum_assignment(-m)
    match = {ca[i]: cb[j] for i, j in zip(rows, cols)}
    n_a = Counter(x for x, _ in pairs)
    n_b = Counter(y for _, y in pairs)

    out = []
    for x in ca:
        y = match.get(x, "")
        shared = counts[(x, y)] if y else 0
        out.append({
            "cluster": x, "n_main": n_a[x], "matched": y, "n_main_other": n_b[y] if y else 0, "shared": shared,
            "share_of_cluster": f"{shared / n_a[x]:.3f}", "share_of_matched": f"{shared / n_b[y]:.3f}" if y else "",
            "runs_replicated": info.get(x, {}).get("runs_replicated", ""),
            "runs_replicated_other": info_o.get(y, {}).get("runs_replicated", ""),
            "countries": info.get(x, {}).get("countries", ""),
            "countries_other": info_o.get(y, {}).get("countries", ""),
        })
    path = kev / f"k_match_{other}_K{k}.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    kept = sum(r["shared"] for r in out)
    print(f"wrote {path} ({len(out)} clusters)")
    print(f"K {k}: {kept} of {len(pairs)} samples ({kept / len(pairs):.1%}) keep a matched main cluster; "
          f"adjusted Rand index {ari(pairs):.3f}")


if __name__ == "__main__":
    main()
