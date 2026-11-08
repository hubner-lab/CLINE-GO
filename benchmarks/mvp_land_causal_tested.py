#!/usr/bin/env python3
"""mvp_land_causal_tested.py -- Job D of docs/manuscript/lab/data_s1_landscapes_LAB_PROMPT.md.

Methods says every causal locus with a MAF of at least 0.01 is among the tested SNPs. Count it per replicate:
  causal loci  = rows of data/mvp/MVP<seed>/loci_freq.tsv with causal TRUE and maf >= 0.01 (the deposit's MAF of the
                 sampled individuals); cross-checked against the `causal` rows of truth_any.tsv (the truth set every
                 scorer reads) and the manifest's n_causal_maf01 (benchmarks/mvp_seeds.tsv)
  tested SNPs  = the (chr, pos) rows of the replicate's GEA p-value tables with a p-value in some trait column. LFMM,
                 EMMAX and RDA must test the same SNPs (asserted), so one set per replicate.
Trees: SS-Mtn and Est-Clines MVP<seed>hd_results at the repo root; SS-Clines (control) _archive/mvp_ssclines/results/.
Host python, stdlib only. Read-only on every tree; writes OUT (default benchmarks/mvp_eval/landscape_ms/data_s1_2026-11)
causal_tested.tsv and causal_tested_missing.tsv (one row per missing locus, empty when none).
Run from the repo root:  python3 -I benchmarks/mvp_land_causal_tested.py [OUT]
"""
import csv
import glob
import math
import os
import sys

ROOT = os.getcwd()
assert os.path.isfile(os.path.join(ROOT, "Snakefile")), "run from the repo root"
OUT = sys.argv[1] if len(sys.argv) > 1 else "benchmarks/mvp_eval/landscape_ms/data_s1_2026-11"
os.makedirs(OUT, exist_ok=True)


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


MAN = {r["seed"]: r for r in rows("benchmarks/mvp_seeds.tsv")}
LAND = {
    "SS-Clines": sorted(r["seed"] for r in MAN.values() if r["arm"] == "primary_ssclines"),
    "SS-Mtn": open("work/mvp_land/seeds_ssmtn.txt").read().split(),
    "Est-Clines": open("work/mvp_land/seeds_estclines.txt").read().split(),
}
for k, v in LAND.items():
    assert len(v) == 600 == len(set(v)), (k, len(v))


def tree(land, s):
    if land == "SS-Clines":
        cand = sorted(glob.glob(f"_archive/mvp_ssclines/results/MVP{s}hd_results")) or \
               sorted(glob.glob(f"_archive/mvp_ssclines/results/MVP{s}_c1_results"))
        assert cand, f"no SS-Clines tree for {s}"
        return cand[0]
    return f"MVP{s}hd_results"


def tested(t, s):
    sets = {}
    for m in ("LFMM", "EMMAX", "RDA"):
        fs = [f for f in glob.glob(f"{t}/GEA/tables/methods/{m}/{m}_pvalues_K*.tsv") if "_sig_" not in f]
        assert len(fs) == 1, (t, m, fs)
        keys = set()
        with open(fs[0], newline="") as f:
            r = csv.reader(f, delimiter="\t")
            head = next(r)
            ic, ip = head.index("chr"), head.index("pos")
            tc = [i for i, h in enumerate(head) if h not in ("SNPID", "chr", "pos")]
            assert tc, (fs[0], head)
            for row in r:
                if any(row[i] not in ("", "NA") for i in tc):
                    keys.add((row[ic], row[ip]))
        sets[m] = keys
    assert sets["LFMM"] == sets["EMMAX"] == sets["RDA"], f"{t}: the three methods test different SNPs"
    return sets["LFMM"]


summary, missing = [], []
for land, seeds in LAND.items():
    n_rep_miss = n_causal = n_miss = 0
    for s in seeds:
        d = f"data/mvp/MVP{s}"
        lf = rows(f"{d}/loci_freq.tsv")
        causal = {(r["chr"], r["pos"]) for r in lf if r["causal"] == "TRUE" and float(r["maf"]) >= 0.01}
        truth = {(r["chr"], r["pos"]) for r in rows(f"{d}/truth_any.tsv") if r["category"] == "causal"}
        assert causal == truth, f"{s}: loci_freq causal at MAF >= 0.01 is not truth_any's causal set"
        assert len(causal) == int(MAN[s]["n_causal_maf01"]), (s, len(causal), MAN[s]["n_causal_maf01"])
        T = tested(tree(land, s), s)
        miss = sorted(causal - T, key=lambda x: (int(x[0]), int(x[1])))
        n_causal += len(causal)
        n_miss += len(miss)
        n_rep_miss += bool(miss)
        maf = {(r["chr"], r["pos"]): r["maf"] for r in lf}
        missing += [(land, s, c, p, maf[(c, p)]) for c, p in miss]
    summary.append((land, len(seeds), n_rep_miss, n_causal, n_miss))
    print(f"{land}: {len(seeds)} replicates, {n_rep_miss} with a missing causal locus; "
          f"{n_causal} causal loci at MAF >= 0.01, {n_miss} not tested", flush=True)

with open(f"{OUT}/causal_tested.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["landscape", "replicates", "replicates_with_missing", "causal_loci", "causal_loci_missing"])
    w.writerows(summary)
with open(f"{OUT}/causal_tested_missing.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["landscape", "seed", "chr", "pos", "maf"])
    w.writerows(missing)
print("done")
