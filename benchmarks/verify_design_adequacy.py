#!/usr/bin/env python3
"""Exercise scripts/design_adequacy.py on synthetic fixtures.

Two designs, both written in the pipeline's own file formats:

  A) "Trifolium-shaped": 58 samples at 11 sites, 18 bioclim predictors driven by
     3 latent gradients. This is the design the dossier records as producing a
     rank-10 environment covariance with condition number ~1e6 and a geometric
     offset that put 99.85% of its mass on environment PC10.

  B) "Adequate": 300 samples at 60 sites, 4 decorrelated predictors.

Run:  python3 benchmarks/verify_design_adequacy.py
"""
import math
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "design_adequacy.py")


def make_fixture(dirpath, n_sites, n_pred, per_site, n_latent, noise, seed,
                 singleton_sites=0):
    rnd = random.Random(seed)
    sites = [f"S{i:02d}" for i in range(n_sites)]
    preds = [f"bio_{i+1}" for i in range(n_pred)]

    # Latent environmental gradients over sites, then predictors as noisy
    # linear combinations of them (how real bioclim variables actually behave).
    latent = [[rnd.gauss(0, 1) for _ in range(n_latent)] for _ in sites]
    load = [[rnd.gauss(0, 1) for _ in range(n_latent)] for _ in preds]
    vals = {}
    for si, st in enumerate(sites):
        row = []
        for pj in range(n_pred):
            v = sum(load[pj][k] * latent[si][k] for k in range(n_latent))
            row.append(v + rnd.gauss(0, noise))
        vals[st] = row

    meta = os.path.join(dirpath, "metadata.tsv")
    clim = os.path.join(dirpath, "climate_present_site.tsv")
    with open(meta, "w") as mh, open(clim, "w") as ch:
        mh.write("site\tsample\tlatitude\tlongitude\n")
        ch.write("sample\t" + "\t".join(preds) + "\n")
        k = 0
        for si, st in enumerate(sites):
            n_here = 1 if si < singleton_sites else per_site
            for _ in range(n_here):
                smp = f"ID{k:04d}"
                k += 1
                mh.write(f"{st}\t{smp}\t{30 + si * 0.1:.4f}\t{35 + si * 0.1:.4f}\n")
                ch.write(smp + "\t" + "\t".join(f"{v:.6f}" for v in vals[st]) + "\n")
    return meta, clim


def run(tag, **kw):
    d = tempfile.mkdtemp(prefix="da_")
    meta, clim = make_fixture(d, **kw)
    out = os.path.join(d, "design_adequacy.tsv")
    r = subprocess.run([sys.executable, SCRIPT, meta, clim, "all", out],
                       capture_output=True, text=True)
    print(f"===== {tag} =====")
    print(r.stderr.rstrip())
    if r.returncode != 0:
        print(f"!! exit {r.returncode}")
        print(r.stdout)
        return
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            print(f"  {f[0]:<32s} {f[1]:>10s}  {f[2]}")
    print()


def make_messy_fixture(dirpath):
    """C) The paths a clean fixture never reaches, all in one design.

    Real climate tables lose rows: download_climate_present.R drops samples whose
    raster cell is NA and records them in climate_na_excluded.tsv. So the climate
    table can (a) be missing samples the metadata has, (b) carry an empty cell,
    (c) carry a sample the metadata does not know, and (d) have one whole site
    whose every row is unusable. Each of those must be COUNTED, and must not
    silently change a per-site n.

    Design: 6 sites x 4 samples, 3 varying predictors + 1 invariant one.
      site S03 — every row has an empty bio_2 -> the whole site drops out
      site S00 — one row has an empty bio_2  -> that site's n must be 3, not 4
      one climate row (ORPHAN) has no metadata entry
      metadata carries 2 samples with no climate row at all
    """
    rnd = random.Random(7)
    sites = [f"S{i:02d}" for i in range(6)]
    preds = ["bio_1", "bio_2", "bio_3", "bio_12"]
    vals = {st: [rnd.gauss(0, 1), rnd.gauss(0, 1), rnd.gauss(0, 1), 5.0] for st in sites}

    meta = os.path.join(dirpath, "metadata.tsv")
    clim = os.path.join(dirpath, "climate_present_site.tsv")
    with open(meta, "w") as mh, open(clim, "w") as ch:
        mh.write("site\tsample\tlatitude\tlongitude\n")
        ch.write("sample\t" + "\t".join(preds) + "\n")
        k = 0
        for si, st in enumerate(sites):
            for row_i in range(4):
                smp = f"ID{k:04d}"
                k += 1
                mh.write(f"{st}\t{smp}\t{30 + si * 0.1:.4f}\t{35 + si * 0.1:.4f}\n")
                v = list(vals[st])
                if st == "S03" or (st == "S00" and row_i == 0):
                    v[1] = ""            # empty cell, as fread would leave an NA
                if st == "S05" and row_i >= 2:
                    continue             # in metadata, absent from climate
                ch.write(smp + "\t" + "\t".join(str(x) if x == "" else f"{x:.6f}"
                                                for x in v) + "\n")
        ch.write("ORPHAN\t" + "\t".join(f"{x:.6f}" for x in vals["S01"]) + "\n")
    return meta, clim


def run_fixture(tag, builder):
    d = tempfile.mkdtemp(prefix="da_")
    meta, clim = builder(d)
    out = os.path.join(d, "design_adequacy.tsv")
    r = subprocess.run([sys.executable, SCRIPT, meta, clim, "all", out],
                       capture_output=True, text=True)
    print(f"===== {tag} =====")
    print(r.stderr.rstrip())
    if r.returncode != 0:
        print(f"!! exit {r.returncode}")
        print(r.stdout)
        return
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            print(f"  {f[0]:<32s} {f[1]:>10s}  {f[2]}")
    print()


if __name__ == "__main__":
    run("A) Trifolium-shaped: 58 samples / 11 sites / 18 predictors",
        n_sites=11, n_pred=18, per_site=6, n_latent=3, noise=0.05, seed=1,
        singleton_sites=2)
    run("B) Adequate: 300 samples / 60 sites / 4 decorrelated predictors",
        n_sites=60, n_pred=4, per_site=5, n_latent=4, noise=1.6, seed=2)
    run_fixture("C) Messy: NA cells, an orphan climate row, an invariant predictor",
                make_messy_fixture)
