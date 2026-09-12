"""Shared helpers for the Python suite.

Two DIFFERENT import mechanisms are needed and one helper cannot cover both:

  * ``scripts/`` is not a package (no ``__init__.py``), so each script is loaded
    by absolute path with ``importlib.util.spec_from_file_location``.
    ``load_script()`` does that.
  * ``workflow/methods/`` IS a package (0-byte ``__init__.py``) and ``gwas.py``
    does ``from .gea import GEA_METHODS``. A file-path load breaks that relative
    import, so the package needs ``workflow/`` on ``sys.path`` and a real
    ``import``. ``import_methods()`` does that.

``load_script()`` must NOT be pointed at ``scripts/gff2topr.py``: that file has no
``if __name__`` guard and calls ``sys.exit(1)`` from its module body, so importing
it aborts the interpreter. It is exercised as a subprocess instead.
"""
import importlib
import importlib.util
import os
import random
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPTS = os.path.join(REPO_ROOT, "scripts")


def load_script(name):
    """Import scripts/<name>.py as a module object, by absolute path."""
    path = os.path.join(SCRIPTS, name + ".py")
    spec = importlib.util.spec_from_file_location("clinego_" + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def script_path(name):
    """Absolute path to scripts/<name>.py, for subprocess invocation."""
    return os.path.join(SCRIPTS, name + ".py")


def import_methods(name):
    """Import workflow/methods/<name> as a package submodule."""
    workflow_dir = os.path.join(REPO_ROOT, "workflow")
    if workflow_dir not in sys.path:
        sys.path.insert(0, workflow_dir)
    return importlib.import_module("methods." + name)


# --------------------------------------------------------------- fixtures
# Lifted verbatim (bar this comment) from benchmarks/verify_design_adequacy.py
# :26 and :79, which builds exactly these designs but only PRINTS the result —
# it asserts nothing and always exits 0. Both builders already take the target
# directory as their first argument, so the tests own the tempdir.

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


def make_messy_fixture(dirpath):
    """The paths a clean fixture never reaches, all in one design.

    site S03 — every row has an empty bio_2 -> the whole site drops out
    site S00 — one row has an empty bio_2  -> that site's n must be 3, not 4
    one climate row (ORPHAN) has no metadata entry
    metadata carries 2 samples with no climate row at all
    bio_12 is invariant across sites
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


def read_table(path):
    """Read a design_adequacy output TSV into {metric: (value, flag, note)}."""
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = (f[1], f[2], f[3] if len(f) > 3 else "")
    return out
