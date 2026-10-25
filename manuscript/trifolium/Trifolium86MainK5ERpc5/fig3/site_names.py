#!/usr/bin/env python3
"""Short site names for the whole Trifolium manuscript (plan ~/.claude/plans/fancy-scribbling-garden.md, Phase 1).

Rule (user, 2026-10-16): sites are numbered 1..11 from south to north by latitude. Writes
  fig3/site_names.tsv                      site_id, short, latitude, longitude  (single source for every figure + text)
  fig3/inputs/clusters_K5.short.tsv        clusters_K5.tsv with `site` replaced by the short name
  fig3/inputs/metadata_climate.short.tsv   metadata_climate.tsv with `site` replaced by the short name
The two input copies let the derived pipeline plot scripts label sites by short name without any code change.
Run from the repo root: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.py
"""
import csv
import os
import sys

R = "Trifolium86MainK5ERpc5_results"
F = "ANALYSIS/Trifolium86MainK5ERpc5/fig3"
META = f"{R}/_intermediate/samples/metadata_climate.tsv"
CLUS = f"{R}/PreStructure/tables/K5/clusters_K5.tsv"


def read(path):
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return rd.fieldnames, list(rd)


mcols, meta = read(META)
coords = {}
for r in meta:
    c = (float(r["latitude"]), float(r["longitude"]))
    if coords.setdefault(r["site"], c) != c:
        sys.exit(f"site {r['site']} has more than one coordinate")
lats = [c[0] for c in coords.values()]
if len(set(lats)) != len(lats):
    sys.exit("latitude tie: the south-to-north rule needs a tiebreak")
order = sorted(coords, key=lambda s: coords[s][0])
short = {s: str(i) for i, s in enumerate(order, 1)}

os.makedirs(f"{F}/inputs", exist_ok=True)
with open(f"{F}/site_names.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["site_id", "short", "latitude", "longitude"])
    for s in order:
        w.writerow([s, short[s], coords[s][0], coords[s][1]])

ccols, clus = read(CLUS)
if {r["site"] for r in clus} != set(short):
    sys.exit("clusters_K5.tsv and metadata_climate.tsv name different sites")
for path, cols, rows in ((f"{F}/inputs/clusters_K5.short.tsv", ccols, clus),
                         (f"{F}/inputs/metadata_climate.short.tsv", mcols, meta)):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({**r, "site": short[r["site"]]})
print(" ".join(f"{s}={short[s]}" for s in order))
