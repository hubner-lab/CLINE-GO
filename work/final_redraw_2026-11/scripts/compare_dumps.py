#!/usr/bin/env python3
"""compare_dumps.py <ref dir> <new dir> <out tsv> [refmember=newmember ...] -- the value check of the final redraw.

Compares two dump_layers.R dumps layer by layer: the reference run (the starting copy, reproducing the published
figure) against the redrawn one. Members are matched by name, or by the given refmember=newmember pairs (a renamed or
re-lettered panel). For every matched layer: same row count; every numeric column compared cell by cell (largest
absolute deviation; NA must match NA); every text column compared cell by cell (cells that differ are counted, with
up to three examples, so a relabel shows up as exactly the strings the job changes). Members present on one side only
are listed. Writes one row per layer and prints the overall largest numeric deviation. Stdlib only.
"""
import csv
import math
import sys
from pathlib import Path

ref, new, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
pairs = dict(a.split("=", 1) for a in sys.argv[4:])


def index(d):
    with open(d / "objects.tsv", newline="") as f:
        return {(r["member"], int(r["layer"])): r for r in csv.DictReader(f, delimiter="\t")}


def table(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    return (rows[0], rows[1:]) if rows else ([], [])


def num(x):
    if x in ("", "NA"):
        return None
    try:
        return float(x)
    except ValueError:
        return "text"


R, N = index(ref), index(new)
rows, worst = [], 0.0
matched_new = set()
for (m, layer), r in sorted(R.items()):
    m2 = pairs.get(m, m)
    n = N.get((m2, layer))
    if n is None:
        rows.append([m, m2, layer, r["geom"], r["rows"], "", "MISSING in new", "", "", ""])
        continue
    matched_new.add((m2, layer))
    h1, d1 = table(ref / r["file"])
    h2, d2 = table(new / n["file"])
    status, dev, txt_diff, examples = "ok", 0.0, [], []
    if len(d1) != len(d2):
        status = f"rows {len(d1)} vs {len(d2)}"
    cols = [c for c in h1 if c in h2]
    only = sorted(set(h1) ^ set(h2))
    if only:
        status += f"; columns on one side: {','.join(only)}"
    for c in cols:
        i1, i2 = h1.index(c), h2.index(c)
        diffs = 0
        for a, b in zip(d1, d2):
            x, y = num(a[i1]), num(b[i2])
            if x is None or y is None or x == "text" or y == "text":
                if a[i1] != b[i2]:
                    diffs += 1
                    if len(examples) < 3:
                        examples.append(f"{c}: {a[i1]!r} -> {b[i2]!r}")
            else:
                if math.isinf(x) or math.isinf(y):
                    d = 0.0 if x == y else math.inf
                else:
                    d = abs(x - y)
                dev = max(dev, d)
        if diffs:
            txt_diff.append(f"{c} ({diffs})")
    worst = max(worst, dev)
    if dev > 0:
        status = "NUMERIC DIFF" if status == "ok" else status + "; NUMERIC DIFF"
    rows.append([m, m2, layer, r["geom"], r["rows"], n["rows"], status, f"{dev:.3g}", "; ".join(txt_diff),
                 " | ".join(examples)])
for (m2, layer), n in sorted(N.items()):
    if (m2, layer) not in matched_new:
        rows.append(["", m2, layer, n["geom"], "", n["rows"], "NEW layer (no reference)", "", "", ""])
with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["ref_member", "new_member", "layer", "geom", "ref_rows", "new_rows", "status", "max_abs_dev",
                "text_columns_differing", "examples"])
    w.writerows(rows)
bad = [r for r in rows if r[6] not in ("ok", "NEW layer (no reference)")]
print(f"{out.name}: {len(rows)} layers, largest numeric deviation {worst:.3g}, not ok: {len(bad)}")
for r in bad:
    print("  ", r[0], r[1], r[2], r[6])
