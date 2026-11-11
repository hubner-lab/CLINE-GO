#!/usr/bin/env python3
"""j4_readings.py <s_random_climate.svg> <s_random_sets.svg> <random_loo.tsv> <random_climate_tau.tsv> <out tsv>
Job 4's full-resolution readings, from the SVG text nodes of the two members (byte-identical to the published ones)
against the tables they were drawn from: badge numbers, colour-bar titles, the tau over each rank plot, tau and p in
each histogram. The histogram p is recomputed from random_loo.tsv as (1 + #{random tau <= adaptive tau}) / (n + 1)."""
import csv, re, sys
from collections import Counter
clim, sets, loo, tau, out = sys.argv[1:6]
def texts(p):
    return [t.replace("&#39;", "'") for t in re.findall(r"<text[^>]*>([^<]*)</text>", open(p).read())]
T18, T17 = texts(clim), texts(sets)
rows = []
nums = Counter(t for t in T18 if re.fullmatch(r"\d+", t))
for k in range(1, 12):
    rows.append(("badge/axis count", str(k), str(nums[str(k)]), "", ""))
rows.append(("colour-bar titles", "Consensus / rank", str(T18.count("Consensus")) + " / " + str(T18.count("rank")), "3 / 3 expected (three offset maps)", ""))
for e, lab in (("gf", "Gradient Forest"), ("geometric", "Geometric offset"), ("rda", "RDA-based offset")):
    r = next(x for x in csv.DictReader(open(tau), delimiter="\t") if x["engine"] == e)
    want = f"Kendall's τ = {float(r['med']):.2f}"
    rows.append(("rank-plot tau", lab, want, "printed" if want in T18 else "NOT PRINTED", f"random_climate_tau.tsv med {r['med']}"))
    L = [x for x in csv.DictReader(open(loo), delimiter="\t") if x["method"] == e]
    a = float(L[0]["adaptive"]); n = len(L); k = sum(float(x["tau"]) <= a + 1e-12 for x in L)
    p = (1 + k) / (n + 1)
    wt = f"Kendall's τ = {a:.2f}".replace("-", "−"); wp = f"p = {p:.3f}"
    rows.append(("histogram tau / p", lab, f"{wt} / {wp}", ("printed" if wt in T17 else "TAU NOT PRINTED") + " / " + ("printed" if wp in T17 else "P NOT PRINTED"),
                 f"random_loo.tsv adaptive {a:.6f}, {k} of {n} random sets at or below"))
with open(out, "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["reading", "item", "value", "status", "source"]); w.writerows(rows)
for r in rows: print("\t".join(r))
print("all texts with 'shift' or 'rank' or 'τ':", [t for t in T18 if re.search("shift|rank|τ|SD|Consensus", t)])
