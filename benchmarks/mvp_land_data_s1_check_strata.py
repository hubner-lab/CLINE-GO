#!/usr/bin/env python3
"""mvp_land_data_s1_check_strata.py -- the checks of Job F5 (data_s1_landscapes_LAB_PROMPT.md): every plotted cell of the
new s_rank_strata_parts_* (data_s1_2026-11/checks/strata_top3.tsv, written by mvp_land_figures.R FIG_MODE=merge) equals
100 minus the published `ranks 4-6` value and the table's `ranks 1-3` value (landscape_ms/tables/t_ranks_strata.tsv,
whose design-cell and scenario rows are median of the three), and the six sets sum to 300 in every cell.
Writes data_s1_2026-11/jobF5_checks.tsv. Run from the repo root."""
import collections, csv, sys
L = "benchmarks/mvp_eval/landscape_ms"; D = f"{L}/data_s1_2026-11"
LAB = {"2/3 methods": "2/3 set", "1/3 methods": "1/3 set", "3/3 methods": "3/3 set", "LFMM": "LFMM2", "RDA": "RDA", "EMMAX": "EMMAX"}
T = {}
for r in csv.DictReader(open(f"{L}/tables/t_ranks_strata.tsv", newline=""), delimiter="\t"):
    if r["stratum_kind"] in ("design cell", "demographic scenario"):
        T[(r["landscape"], r["stratum_kind"], r["stratum"], r["set"], r["statistic"])] = float(r["value"])
P = list(csv.DictReader(open(f"{D}/checks/strata_top3.tsv", newline=""), delimiter="\t"))
d1 = d2 = d3 = 0.0; sums = collections.defaultdict(float); n = 0
for r in P:
    k = (r["landscape"], r["stratum_kind"], r["stratum"].replace("\n", " "), LAB[r["panel"]])
    v = float(r["top3_pct"])
    d1 = max(d1, abs(v - (100 - T[k + ("ranks 4–6",)])))
    d2 = max(d2, abs(v - T[k + ("ranks 1–3",)]))
    d3 = max(d3, abs(float(r["half_pct"]) - T[k + ("ranks 4–6",)]))
    sums[k[:3]] += v; n += 1
dsum = max(abs(s - 300) for s in sums.values())
cells = len(sums)
rows = [("plotted cells", "", n, ""),
        ("top3 = 100 - published ranks 4-6", "max |diff|", f"{d1:.3g}", "OK" if d1 < 1e-9 else "FAIL"),
        ("top3 = table ranks 1-3", "max |diff|", f"{d2:.3g}", "OK" if d2 < 1e-9 else "FAIL"),
        ("published key's ranks 4-6 (redrawn) = table ranks 4-6", "max |diff|", f"{d3:.3g}", "OK" if d3 < 1e-9 else "FAIL"),
        (f"six sets sum to 300 in each of {cells} cells", "max |sum - 300|", f"{dsum:.3g}", "OK" if dsum < 1e-9 else "FAIL")]
assert n == 3 * 6 * (12 + 5), n
with open(f"{D}/jobF5_checks.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["check", "measure", "value", "result"]); w.writerows(rows)
for r in rows: print(*r, sep="\t")
sys.exit(1 if any(r[3] == "FAIL" for r in rows) else 0)
