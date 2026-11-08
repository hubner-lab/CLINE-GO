#!/usr/bin/env python3
"""mvp_land_data_s1_check_tables.py -- the checks of Job E (docs/manuscript/lab/data_s1_landscapes_LAB_PROMPT.md).

Compares the rerun of benchmarks/mvp_land_tables.R (NEW = data_s1_2026-11/tables) with the published landscapes tables
(PUB = landscape_ms/tables) as strings, cell by cell:
  t_ranks_omnibus   every cell equal
  t_ranks           value, lo, hi, n equal; NEW has no p / p_holm column
  t_ranks.md        no p in a set's row; every ranking ends with the whole-table line
  numbers.tsv       every row equal on the published columns other than p_perm; NEW has no p_perm column, and the
                    published rows that carried a p_perm are listed (dropped by design)
  every other .tsv  equal in every shared column and row (gate_* written by the run are compared too where both exist;
                    gate_holm_change.tsv and gate_new_landscapes.tsv are not written any more by design)
Writes NEW/../jobE_checks.tsv. Run from the repo root: python3 -I benchmarks/mvp_land_data_s1_check_tables.py
"""
import csv
import os
import re
import sys

L = "benchmarks/mvp_eval/landscape_ms"
PUB, NEW = f"{L}/tables", f"{L}/data_s1_2026-11/tables"
OUT = f"{L}/data_s1_2026-11/jobE_checks.tsv"


def read(p):
    with open(p, newline="") as f:
        r = list(csv.reader(f, delimiter="\t"))
    return r[0], r[1:]


checks = []


def add(item, what, result, detail=""):
    checks.append((item, what, result, detail))
    print(f"{result:5s} {item}: {what} {detail}", flush=True)


def compare(name, cols=None, drop=()):
    hp, rp = read(f"{PUB}/{name}")
    hn, rn = read(f"{NEW}/{name}")
    shared = [c for c in hp if c in hn and c not in drop] if cols is None else cols
    ip = [hp.index(c) for c in shared]
    inn = [hn.index(c) for c in shared]
    if len(rp) != len(rn):
        add(name, f"rows {len(rp)} vs {len(rn)}", "FAIL")
        return
    bad = []
    for i, (a, b) in enumerate(zip(rp, rn)):
        for c, x, y in zip(shared, ip, inn):
            if a[x] != b[y]:
                bad.append(f"row {i + 2} {c}: {a[x]} -> {b[y]}")
    add(name, f"{len(rp)} rows x {len(shared)} columns ({', '.join(shared)}) equal as strings",
        "OK" if not bad else "FAIL", "; ".join(bad[:10]) + (f" (+{len(bad) - 10})" if len(bad) > 10 else ""))
    only_p = [c for c in hp if c not in hn]
    only_n = [c for c in hn if c not in hp]
    if only_p or only_n:
        add(name, "columns only in published / only in new", "INFO", f"{only_p} / {only_n}")


# omnibus: every cell
compare("t_ranks_omnibus.tsv")
# t_ranks: value lo hi n, and no p columns
compare("t_ranks.tsv", cols=["landscape", "offset_method", "set", "statistic", "value", "lo", "hi", "n"])
hn, _ = read(f"{NEW}/t_ranks.tsv")
add("t_ranks.tsv", "no p or p_holm column", "OK" if not {"p", "p_holm"} & set(hn) else "FAIL", str(hn))
# t_ranks.md
md = open(f"{NEW}/t_ranks.md").read()
sections = re.split(r"^## ", md, flags=re.M)[1:]
bad_p = [l for l in md.splitlines() if l.startswith("| ") and re.search(r"\bp\s*[=<]|P\s*[=<]|\bp_holm", l)]
add("t_ranks.md", "no p in any table cell", "OK" if not bad_p else "FAIL", "; ".join(bad_p[:3]))
ends = [s.strip().splitlines()[-1].startswith("Whole table") for s in sections]
add("t_ranks.md", f"each of {len(sections)} rankings ends with the whole-table line", "OK" if all(ends) else "FAIL")
# numbers.tsv: all columns but p_perm; rows that carried a p_perm
hp, rp = read(f"{PUB}/numbers.tsv")
hn, rn = read(f"{NEW}/numbers.tsv")
add("numbers.tsv", "no p_perm column", "OK" if "p_perm" not in hn else "FAIL")
ipp = hp.index("p_perm")
with_p = [r for r in rp if r[ipp] != ""]
add("numbers.tsv", "published rows with a p_perm value (the dropped per-set tests)", "INFO", str(len(with_p)))
cols = [c for c in hp if c != "p_perm"]
SUFFIX = ", permutation P (10,000)"          # the per-set test's name in the figure label, dropped with the test
ifg = hp.index("figure")
relabel = sum(SUFFIX in r[ifg] for r in rp)
add("numbers.tsv", f"published figure labels ending in '{SUFFIX}' (the label loses it with the test)", "INFO", str(relabel))
key = lambda h, r: tuple(r[h.index(c)].replace(SUFFIX, "") if c == "figure" else r[h.index(c)] for c in cols)
sp, sn = sorted(key(hp, r) for r in rp), sorted(key(hn, r) for r in rn)
add("numbers.tsv", f"{len(rp)} published rows vs {len(rn)} new, all columns but p_perm, as a multiset, that suffix removed",
    "OK" if sp == sn else "FAIL", "" if sp == sn else f"{len(set(sp) ^ set(sn))} rows differ")
# gate_ssclines: the published gate minus its raw-p row (removed with the per-set tests), every cell
hp, rp = read(f"{PUB}/gate_ssclines.tsv")
hn, rn = read(f"{NEW}/gate_ssclines.tsv")
removed = [r[0] for r in rp if r[0] not in {x[0] for x in rn}]
add("gate_ssclines.tsv", "published gate rows absent from the rerun", "OK" if removed == ["raw per-set p vs Table_S1.tsv (a)"]
    else "FAIL", str(removed))
same = [r for r in rp if r[0] not in removed] == rn and hp == hn
add("gate_ssclines.tsv", f"the other {len(rn)} rows equal the published gate in every cell", "OK" if same else "FAIL")
# every other tsv
done = {"t_ranks_omnibus.tsv", "t_ranks.tsv", "numbers.tsv", "gate_ssclines.tsv"}
for f in sorted(os.listdir(NEW)):
    if f.endswith(".tsv") and f not in done:
        if os.path.exists(f"{PUB}/{f}"):
            compare(f)
        else:
            add(f, "not in the published tables", "INFO")
for f in sorted(os.listdir(PUB)):
    if f.endswith(".tsv") and not os.path.exists(f"{NEW}/{f}"):
        add(f, "published, not written by the rerun", "INFO")
for f in sorted(os.listdir(NEW)):
    if f.endswith(".md") and f != "t_ranks.md":
        same = open(f"{NEW}/{f}").read() == open(f"{PUB}/{f}").read()
        add(f, "byte-identical to the published .md", "OK" if same else "DIFF")
gate = read(f"{NEW}/gate_ssclines.tsv")[1]
nd = sum(int(r[4]) for r in gate)
add("gate_ssclines.tsv", f"{len(gate)} items, cells differing", "OK" if nd == 0 else "FAIL", str(nd))

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["item", "check", "result", "detail"])
    w.writerows(checks)
sys.exit(1 if any(c[2] == "FAIL" for c in checks) else 0)
