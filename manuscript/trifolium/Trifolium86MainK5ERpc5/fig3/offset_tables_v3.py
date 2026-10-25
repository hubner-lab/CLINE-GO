#!/usr/bin/env python3
"""Offset tables for the merged manuscript pack (plan idempotent-conjuring-noodle.md, Phase 5).

v3 rules (user 2026-10-17): consensus = all three methods agree (mean of the GF, geometric and RDA ranks, NO_LINEAR=1),
no linear-family rank; predictor change vs rank per method (DELTA_CORR=1); sites by short name 1..11 (S->N), raw id kept
in a second column. Reads fig3/offset_minou/paper_pack/ (built by fig3/run_offset_minou.sh) and the v2 tables for the
random-SNP test, which does not depend on the consensus rule (asserted against random_loo.tsv of the v3 build).
Writes manuscript_fig3/tables/offset_Table*.tsv and prints them as markdown for the text.
Run: python3 -I fig3/offset_tables_v3.py <analysis dir>
"""
import csv
import sys
from pathlib import Path

AN = Path(sys.argv[1])
PP = AN / "fig3/offset_minou/paper_pack"
V2 = AN / "manuscript_offset_v2/tables"
OUT = AN / "manuscript_fig3/tables"
OUT.mkdir(parents=True, exist_ok=True)


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


SN = {r["site_id"]: r["short"] for r in rows(AN / "fig3/site_names.tsv")}
assert len(SN) == 11


def g(x, nd):
    return f"{float(x):.{nd}g}"


def write(name, header, body):
    with open(OUT / f"{name}.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(body)
    print(f"\n### {name}\n| " + " | ".join(header) + " |\n|" + "---|" * len(header))
    for r in body:
        print("| " + " | ".join(str(v) for v in r) + " |")


# Table 1: consensus ranks (order = consensus score, as written by build_paper_pack.R)
cs = rows(PP / "consensus_sites.tsv")
assert "linear_rank" not in cs[0], "consensus_sites.tsv is not the NO_LINEAR build"
write("offset_Table1_consensus_ranks",
      ["site", "site_id", "rank_gf", "rank_geometric", "rank_rda", "consensus_score", "consensus_rank", "category",
       "bio_1_present", "bio_3_present", "d_bio_3", "d_bio_5"],
      [[SN[r["site"]], r["site"], r["rank_gf"], r["rank_geometric"], r["rank_rda"], f"{float(r['consensus_score']):.2f}",
        r["consensus_rank"], r["category"], g(r["bio_1_present"], 3), g(r["bio_3_present"], 3),
        f"{float(r['d_bio_3']):.2f}", f"{float(r['d_bio_5']):.2f}"] for r in cs])

# Table S: random-SNP test, unchanged from v2 (checked against the v3 random_loo.tsv)
t1 = rows(V2 / "TableS1_random_snp_test.tsv")
loo = rows(PP / "random_loo.tsv")
for r in t1:
    adapt = {round(float(x["adaptive"]), 3) for x in loo if x.get("method") == r["method"] and "adaptive" in x}
    assert not adapt or round(float(r["adaptive_tau"]), 3) in adapt, r
write("offset_TableS_random_snp_test", list(t1[0].keys()), [list(r.values()) for r in t1])

# Table S: predictor change vs rank, per method
dt = rows(PP / "delta_rank_tau_sep.tsv")
assert len(dt) == 9
write("offset_TableS_predictor_change_vs_rank", ["method", "predictor", "kendall_tau", "p_exact"],
      [[r["method"], r["predictor"], g(r["kendall_tau"], 3), g(r["p_exact"], 3)] for r in dt])

# Table S: sNMF cluster per site and K, sites by short name
kc = rows(V2 / "TableS3_k_clusters.tsv")
kc.sort(key=lambda r: int(SN[r["site"]]))
write("offset_TableS_k_clusters", ["site", "site_id"] + [k for k in kc[0] if k != "site"],
      [[SN[r["site"]], r["site"]] + [v for k, v in r.items() if k != "site"] for r in kc])
