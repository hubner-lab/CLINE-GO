#!/usr/bin/env python3
"""tables/literature_raw.tsv for the 15 PC5 regions (plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7).

  carried  every row of the PC3 run's literature_raw.tsv (ANALYSIS/Trifolium86MainK5/regionsER/tables/) whose gene is a
           gene of a PC5 region; region_id rewritten to the PC5 region of that gene (two PC3 regions changed id).
           Rows of genes in the two dropped chr-6 regions are left out and counted.
  new      the Phase 7 literature agents' rows, gene ids without ".t1", region_id added:
             logs/agent_batch_v2.tsv        the new genes sent to literature by write_curation.py (lit = TRUE);
             logs/agent_batch_biotic_A.tsv  defence-class genes (NB-LRR arrays of 2_176.1M), and
             logs/agent_batch_biotic_B.tsv  (NB-LRRs of 5_111M, GH18 chitinase array of 5_5.6M, VQ g26559 of 4_8.1M):
             user request 2026-10-16 to test whether disease-resistance genes have abiotic roles instead of
             rejecting them as biotic by rule. Their curation calls are not changed here.
Then lit_check.py re-verifies every row (carried ones included) against NCBI.
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/carry_literature.py
"""
import csv
import os

B = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.join(B, "../../Trifolium86MainK5/regionsER/tables/literature_raw.tsv")


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


reg = {r["gene_id"].replace(".t1", ""): r["region_id"] for r in rows(os.path.join(B, "tables/curation.tsv"))}
lit_new = {r["gene_id"].replace(".t1", "") for r in rows(os.path.join(B, "tables/curation.tsv"))
           if r["lit"] == "TRUE" and r["source"] != "carried"}
v1 = rows(V1)
cols = list(v1[0])
out, n_drop = [], 0
for r in v1:
    if r["gene_id"] not in reg:
        n_drop += 1
        continue
    out.append({**r, "region_id": reg[r["gene_id"]]})
n_car = len(out)
BIOTIC = {"g15018", "g15021", "g15022", "g15023", "g15024", "g15036", "g15041", "g15043", "g15046",
          "g38791", "g38792", "g32849", "g32850", "g32851", "g32852", "g26559"}
ag = []
for f in ("agent_batch_v2.tsv", "agent_batch_biotic_A.tsv", "agent_batch_biotic_B.tsv"):
    ag += rows(os.path.join(B, "logs", f))
assert {r["gene_id"] for r in ag} <= lit_new | BIOTIC, {r["gene_id"] for r in ag} - lit_new - BIOTIC
for r in ag:
    assert r["role"] in ("key", "support", "conflict", "rejected"), r
    out.append({"region_id": reg[r["gene_id"]], **{k: r[k] for k in cols if k != "region_id"}})
with open(os.path.join(B, "tables/literature_raw.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
print(f"INFO: {len(v1)} PC3 rows: carried {n_car}, left out with the dropped regions {n_drop}; "
      f"agent rows {len(ag)} for {len({r['gene_id'] for r in ag})} genes ({len(lit_new)} new lit genes + {len(BIOTIC)} defence-class)")
