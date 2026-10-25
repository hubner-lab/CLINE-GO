#!/usr/bin/env python3
"""08_microsynteny.py -- Arabidopsis microsynteny of the genes in regions 1 and 12.

Run with `python3 -I` (stdlib only; gi_common.py is loaded from this script's directory).

Usage:
  python3 -I 08_microsynteny.py <data_dir> <regions> [frac=0.90] [agi_window=100] [min_support=2]
  (<regions>: space- or comma-separated region ids, e.g. "1_8649433-8911055 5_5557888-5883272")

Rule (fixed before the run; parameters in config.sh):
  1. Genes of a region are ordered by start position.
  2. Candidates of a gene: Arabidopsis loci (UP000006548 hits collapsed to AGI groups, 07) with a bitscore
     >= frac x the gene's best locus bitscore. Only nuclear AGI codes (AT1G-AT5G) count; plastid/mitochondrial
     codes (ATCG/ATMG) are dropped.
  3. Two candidates are linked when they belong to two different genes of the same region, are different
     Arabidopsis loci (no shared AGI), lie on the same Arabidopsis chromosome and their AGI numbers differ by
     <= agi_window (AGI codes step by about 10 per gene, so 100 is about 10 genes).
  4. A candidate's support is the number of OTHER region genes with at least one candidate linked to it.
     It is "syntenic" when support >= min_support. Requiring a different locus stops a tandem array whose
     copies all hit one Arabidopsis gene from supporting itself. The number of distinct Arabidopsis loci behind
     that support is reported as well (support_loci_n): Trifolium tandem copies that hit one Arabidopsis locus
     count as several supporting genes but as one locus.
  5. Tie-break of a gene: the syntenic candidate with the largest support; then the higher bitscore.
     No syntenic candidate -> no synteny call.
  6. Collinear runs (reported for reading): connected components of linked candidates (one node per
     gene x AGI); a run lists its genes in Trifolium order and states whether the AGI numbers are monotonic.
Writes <data_dir>/work/microsynteny_candidates.tsv and <data_dir>/work/microsynteny_runs.tsv
"""
import os
import sys
from collections import defaultdict

sys.dont_write_bytecode = True  # keep the scripts directory free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gi_common as gc  # noqa: E402


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    data = sys.argv[1]
    regions = [r for r in sys.argv[2].replace(",", " ").split() if r]
    frac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.90
    window = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    min_support = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    work = os.path.join(data, "work")
    genes = gc.read_tsv(os.path.join(work, "genes.tsv"))
    loci = defaultdict(list)
    for r in gc.read_tsv(os.path.join(work, "arabidopsis_loci.tsv")):
        loci[r["transcript_id"]].append(r)

    cand_rows, run_rows = [], []
    for region in regions:
        rg = sorted([g for g in genes if g["region_id"] == region], key=lambda g: int(g["rank_in_region"]))
        if not rg:
            sys.exit(f"[08] ERROR: no genes for region {region}")
        rank = {g["transcript_id"]: int(g["rank_in_region"]) for g in rg}
        cands = {}
        for g in rg:
            q = g["transcript_id"]
            ls = sorted(loci.get(q, []), key=lambda r: -float(r["bitscore"]))
            if not ls:
                cands[q] = []
                continue
            top = float(ls[0]["bitscore"])
            keep = []
            for r in ls:
                if float(r["bitscore"]) < frac * top:
                    break
                nuc = [a for a in r["agis"].split(";") if a and gc.is_nuclear(a)]
                if nuc:
                    keep.append(dict(r, nuc=nuc, is_best=(r is ls[0])))
            cands[q] = keep

        def linked(c, d):
            if set(c["nuc"]) & set(d["nuc"]):
                return False
            for a in c["nuc"]:
                ca, na = gc.agi_pos(a)
                for b in d["nuc"]:
                    cb, nb = gc.agi_pos(b)
                    if ca == cb and abs(na - nb) <= window:
                        return True
            return False

        # support per candidate
        for g in rg:
            q = g["transcript_id"]
            best_sel = None
            recs = []
            for c in cands[q]:
                sup, sup_loci = [], set()
                for h in rg:
                    p = h["transcript_id"]
                    if p == q:
                        continue
                    hits = [d for d in cands[p] if linked(c, d)]
                    if hits:
                        sup.append(f"{p.split('.')[0]}:{'/'.join(sorted({a for d in hits for a in d['nuc']}))}")
                        sup_loci.update(d["locus"] for d in hits)
                recs.append((c, sup, sup_loci))
            syn = [(c, s) for c, s, _ in recs if len(s) >= min_support]
            if syn:
                best_sel = sorted(syn, key=lambda cs: (-len(cs[1]), -float(cs[0]["bitscore"])))[0][0]
            for c, sup, sup_loci in recs:
                cand_rows.append({
                    "region_id": region, "rank_in_region": rank[q], "transcript_id": q, "candidate_locus": c["locus"],
                    "gene_symbol": c["gene_symbol"], "protein_name": c["protein_name"][:60],
                    "bitscore": c["bitscore"], "frac_of_best_bitscore": c["frac_of_best_bitscore"],
                    "best_hit": "yes" if c["is_best"] else "no", "n_candidates": len(cands[q]),
                    "support_n": len(sup), "support_loci_n": len(sup_loci), "supporting_genes": "; ".join(sup),
                    "syntenic": "yes" if len(sup) >= min_support else "no",
                    "synteny_choice": "yes" if c is best_sel else "no"})
            if not cands[q]:
                cand_rows.append({"region_id": region, "rank_in_region": rank[q], "transcript_id": q,
                                  "candidate_locus": "no Arabidopsis hit", "n_candidates": 0, "support_n": 0,
                                  "support_loci_n": 0,
                                  "syntenic": "no", "synteny_choice": "no"})

        # collinear runs: components over (gene, AGI) nodes
        nodes = [(g["transcript_id"], a, c) for g in rg for c in cands[g["transcript_id"]] for a in c["nuc"]]
        parent = list(range(len(nodes)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                (qi, ai, ci), (qj, aj, cj) = nodes[i], nodes[j]
                if qi != qj and linked(ci, cj):
                    ca, na = gc.agi_pos(ai)
                    cb, nb = gc.agi_pos(aj)
                    if ca == cb and abs(na - nb) <= window and not set(ci["nuc"]) & set(cj["nuc"]):
                        parent[find(i)] = find(j)
        comp = defaultdict(list)
        for i in range(len(nodes)):
            comp[find(i)].append(nodes[i])
        k = 0
        for members in sorted(comp.values(), key=lambda m: min(rank[q] for q, _, _ in m)):
            gs = sorted({q for q, _, _ in members}, key=lambda q: rank[q])
            if len(gs) < 2:
                continue
            k += 1
            order = sorted(members, key=lambda m: (rank[m[0]], gc.agi_pos(m[1])[1]))
            nums = [gc.agi_pos(a)[1] for _, a, _ in order]
            mono = all(x <= y for x, y in zip(nums, nums[1:])) or all(x >= y for x, y in zip(nums, nums[1:]))
            chrom = gc.agi_pos(order[0][1])[0]
            run_rows.append({
                "region_id": region, "run": f"{region}#{k}", "arabidopsis_chr": chrom,
                "agi_range": f"{min(a for _, a, _ in members)}-{max(a for _, a, _ in members)}",
                "n_genes": len(gs), "genes_in_order": " ".join(q.split(".")[0] for q in gs),
                "anchors_in_order": " ".join(f"{q.split('.')[0]}={a}" for q, a, _ in order),
                "agi_order_monotonic": "yes" if mono else "no"})

    cand_cols = ["region_id", "rank_in_region", "transcript_id", "candidate_locus", "gene_symbol", "protein_name",
                 "bitscore", "frac_of_best_bitscore", "best_hit", "n_candidates", "support_n", "support_loci_n",
                 "supporting_genes", "syntenic", "synteny_choice"]
    gc.write_tsv(os.path.join(work, "microsynteny_candidates.tsv"), cand_cols, cand_rows)
    run_cols = ["region_id", "run", "arabidopsis_chr", "agi_range", "n_genes", "genes_in_order", "anchors_in_order",
                "agi_order_monotonic"]
    gc.write_tsv(os.path.join(work, "microsynteny_runs.tsv"), run_cols, run_rows)
    print(f"[08] regions {len(regions)}; candidate rows {len(cand_rows)}; runs with >= 2 genes {len(run_rows)} "
          f"(frac {frac}, AGI window {window}, min support {min_support})")


if __name__ == "__main__":
    main()
