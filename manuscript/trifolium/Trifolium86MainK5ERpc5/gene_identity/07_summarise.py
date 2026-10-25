#!/usr/bin/env python3
"""07_summarise.py -- best hits per gene and database, PANTHER of the best Arabidopsis hit, CD-Search domains.

Run with `python3 -I` (stdlib only; gi_common.py is loaded from this script's directory).

Usage:
  python3 -I 07_summarise.py <data_dir> <out_dir> [top_n=5] [tie_frac=0.90] \
      [region_gene_evidence.tsv evidence_rgt12.tsv]

Rules
  * Per query and subject the best HSP (highest bitscore, then lowest E-value) gives % identity, alignment
    length, E-value and bitscore. Query cover is BLAST's qcovs (all HSPs of that subject); subject cover is the
    merged length of all HSP intervals on the subject divided by its length (not a BLAST field).
  * UniProt entries are collapsed to loci before ranking. Arabidopsis: entries sharing any AGI code
    (ordered locus name) form one locus; Medicago / soybean: the ordered locus name, else the ORF name, else
    the accession. A locus is represented by its highest-scoring entry (ties: reviewed first, then accession).
  * Loci are ranked by bitscore (E-values depend on database size and are not compared across databases).
  * The best Arabidopsis hit is taken from the reference proteome UP000006548 (which contains the Swiss-Prot
    entries); the Swiss-Prot-only best hit is reported separately. "Near-tie" loci are those with a bitscore
    >= tie_frac x the best locus bitscore.
  * PANTHER family:subfamily comes from the UniProt cross-reference of the representative entry.
  * CD-Search: specific hits (accession, short name, query range, E-value); superfamily-only rows listed apart.
  * Medicago hits that are published family members (06: MtBURP / MtYUC) carry that name.
  * When the curated tables are given, the curated identity is carried (region 12: evidence_rgt12.tsv, else
    region_gene_evidence.tsv) and every AGI code it cites is located among the re-run's Arabidopsis loci
    (best / near-tie with its bitscore fraction / lower / not hit). This is a lookup, not a verdict.

Reads  <data_dir>/{work/genes.tsv, blast/*.blastp.tsv, uniprot/*.annot.tsv, cdsearch/cdsearch_hits.tsv,
       work/medicago_family_names.tsv}
Writes <data_dir>/work/blast_top_hits.tsv (top_n loci per gene and database, long format),
       <data_dir>/work/arabidopsis_loci.tsv (every Arabidopsis locus hit, for 08/10),
       <out_dir>/gene_identity_all_genes.tsv
"""
import os
import sys

sys.dont_write_bytecode = True  # keep the scripts directory free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gi_common as gc  # noqa: E402


def fmt(x, nd=1):
    return "" if x is None or x == "" else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    data, out_dir = sys.argv[1], sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    tie_frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.90
    work = os.path.join(data, "work")
    curated = {}
    if len(sys.argv) > 6:
        for r in gc.read_tsv(sys.argv[5]):
            curated[r["gene_id"]] = r["identity"]
        for r in gc.read_tsv(sys.argv[6]):
            if r["identity"]:
                curated[r["transcript"] or r["gene_id"] + ".t1"] = r["identity"]

    genes = gc.read_tsv(os.path.join(work, "genes.tsv"))
    annot = {s: gc.load_annot(os.path.join(data, "uniprot", f"{s}.annot.tsv")) for s in gc.SETS}
    blast = {s: gc.read_blast(os.path.join(data, "blast", f"{s}.blastp.tsv")) for s in gc.SETS}
    fam_name = {}
    for r in gc.read_tsv(os.path.join(work, "medicago_family_names.tsv")):
        for a in r["uniprot_accessions"].split(";"):
            if a:
                fam_name[a] = r["published_name"]

    cds_path = os.path.join(data, "cdsearch", "cdsearch_hits.tsv")
    cds_spec, cds_other, cds_ran = {}, {}, os.path.exists(cds_path)
    if cds_ran:
        for r in gc.read_tsv(cds_path):
            label = (f"{r['accession']}:{r['short_name']}:{r['from']}-{r['to']}:{float(r['e_value']):.1e}"
                     + (f":incomplete_{r['incomplete'].strip()}" if r["incomplete"].strip() not in ("-", "") else ""))
            (cds_spec if r["hit_type"] == "specific" else cds_other).setdefault(r["query_id"], []).append(label)

    long_rows, at_rows, out_rows = [], [], []
    for g in genes:
        q = g["transcript_id"]
        best = {}
        for s in gc.SETS:
            loci = gc.collapse_loci(s, blast[s].get(q, []), annot[s])
            best[s] = loci
            top_bits = loci[0]["rep"]["bitscore"] if loci else None
            for rank, L in enumerate(loci[:top_n], 1):
                h, row = L["rep"], L["row"]
                long_rows.append({
                    "transcript_id": q, "database": s, "rank": rank, "locus": L["locus"],
                    "accession": h["sseqid"], "reviewed": row.get("Reviewed", ""),
                    "gene_symbol": gc.gene_symbol(row), "protein_name": gc.first_protein_name(row.get("Protein names")),
                    "published_family_name": fam_name.get(h["sseqid"], "") if s == "mt_proteome" else "",
                    "panther": gc.panther_label(row.get("PANTHER")), "pident": fmt(h["pident"]),
                    "align_len": h["length"], "qcovs": fmt(h["qcovs"], 0), "scov": fmt(100 * h["scov"], 1),
                    "evalue": f"{h['evalue']:.1e}", "bitscore": fmt(h["bitscore"]), "qlen": h["qlen"],
                    "slen": h["slen"], "frac_of_best_bitscore": fmt(h["bitscore"] / top_bits, 3)})
            if s == "at_proteome":
                for L in loci:
                    h, row = L["rep"], L["row"]
                    at_rows.append({
                        "transcript_id": q, "locus": L["locus"], "agis": ";".join(L["agis"]),
                        "accession": h["sseqid"], "all_accessions": ";".join(L["accessions"]),
                        "gene_symbol": gc.gene_symbol(row),
                        "protein_name": gc.first_protein_name(row.get("Protein names")),
                        "panther": gc.panther_label(row.get("PANTHER")), "pident": fmt(h["pident"]),
                        "qcovs": fmt(h["qcovs"], 0), "scov": fmt(100 * h["scov"], 1), "evalue": f"{h['evalue']:.1e}",
                        "bitscore": fmt(h["bitscore"]), "slen": h["slen"],
                        "frac_of_best_bitscore": fmt(h["bitscore"] / top_bits, 3)})

        def bestcols(prefix, s, agi=False, famname=False):
            loci = best[s]
            if not loci:
                return {f"{prefix}_accession": "no hit (E <= cutoff)"}
            L = loci[0]
            h, row = L["rep"], L["row"]
            d = {f"{prefix}_accession": h["sseqid"],
                 f"{prefix}_{'agi' if agi else 'locus'}": L["locus"],
                 f"{prefix}_gene_name": gc.gene_symbol(row),
                 f"{prefix}_protein_name": gc.first_protein_name(row.get("Protein names")),
                 f"{prefix}_pident": fmt(h["pident"]), f"{prefix}_qcov": fmt(h["qcovs"], 0),
                 f"{prefix}_evalue": f"{h['evalue']:.1e}", f"{prefix}_bitscore": fmt(h["bitscore"])}
            if famname:
                d[f"{prefix}_published_name"] = fam_name.get(h["sseqid"], "")
            return d

        rec = {"region_id": g["region_id"], "gene": q, "chr": g["chr"], "start": g["start"], "end": g["end"],
               "strand": g["strand"], "protein_length": g["protein_length"], "eggnog_desc": g["eggnog_desc"],
               "eggnog_seed_ortholog": g["em_target"]}
        rec.update(bestcols("at", "at_proteome", agi=True))
        if best["at_proteome"]:
            L = best["at_proteome"][0]
            rec["at_panther"] = gc.panther_label(L["row"].get("PANTHER"))
            rec["at_scov"] = fmt(100 * L["rep"]["scov"], 1)
            top = L["rep"]["bitscore"]
            ties = [M for M in best["at_proteome"] if M["rep"]["bitscore"] >= tie_frac * top]
            rec["at_near_tie_loci"] = "; ".join(
                f"{M['locus']} {gc.gene_symbol(M['row']) or gc.first_protein_name(M['row'].get('Protein names'))[:40]}"
                f" {M['rep']['bitscore']:.0f}" for M in ties[1:])
        rec.update(bestcols("sp", "at_swissprot", agi=True))
        rec.update(bestcols("mt", "mt_proteome", famname=True))
        rec.update(bestcols("gm", "gm_proteome"))
        cur = curated.get(q, "")
        rec["curated_identity"] = cur
        looked = []
        for a in gc.agis(cur):
            L = next((M for M in best["at_proteome"] if a in M["agis"]), None)
            if L is None:
                looked.append(f"{a} not hit")
            else:
                fr = L["rep"]["bitscore"] / best["at_proteome"][0]["rep"]["bitscore"]
                looked.append(f"{a} " + ("best" if L is best["at_proteome"][0]
                                          else f"{'near-tie' if fr >= tie_frac else 'lower'} {fr:.3f}"))
        rec["curated_agis_in_rerun"] = "; ".join(looked)
        rec["cdd_specific_hits"] = "; ".join(cds_spec.get(q, [])) if cds_ran else "CD-Search not run"
        rec["cdd_other_hits"] = "; ".join(cds_other.get(q, [])) if cds_ran else "CD-Search not run"
        out_rows.append(rec)

    long_cols = ["transcript_id", "database", "rank", "locus", "accession", "reviewed", "gene_symbol", "protein_name",
                 "published_family_name", "panther", "pident", "align_len", "qcovs", "scov", "evalue", "bitscore",
                 "qlen", "slen", "frac_of_best_bitscore"]
    gc.write_tsv(os.path.join(work, "blast_top_hits.tsv"), long_cols, long_rows)
    at_cols = ["transcript_id", "locus", "agis", "accession", "all_accessions", "gene_symbol", "protein_name", "panther",
               "pident", "qcovs", "scov", "evalue", "bitscore", "slen", "frac_of_best_bitscore"]
    gc.write_tsv(os.path.join(work, "arabidopsis_loci.tsv"), at_cols, at_rows)
    cols = ["region_id", "gene", "chr", "start", "end", "strand", "protein_length", "eggnog_desc",
            "eggnog_seed_ortholog",
            "at_accession", "at_agi", "at_gene_name", "at_protein_name", "at_panther", "at_pident", "at_qcov",
            "at_scov", "at_evalue", "at_bitscore", "at_near_tie_loci",
            "sp_accession", "sp_agi", "sp_gene_name", "sp_protein_name", "sp_pident", "sp_qcov", "sp_evalue",
            "sp_bitscore",
            "mt_accession", "mt_locus", "mt_published_name", "mt_protein_name", "mt_pident", "mt_qcov", "mt_evalue",
            "mt_bitscore",
            "gm_accession", "gm_locus", "gm_protein_name", "gm_pident", "gm_qcov", "gm_evalue", "gm_bitscore",
            "cdd_specific_hits", "cdd_other_hits", "curated_identity", "curated_agis_in_rerun"]
    os.makedirs(out_dir, exist_ok=True)
    gc.write_tsv(os.path.join(out_dir, "gene_identity_all_genes.tsv"), cols, out_rows)
    n_at = sum(1 for r in out_rows if r.get("at_agi"))
    print(f"[07] genes {len(out_rows)}; with an Arabidopsis hit {n_at}; "
          f"Swiss-Prot {sum(1 for r in out_rows if r.get('sp_agi') or r.get('sp_gene_name'))}; "
          f"Medicago {sum(1 for r in out_rows if r.get('mt_locus'))}; soybean {sum(1 for r in out_rows if r.get('gm_locus'))}; "
          f"CD-Search specific {sum(1 for r in out_rows if r['cdd_specific_hits'] and r['cdd_specific_hits'] != 'CD-Search not run')}")


if __name__ == "__main__":
    main()
