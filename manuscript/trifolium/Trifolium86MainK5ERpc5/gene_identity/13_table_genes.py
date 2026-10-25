#!/usr/bin/env python3
"""13_table_genes.py -- the genes table of the supplement (slug t_genes): one row per gene of the 15 candidate
regions, from the re-run's gene_identity_all_genes.tsv plus the distance to the nearest associated SNP.

Run with `python3 -I` (stdlib only; gi_common.py is loaded from this script's directory).

Usage:
  python3 -I 13_table_genes.py <data_dir> <gene_identity_all_genes.tsv> <out_dir> \
      [n_genes=353] [n_regions=15] [n_snps=30] [window_bp=130811]

Columns (what the settled Methods says the table reports):
  * region: 1-15 in genomic order (chromosome, then start), with its coordinates from regions_combined.tsv;
  * gene and transcript (BRAKER), position, strand and protein length; strand is the BRAKER GTF strand,
    asserted equal to the pipeline's by 02_prepare_queries.py;
  * the eggNOG description written in the genome annotation (GFF);
  * CD-Search domains: specific hits and superfamily-only hits (accession:short name:range:E-value);
  * the best BLASTP hit in each proteome (Arabidopsis UP000006548, Medicago truncatula UP000002051, soybean
    UP000008827): locus, name, UniProt accession, % identity, query coverage (qcovs), E-value and bit score; the
    PANTHER family:subfamily and the near-tie loci (bit score >= 0.90 x the best) of the best Arabidopsis hit;
    the published MtBURP / MtYUC name of the best Medicago hit;
  * the nearest associated SNP (chr:pos) and the distance to it in bp, 0 when the SNP lies in the gene.
Not in the table: the Swiss-Prot-only hit, the eggNOG seed ortholog and the curated identity (internal screen).
Values are copied as the re-run wrote them, with one exception: a description of "None" (the pipeline's text for a gene
the annotation gives no eggNOG description) is written empty.
The literature keys are added after the gene-literature pass (plan abundant-sniffing-goose, Phase 3).

Checks (the script stops if one fails): 353 genes, unique, all transcripts ".t1"; 15 regions, each holding at least
one associated SNP and every gene of the table in exactly one region; 30 associated SNPs, chromosome names as in
the gene table (1-6, no "Chr_"); every gene's chromosome carries an associated SNP and every gene lies within
window_bp of its nearest one (the regions extend each associated SNP by that distance).

Reads  <data_dir>/lab/{selected_snps.tsv, regions_combined.tsv}, gene_identity_all_genes.tsv
Writes <out_dir>/t_genes.tsv (every column, full precision) and <out_dir>/t_genes.md (a reading copy)
"""
import os
import sys

sys.dont_write_bytecode = True  # keep the scripts directory free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gi_common as gc  # noqa: E402

COLS = ["region", "region_chr", "region_start", "region_end", "gene", "transcript", "chr", "start", "end",
        "strand", "protein_length_aa", "gff_description", "cdsearch_specific_domains",
        "cdsearch_superfamily_domains",
        "at_locus", "at_gene_name", "at_protein_name", "at_uniprot", "at_identity_pct", "at_query_coverage_pct",
        "at_evalue", "at_bitscore", "at_panther", "at_near_tie_loci",
        "mt_locus", "mt_published_name", "mt_protein_name", "mt_uniprot", "mt_identity_pct",
        "mt_query_coverage_pct", "mt_evalue", "mt_bitscore",
        "gm_locus", "gm_protein_name", "gm_uniprot", "gm_identity_pct", "gm_query_coverage_pct", "gm_evalue",
        "gm_bitscore",
        "nearest_associated_snp", "distance_to_nearest_associated_snp_bp"]
# gene_identity_all_genes.tsv column -> t_genes column (renamed only; values are copied as written)
RENAME = {"chr": "chr", "start": "start", "end": "end", "strand": "strand", "protein_length": "protein_length_aa",
          "eggnog_desc": "gff_description", "cdd_specific_hits": "cdsearch_specific_domains",
          "cdd_other_hits": "cdsearch_superfamily_domains",
          "at_agi": "at_locus", "at_gene_name": "at_gene_name", "at_protein_name": "at_protein_name",
          "at_accession": "at_uniprot", "at_pident": "at_identity_pct", "at_qcov": "at_query_coverage_pct",
          "at_evalue": "at_evalue", "at_bitscore": "at_bitscore", "at_panther": "at_panther",
          "at_near_tie_loci": "at_near_tie_loci",
          "mt_locus": "mt_locus", "mt_published_name": "mt_published_name", "mt_protein_name": "mt_protein_name",
          "mt_accession": "mt_uniprot", "mt_pident": "mt_identity_pct", "mt_qcov": "mt_query_coverage_pct",
          "mt_evalue": "mt_evalue", "mt_bitscore": "mt_bitscore",
          "gm_locus": "gm_locus", "gm_protein_name": "gm_protein_name", "gm_accession": "gm_uniprot",
          "gm_pident": "gm_identity_pct", "gm_qcov": "gm_query_coverage_pct", "gm_evalue": "gm_evalue",
          "gm_bitscore": "gm_bitscore"}


def fail(msg):
    sys.exit(f"[13] STOP: {msg}")


def short_domains(text):
    """'pfam03181:BURP:89-301:1.2e-118; ...' -> 'BURP; ...' (short names, order kept)."""
    names = [d.split(":")[1] for d in text.split("; ") if d]
    return ", ".join(dict.fromkeys(names))


def hit_md(locus, name, pident):
    if not locus:
        return "no hit"
    label = f"{locus} {name}" if name and name != locus else locus
    return f"{label} ({float(pident):.0f}%)"


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    data, gi_path, out_dir = sys.argv[1:4]
    n_genes = int(sys.argv[4]) if len(sys.argv) > 4 else 353
    n_regions = int(sys.argv[5]) if len(sys.argv) > 5 else 15
    n_snps = int(sys.argv[6]) if len(sys.argv) > 6 else 30
    window = int(sys.argv[7]) if len(sys.argv) > 7 else 130811

    genes = gc.read_tsv(gi_path)
    snps = gc.read_tsv(os.path.join(data, "lab", "selected_snps.tsv"))
    regions = gc.read_tsv(os.path.join(data, "lab", "regions_combined.tsv"))

    # ---- inputs
    if len(genes) != n_genes:
        fail(f"{len(genes)} genes, expected {n_genes}")
    if len({g["gene"] for g in genes}) != n_genes or not all(g["gene"].endswith(".t1") for g in genes):
        fail("gene ids are not unique '.t1' transcripts")
    if len(snps) != n_snps or len({s["SNPID"] for s in snps}) != n_snps:
        fail(f"{len(snps)} associated SNPs, expected {n_snps}, all distinct")
    for s in snps:
        if s["SNPID"] != f"{s['chr']}:{s['pos']}" or s["chr"].lower().startswith("chr"):
            fail(f"associated SNP id or chromosome name: {s}")
    if len(regions) != n_regions:
        fail(f"{len(regions)} regions, expected {n_regions}")
    gene_chrs = {g["chr"] for g in genes}
    if any(c.lower().startswith("chr") for c in gene_chrs) or not gene_chrs <= {s["chr"] for s in snps}:
        fail(f"gene chromosomes {sorted(gene_chrs)} vs SNP chromosomes {sorted({s['chr'] for s in snps})}")

    # ---- regions 1-15 in genomic order; each holds an associated SNP
    order = sorted(regions, key=lambda r: (int(r["chr"]), int(r["start"])))
    num = {r["region_id"]: i + 1 for i, r in enumerate(order)}
    reg = {r["region_id"]: r for r in regions}
    for r in regions:
        inside = [s for s in snps if s["chr"] == r["chr"] and int(r["start"]) <= int(s["pos"]) <= int(r["end"])]
        if not inside:
            fail(f"region {r['region_id']} holds no associated SNP")
    for name, rid in (("region 1", "1_8649433-8911055"), ("region 12", "5_5557888-5883272"),
                      ("region 14", "5_111174564-111436186")):
        if num.get(rid) != int(name.split()[1]):
            fail(f"{name} is not {rid} (the numbering of the regions curation)")

    # ---- one row per gene
    rows = []
    for g in genes:
        r = reg.get(g["region_id"])
        if r is None or not (r["chr"] == g["chr"] and int(r["start"]) <= int(g["end"]) and int(g["start"]) <= int(r["end"])):
            fail(f"{g['gene']} does not overlap its region {g['region_id']}")
        start, end = int(g["start"]), int(g["end"])
        best = None
        for s in snps:
            if s["chr"] != g["chr"]:
                continue
            p = int(s["pos"])
            d = 0 if start <= p <= end else min(abs(p - start), abs(p - end))
            if best is None or d < best[0] or (d == best[0] and p < best[1]):
                best = (d, p, s["SNPID"])
        if best[0] > window:
            fail(f"{g['gene']} lies {best[0]} bp from its nearest associated SNP, beyond {window} bp")
        row = {"region": num[g["region_id"]], "region_chr": r["chr"], "region_start": r["start"],
               "region_end": r["end"], "gene": g["gene"][:-3], "transcript": g["gene"],
               "nearest_associated_snp": best[2], "distance_to_nearest_associated_snp_bp": best[0]}
        for src, dst in RENAME.items():
            row[dst] = g[src]
        if row["gff_description"] == "None":
            row["gff_description"] = ""
        rows.append(row)
    rows.sort(key=lambda x: (x["region"], int(x["start"]), x["gene"]))
    if [len([x for x in rows if x["region"] == i]) == 0 for i in range(1, n_regions + 1)].count(True):
        fail("a region has no gene")
    gc.write_tsv(os.path.join(out_dir, "t_genes.tsv"), COLS, rows)

    # ---- reading copy
    n_in = sum(1 for x in rows if x["distance_to_nearest_associated_snp_bp"] == 0)
    n_nodesc = sum(1 for x in rows if not x["gff_description"])
    md = ["# t_genes", "",
          f"One row per gene of the {n_regions} candidate regions ({n_genes} genes), regions numbered 1–{n_regions} in "
          "genomic order. Description: the eggNOG description in the genome annotation. Domains: NCBI CD-Search "
          "specific hits (short names; superfamily-only hits are in the TSV). Best BLASTP hit (BLAST+ 2.17.0, "
          "E ≤ 10⁻⁵, UniProt 2026_03) in the reference proteomes of *Arabidopsis thaliana*, *Medicago truncatula* "
          "and soybean (*Glycine max*), ranked by bit score, with its percent identity; the published MtBURP or MtYUC "
          "name is given for Medicago where one exists. PANTHER: family:subfamily of the best Arabidopsis hit. "
          "Distance: to the nearest of the 30 associated SNPs, 0 when the SNP lies in the gene. "
          "Values: `t_genes.tsv`, which also holds the UniProt accessions, query coverage, E-values, bit scores, "
          "near-tie loci and superfamily domains; built by "
          "`scripts/trifolium/gene_identity/13_table_genes.py`.", "",
          f"{n_in} genes contain an associated SNP; {n_nodesc} have no description in the annotation (left empty). "
          "The literature column follows the gene-literature pass.", "",
          "| Region | Gene | Position | Description | Domains | Arabidopsis best hit | PANTHER | Medicago best hit "
          "| Soybean best hit | Nearest associated SNP, kb |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for x in rows:
        mt_name = x["mt_published_name"] or x["mt_locus"]
        md.append("| " + " | ".join([
            str(x["region"]), x["gene"], f"{x['chr']}:{int(x['start']):,}–{int(x['end']):,} ({x['strand']})",
            x["gff_description"].replace("|", "/"), short_domains(x["cdsearch_specific_domains"]),
            hit_md(x["at_locus"], x["at_gene_name"], x["at_identity_pct"]), x["at_panther"],
            hit_md(mt_name, "", x["mt_identity_pct"]), hit_md(x["gm_locus"], "", x["gm_identity_pct"]),
            f"{x['distance_to_nearest_associated_snp_bp'] / 1000:.1f}"]) + " |")
    with open(os.path.join(out_dir, "t_genes.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"[13] t_genes: {len(rows)} genes in {n_regions} regions; {n_in} contain an associated SNP; "
          f"largest distance {max(x['distance_to_nearest_associated_snp_bp'] for x in rows):,} bp")


if __name__ == "__main__":
    main()
