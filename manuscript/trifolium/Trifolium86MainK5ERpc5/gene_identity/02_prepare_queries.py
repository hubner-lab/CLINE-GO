#!/usr/bin/env python3
"""02_prepare_queries.py -- build the 353-gene query set and check it against the gene models.

Run with `python3 -I` (isolated mode; stdlib only).

Usage:
  python3 -I 02_prepare_queries.py <data_dir> <region_gene_evidence.tsv> [expected_n_genes=353]

Reads  <data_dir>/lab/{genes_per_region.tsv, regions_combined.tsv, Trifolium_genes.gtf, braker_all_proteins.faa}
Writes <data_dir>/work/genes.tsv        one row per gene: region, position rank, model checks, eggNOG fields
       <data_dir>/work/queries.faa      the Trifolium proteins (terminal stop already removed by gffread)
       <data_dir>/work/model_checks.txt summary of every assertion and flag

Checks (hard failures stop the script; flags are reported, never used to drop a gene):
  * the pipeline table holds exactly the expected number of distinct transcripts;
  * each gene falls in exactly one of the combined regions (coordinate overlap, chr N vs Chr_N);
  * that region equals the curated region_id in region_gene_evidence.tsv for every gene;
  * the BRAKER transcript span equals gene_start/gene_end of the pipeline table;
  * the protein length equals the CDS length / 3 (minus the terminal stop when a stop codon is annotated).
Flags: no start codon feature, no stop codon feature, protein not starting with M, in-frame stop ('*').
"""
import csv
import os
import re
import sys
from collections import defaultdict

CHR_PREFIX = "Chr_"


def die(msg):
    sys.exit("[02] ERROR: " + msg)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    if len(sys.argv) < 3:
        die(__doc__)
    data, evidence_path = sys.argv[1], sys.argv[2]
    expected_n = int(sys.argv[3]) if len(sys.argv) > 3 else 353
    lab, work = os.path.join(data, "lab"), os.path.join(data, "work")
    os.makedirs(work, exist_ok=True)
    log = []

    # 1. distinct transcripts of the pipeline table, with their (identical) per-row fields
    gpr = read_tsv(os.path.join(lab, "genes_per_region.tsv"))
    genes = {}
    for r in gpr:
        tid = r["gene_id"]
        rec = {k: r[k] for k in ("chr", "gene_start", "gene_end", "strand", "em_target", "em_desc",
                                 "em_PFAMs", "em_tcov")}
        if tid in genes and genes[tid] != rec:
            die(f"{tid}: per-trait rows disagree on gene fields")
        genes[tid] = rec
    log.append(f"pipeline rows {len(gpr)}, distinct transcripts {len(genes)} (expected {expected_n})")
    if len(genes) != expected_n:
        die(f"expected {expected_n} distinct transcripts, found {len(genes)}")

    # 2. assign each gene to one combined region by coordinate overlap
    regions = read_tsv(os.path.join(lab, "regions_combined.tsv"))
    for tid, g in genes.items():
        hits = [rg["region_id"] for rg in regions
                if rg["chr"] == g["chr"] and int(g["gene_start"]) <= int(rg["end"])
                and int(g["gene_end"]) >= int(rg["start"])]
        if len(hits) != 1:
            die(f"{tid}: overlaps {len(hits)} combined regions {hits}")
        g["region_id"] = hits[0]
    log.append(f"combined regions {len(regions)}; every gene overlaps exactly one")

    # 3. compare with the curated evidence table (gene ids there carry the .t1 suffix too)
    ev = {r["gene_id"]: r for r in read_tsv(evidence_path)}
    if set(ev) != set(genes):
        die(f"gene sets differ: only pipeline {sorted(set(genes) - set(ev))[:5]}, "
            f"only curated {sorted(set(ev) - set(genes))[:5]}")
    bad = [t for t in genes if ev[t]["region_id"] != genes[t]["region_id"]]
    if bad:
        die(f"region assignment differs from the curated table for {bad[:5]}")
    log.append("gene set and region assignment equal the curated region_gene_evidence.tsv (353/353)")

    # 4. BRAKER models: span, CDS length, start/stop codon features
    model = defaultdict(lambda: {"chr": None, "min": None, "max": None, "strand": None,
                                 "cds_len": 0, "n_cds": 0, "start_codon": False, "stop_codon": False})
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    with open(os.path.join(lab, "Trifolium_genes.gtf")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            m = tid_re.search(f[8])
            tid = m.group(1) if m else (f[8].strip() if f[2] == "transcript" else None)
            if tid not in genes:
                continue
            s, e = int(f[3]), int(f[4])
            mo = model[tid]
            mo["chr"], mo["strand"] = f[0], f[6]
            mo["min"] = s if mo["min"] is None else min(mo["min"], s)
            mo["max"] = e if mo["max"] is None else max(mo["max"], e)
            if f[2] == "CDS":
                mo["cds_len"] += e - s + 1
                mo["n_cds"] += 1
            elif f[2] == "start_codon":
                mo["start_codon"] = True
            elif f[2] == "stop_codon":
                mo["stop_codon"] = True
    missing = [t for t in genes if t not in model]
    if missing:
        die(f"transcripts absent from the BRAKER GTF: {missing[:5]}")
    span_bad = []
    for tid, g in genes.items():
        mo = model[tid]
        if (mo["chr"] != CHR_PREFIX + g["chr"] or mo["min"] != int(g["gene_start"])
                or mo["max"] != int(g["gene_end"]) or mo["strand"] != g["strand"]):
            span_bad.append(tid)
    if span_bad:
        die(f"BRAKER span differs from the pipeline table for {span_bad[:5]}")
    log.append("BRAKER transcript span and strand equal gene_start/gene_end/strand for all genes")

    # 5. proteins: subset, length check, flags
    prot_all = read_fasta(os.path.join(lab, "braker_all_proteins.faa"))
    flags_count = defaultdict(int)
    with open(os.path.join(work, "queries.faa"), "w") as out:
        for tid in sorted(genes, key=lambda t: (genes[t]["region_id"], int(genes[t]["gene_start"]))):
            if tid not in prot_all:
                die(f"{tid}: no protein from gffread")
            p = prot_all[tid]
            mo = model[tid]
            expect = mo["cds_len"] // 3 - (1 if mo["stop_codon"] else 0)
            if len(p) != expect:
                die(f"{tid}: protein length {len(p)} != CDS-derived {expect}")
            flags = []
            if not mo["start_codon"]:
                flags.append("no_start_codon_feature")
            if not mo["stop_codon"]:
                flags.append("no_stop_codon_feature")
            if not p.startswith("M"):
                flags.append("not_M_start")
            if "*" in p or "." in p:
                flags.append("inframe_stop")
            for fl in flags:
                flags_count[fl] += 1
            genes[tid].update(protein_length=len(p), n_cds=mo["n_cds"], cds_length=mo["cds_len"],
                              model_flags=",".join(flags))
            out.write(f">{tid}\n")
            for i in range(0, len(p), 60):
                out.write(p[i:i + 60] + "\n")
    log.append("protein length equals the CDS length / 3 (minus the stop) for all genes")
    log.append("model flags: " + (", ".join(f"{k} {v}" for k, v in sorted(flags_count.items())) or "none"))

    # 6. identical spans (two BRAKER genes over one locus)
    by_span = defaultdict(list)
    for tid, g in genes.items():
        by_span[(g["chr"], g["gene_start"], g["gene_end"], g["strand"])].append(tid)
    for span, tids in sorted(by_span.items()):
        if len(tids) > 1:
            desc = "; ".join(f"{t} {model[t]['n_cds']} CDS {genes[t]['protein_length']} aa" for t in sorted(tids))
            log.append(f"identical span {span[0]}:{span[1]}-{span[2]}{span[3]}: {desc}")

    # 7. position rank within region and the gene table
    by_region = defaultdict(list)
    for tid, g in genes.items():
        by_region[g["region_id"]].append(tid)
    cols = ["region_id", "rank_in_region", "transcript_id", "gene", "chr", "start", "end", "strand",
            "n_cds", "cds_length", "protein_length", "model_flags", "em_target", "em_tcov", "em_PFAMs",
            "eggnog_desc"]
    with open(os.path.join(work, "genes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for rid in sorted(by_region, key=lambda r: (r.split("_")[0], int(r.split("_")[1].split("-")[0]))):
            tids = sorted(by_region[rid], key=lambda t: (int(genes[t]["gene_start"]), int(genes[t]["gene_end"]), t))
            for rank, tid in enumerate(tids, 1):
                g = genes[tid]
                w.writerow([rid, rank, tid, tid.rsplit(".t", 1)[0], g["chr"], g["gene_start"], g["gene_end"],
                            g["strand"], g["n_cds"], g["cds_length"], g["protein_length"], g["model_flags"],
                            g["em_target"], g["em_tcov"], g["em_PFAMs"], ev[tid]["desc"]])
    log.append(f"regions: " + ", ".join(f"{r} {len(by_region[r])}" for r in sorted(by_region)))

    with open(os.path.join(work, "model_checks.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n".join("[02] " + x for x in log))


if __name__ == "__main__":
    main()
