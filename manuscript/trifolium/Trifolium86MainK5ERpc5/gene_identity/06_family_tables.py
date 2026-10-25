#!/usr/bin/env python3
"""06_family_tables.py -- published Medicago family names (MtBURP, MtYUC) mapped to UniProt accessions.

Run with `python3 -I` (stdlib only).

Usage:
  python3 -I 06_family_tables.py <data_dir> [li2016_burp_figure_classes.tsv]

Sources (open-access full text, fetched as JATS XML from PubMed Central through NCBI E-utilities efetch):
  BURP  Li Y et al. 2016, Front Plant Sci 7:485, doi:10.3389/fpls.2016.00485, PMID 27148311, PMC4829796
        Table 1: MtBURP1-39 -> Mt4.0 locus (MedtrNgNNNNNN.v) and protein length.
  YUC   Shao A et al. 2023, Front Plant Sci 14:1268027, doi:10.3389/fpls.2023.1268027, PMID 37701802, PMC10494245
        Table 1: MtYUC1-15 -> GenBank protein id (AES... / KEH...) and ORF length (aa).
Mapping to UniProt (UP000002051, the BLAST database of step 04):
  MtBURP: Medtr locus -> UniProt entries whose gene names carry MTR_<locus> or Medtr<locus> (the UniProt
          ordered-locus names of the Mt4.0 annotation use the MTR_ prefix for the same locus ids).
  MtYUC:  GenBank protein id -> UniProtKB through the UniProt ID-mapping service
          (from EMBL-GenBank-DDBJ_CDS), kept only if the entry is in the downloaded Medicago proteome.
  Published and UniProt protein lengths are compared and reported (a mismatch means a different gene-model
  version of the same locus).
Class / clade membership from the article text: each statement is stored below verbatim and the script fails
if it is not found in the downloaded XML. Members not named in any statement get "not stated in text".
Li 2016 draws the full subfamily membership only in figures (Fig. 2 tree arcs, Fig. 1 label colours). Those were
read by eye from the publisher PDF into li2016_burp_figure_classes.tsv (a hand-made input, like the claim spec);
when given, the figure class is added in its own column and must agree with every class stated in the text
(the script stops otherwise), which checks the reading against the 23 members the text does name.

Also maps the eggNOG seed orthologs (em_target) of the genes in regions 1 and 12 that are Medicago GenBank
proteins (3880.<id>) to UniProt, to show which published family member each curated "best hit" was.

Writes <data_dir>/pmc/PMC*.xml (downloads), <data_dir>/work/medicago_family_names.tsv,
       <data_dir>/work/eggnog_seed_medicago.tsv, <data_dir>/work/family_sources.tsv
"""
import csv
import datetime
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
IDMAP = "https://rest.uniprot.org/idmapping"
TOOL = "trifolium-gene-identity"

PAPERS = {
    "BURP": {"pmcid": "4829796", "pmid": "27148311", "doi": "10.3389/fpls.2016.00485",
             "citation": "Li Y, Chen X, Chen Z, Cai R, Zhang H, Xiang Y (2016) Front Plant Sci 7:485"},
    "YUC": {"pmcid": "10494245", "pmid": "37701802", "doi": "10.3389/fpls.2023.1268027",
            "citation": "Shao A, Fan S, Xu X, Wang W, Fu J (2023) Front Plant Sci 14:1268027"},
}

# Class / clade statements: (family, class label, members, verbatim sentence from the article text).
STATEMENTS = [
    ("BURP", "RD22-like", ["MtBURP5", "MtBURP10", "MtBURP11", "MtBURP19"],
     "Among the 10 up-regulated genes, four genes (MtBURP5, -10, -11, and -19) belong to RD22-like"),
    ("BURP", "RD22-like", ["MtBURP11"], "MtBURP11 from RD22-like subfamily"),
    ("BURP", "RD22-like", ["MtBURP3", "MtBURP7"], "The highest expression levels of MtBURP3 and -7 of RD22-like subfamily"),
    ("BURP", "PG1beta-like", ["MtBURP25", "MtBURP26"], "two genes (MtBURP25 and -26) to PG1β"),
    ("BURP", "PG1beta-like", ["MtBURP33", "MtBURP28"], "only two genes (MtBURP33 and -28) both included in PG1β"),
    ("BURP", "BNM2-like", ["MtBURP16"], "one genes (MtBURP16) to BNM2-like"),
    ("BURP", "USP-like", ["MtBURP2", "MtBURP15", "MtBURP39"], "others (MtBURP2, -15 and -39) to USP-like"),
    ("BURP", "USP-like", ["MtBURP15"], "MtBURP15 from USP-like subfamily"),
    ("YUC", "clade I-1", ["MtYUC3", "MtYUC4", "MtYUC12"], "MtYUC3, 4, 12 in clade I-1"),
    ("YUC", "clade I-2 (relatively close to AtYUC6)", ["MtYUC1", "MtYUC8", "MtYUC11"],
     "MsYUC7/9 and MtYUC1/8/11, belonging to clade I-2, were relatively closely related to AtYUC6"),
    ("YUC", "clade I-4 (close to AtYUC1/AtYUC4)", ["MtYUC2"],
     "MtYUC2 showed a close relationship with AtYUC1 and AtYUC4, which were clustered into clade I-4"),
    ("YUC", "clade II-4 (close to AtYUC10)", ["MtYUC10", "MtYUC15"],
     "MtYUC10, MtYUC15, and MsYUC8 were closely related to AtYUC10, which belongs to clade II-4"),
]
# Table 2 of Li 2016 groups paralogous pairs under a class label in its first column ("RD22", "BNM2", ...);
# those pairs are read from the table itself (see parse_li_table2).

SYN_REGIONS = ("1_8649433-8911055", "5_5557888-5883272")


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def http(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"User-Agent": TOOL + "/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read().decode("utf-8", "replace")


def cell_text(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def table_rows(root, label):
    for tw in root.iter("table-wrap"):
        lab = tw.find("label")
        if lab is not None and cell_text(lab) == label:
            return [[cell_text(c) for c in tr] for tr in tw.iter("tr")]
    sys.exit(f"[06] ERROR: {label} not found")


def article_text(root):
    return re.sub(r"\s+", " ", " ".join("".join(p.itertext()) for p in root.iter("p")))


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def uniprot_idmap(ids, from_db="EMBL-GenBank-DDBJ_CDS"):
    body = urllib.parse.urlencode({"from": from_db, "to": "UniProtKB", "ids": ",".join(ids)}).encode()
    job = json.loads(http(f"{IDMAP}/run", body))["jobId"]
    for _ in range(60):
        time.sleep(3)
        st = json.loads(http(f"{IDMAP}/status/{job}"))
        if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
            break
        if st.get("jobStatus") not in ("RUNNING", "NEW"):
            sys.exit(f"[06] ERROR: ID mapping job {job}: {st}")
    tsv = http(f"{IDMAP}/uniprotkb/results/stream/{job}?format=tsv&fields=accession,gene_names,gene_oln,length")
    out = {}
    lines = tsv.splitlines()
    for ln in lines[1:]:
        f = ln.split("\t")
        out.setdefault(f[0], []).append({"accession": f[1], "gene_names": f[2], "oln": f[3], "length": f[4]})
    return out, job


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    data = sys.argv[1]
    fig = {}
    if len(sys.argv) > 2:
        for r in read_tsv(sys.argv[2]):
            fig[r["published_name"]] = f"{r['subfamily_fig2']} (Fig. 2 arc; Fig. 1 colour {r['fig1_colour']})"
    pmc, work = os.path.join(data, "pmc"), os.path.join(data, "work")   # pmc/: downloaded text, untrusted
    os.makedirs(pmc, exist_ok=True)
    os.makedirs(work, exist_ok=True)

    roots, sources = {}, []
    for fam, p in PAPERS.items():
        path = os.path.join(pmc, f"PMC{p['pmcid']}.xml")
        url = f"{EFETCH}?db=pmc&id={p['pmcid']}&rettype=xml&tool={TOOL}"
        xml = http(url)
        with open(path, "w") as fh:
            fh.write(xml)
        roots[fam] = ET.fromstring(xml)
        sources.append([fam, p["citation"], p["doi"], p["pmid"], "PMC" + p["pmcid"], url, now()])

    # --- statements: verify verbatim presence, collect class per member
    cls = {}
    for fam, label, members, quote in STATEMENTS:
        txt = article_text(roots[fam])
        if quote not in txt:
            sys.exit(f"[06] ERROR: statement not found verbatim in the {fam} article: {quote!r}")
        for m in members:
            cls.setdefault(m, []).append((label, f"text: \"{quote}\""))
    # Li 2016 Table 2: class label in column 1 carries down over the paralogous-pair rows
    t2 = table_rows(roots["BURP"], "Table 2")
    group = ""
    for row in t2[1:]:
        if row and row[0]:
            group = row[0]
        pair = next((c for c in row if c.startswith("MtBURP")), "")
        for m in re.findall(r"MtBURP\d+", pair):
            label = {"RD22": "RD22-like", "BNM2": "BNM2-like", "USP": "USP-like", "PG1β": "PG1beta-like"}.get(group, group)
            cls.setdefault(m, []).append((label, f"Table 2 group '{group}', pair {pair}"))

    # --- Medicago UniProt annotation (the BLAST database's entries)
    mt = read_tsv(os.path.join(data, "uniprot", "mt_proteome.annot.tsv"))
    mt_acc = {r["Entry"]: r for r in mt}
    by_locus = {}
    for r in mt:
        for tok in re.findall(r"\b(?:MTR_|Medtr)(\d+g\d+)\b", r["Gene Names"] + " " + r["Gene Names (ordered locus)"]):
            by_locus.setdefault(tok, set()).add(r["Entry"])

    rows = []
    # BURP: Table 1 of Li 2016
    t1 = table_rows(roots["BURP"], "Table 1")
    for row in t1:
        if not row or not row[0].startswith("MtBURP"):
            continue
        name, seqid, length = row[0], row[1], row[2]
        m = re.match(r"Medtr(\d+g\d+)", seqid)
        locus = m.group(1) if m else ""
        accs = sorted(by_locus.get(locus, [])) if locus else []
        ulen = ";".join(mt_acc[a]["Length"] for a in accs)
        rows.append(["BURP", name, seqid, length, ";".join(accs),
                     ";".join(mt_acc[a]["Gene Names (ordered locus)"] or mt_acc[a]["Gene Names"] for a in accs),
                     ulen, "yes" if accs and length in ulen.split(";") else ("no" if accs else "unmapped"),
                     " | ".join(sorted({c for c, _ in cls.get(name, [])})) or "not stated in text",
                     " | ".join(s for _, s in cls.get(name, [])), fig.get(name, "")])
    # figure reading must agree with every class the text states
    for r in rows:
        if r[9] and r[8] != "not stated in text":
            fig_cls = r[10].split(" (")[0]
            if r[10] and any(fig_cls != c for c in r[8].split(" | ")):
                sys.exit(f"[06] ERROR: figure class {r[10]!r} disagrees with the text class {r[8]!r} for {r[1]}")
    n_checked = sum(1 for r in rows if r[10] and r[8] != "not stated in text")
    # YUC: Table 1 of Shao 2023 (Mt rows only)
    t1y = table_rows(roots["YUC"], "Table 1")
    yuc = [(r[1], r[0], r[2]) for r in t1y if len(r) > 2 and r[1].startswith("MtYUC")]
    mapped, job = uniprot_idmap([pid for _, pid, _ in yuc])
    for name, pid, length in yuc:
        hits = [h for h in mapped.get(pid, []) if h["accession"] in mt_acc]
        ulen = ";".join(h["length"] for h in hits)
        rows.append(["YUC", name, pid, length, ";".join(h["accession"] for h in hits),
                     ";".join(h["oln"] or h["gene_names"] for h in hits), ulen,
                     "yes" if hits and length in ulen.split(";") else ("no" if hits else "unmapped"),
                     " | ".join(sorted({c for c, _ in cls.get(name, [])})) or "not stated in text",
                     " | ".join(s for _, s in cls.get(name, [])), ""])
    cols = ["family", "published_name", "published_id", "published_length_aa", "uniprot_accessions",
            "uniprot_locus", "uniprot_length_aa", "length_match", "class_from_text", "class_source",
            "class_from_figure_read_by_eye"]
    with open(os.path.join(work, "medicago_family_names.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)

    # --- eggNOG seed orthologs (Medicago GenBank ids) of the region-1 and region-12 genes
    genes = [g for g in read_tsv(os.path.join(work, "genes.tsv")) if g["region_id"] in SYN_REGIONS]
    seeds = {}
    for g in genes:
        if g["em_target"].startswith("3880."):
            seeds[g["transcript_id"]] = g["em_target"].split(".", 1)[1].split(".")[0]
    smap, sjob = uniprot_idmap(sorted(set(seeds.values())))
    name_of = {}
    for r in rows:
        for a in r[4].split(";"):
            if a:
                name_of[a] = r[1]
    with open(os.path.join(work, "eggnog_seed_medicago.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["transcript_id", "em_target", "uniprot_accession", "uniprot_locus", "published_family_name"])
        for tid, pid in sorted(seeds.items(), key=lambda kv: int(re.sub(r"\D", "", kv[0].split(".")[0]))):
            hits = [h for h in smap.get(pid, []) if h["accession"] in mt_acc] or smap.get(pid, [])
            w.writerow([tid, "3880." + pid, ";".join(h["accession"] for h in hits),
                        ";".join(h["oln"] or h["gene_names"] for h in hits),
                        ";".join(sorted({name_of[h["accession"]] for h in hits if h["accession"] in name_of}))])

    with open(os.path.join(work, "family_sources.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "citation", "doi", "pmid", "pmcid", "fetched_from", "fetched_utc"])
        w.writerows(sources)
        w.writerow(["idmapping", "UniProt ID mapping EMBL-GenBank-DDBJ_CDS -> UniProtKB", "", "", "",
                    f"{IDMAP} jobs {job}, {sjob}", now()])
    n_burp = sum(r[0] == "BURP" for r in rows)
    n_yuc = sum(r[0] == "YUC" for r in rows)
    print(f"[06] MtBURP rows {n_burp} (mapped {sum(r[0]=='BURP' and r[4]!='' for r in rows)}), "
          f"MtYUC rows {n_yuc} (mapped {sum(r[0]=='YUC' and r[4]!='' for r in rows)}); "
          f"eggNOG Medicago seeds mapped for {len(seeds)} genes; figure classes given for {len(fig)} MtBURPs, "
          f"{n_checked} of them checked against a text-stated class (all agree)")


if __name__ == "__main__":
    main()
