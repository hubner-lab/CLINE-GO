"""gi_common.py -- shared helpers for the gene-identity scripts (stdlib only).

Imported by 07-10 through an explicit sys.path entry for this scripts directory (the scripts run under
`python3 -I`, which does not put the script directory on sys.path).
"""
import csv
import re
from collections import defaultdict

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
              "send", "evalue", "bitscore", "qlen", "slen", "qcovs", "qcovhsp"]
SETS = ("at_proteome", "at_swissprot", "mt_proteome", "gm_proteome")
# Col-0 AGI codes only: "AtNgNNNNN" not preceded by a word character (excludes e.g. "AXX17_At1g32710",
# an ecotype-specific locus name). Organellar codes (AtCg / AtMg) are parsed but flagged.
AGI_RE = re.compile(r"(?<![\w])At([1-5CM])g(\d{5})\b", re.I)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") if isinstance(r, dict) else r[i] for i, c in enumerate(cols)])


def agis(text):
    """Sorted, upper-case Col-0 AGI codes in a UniProt gene-name string."""
    return sorted({f"AT{m.group(1).upper()}G{m.group(2)}" for m in AGI_RE.finditer(text or "")})


def is_nuclear(agi):
    return agi[2] in "12345"


def agi_pos(agi):
    """(chromosome, number) of a nuclear AGI code, e.g. AT4G32890 -> ('4', 32890)."""
    return agi[2], int(agi[4:])


def panther(xref):
    """UniProt xref_panther 'PTHR45658;PTHR45658:SF18;' -> (families, subfamilies) as sorted lists."""
    ids = [x for x in (xref or "").split(";") if x]
    return sorted({x for x in ids if ":" not in x}), sorted({x for x in ids if ":" in x})


def panther_label(xref):
    fam, sub = panther(xref)
    return ";".join(sub) if sub else ";".join(fam)


def first_protein_name(name):
    """UniProt 'Protein names' holds the recommended name first, alternatives in parentheses."""
    depth, out = 0, []
    for ch in name or "":
        if ch == "(" and out and out[-1] == " ":
            depth += 1
        if depth == 0:
            out.append(ch)
        if ch == ")" and depth:
            depth -= 1
    return "".join(out).strip()


def gene_symbol(row):
    """Primary gene name: first token of UniProt 'Gene Names' unless that token is a locus/ORF name."""
    names = (row.get("Gene Names") or "").split()
    oln = set((row.get("Gene Names (ordered locus)") or "").replace("/", " ").replace(";", " ").split())
    orf = set((row.get("Gene Names (ORF)") or "").replace("/", " ").replace(";", " ").split())
    if names and names[0] not in oln and names[0] not in orf and not AGI_RE.fullmatch(names[0]):
        return names[0]
    return ""


def species_locus(set_name, row):
    """Locus key used to collapse UniProt entries of one gene: AGI group (Arabidopsis), else the ordered
    locus name (Medicago MTR_...), else the ORF name (Medicago MtrunA17_...), else the EnsemblPlants
    protein id (soybean KRH...), else an NCBI LOC gene name, else the accession."""
    if set_name.startswith("at_"):
        a = agis(row.get("Gene Names (ordered locus)")) or agis(row.get("Gene Names"))
        if a:
            return "/".join(a)
    for field in ("Gene Names (ordered locus)", "Gene Names (ORF)"):
        tok = (row.get(field) or "").split()
        if tok:
            return tok[0]
    ens = [t for t in (row.get("EnsemblPlants") or "").split(";") if t]
    if ens:
        return ens[0]
    loc = [t for t in (row.get("Gene Names") or "").split() if t.startswith("LOC")]
    if loc:
        return loc[0]
    return row["Entry"]


def merged_len(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def read_blast(path):
    """Per (query, subject): best HSP (max bitscore, then min E-value) plus coverage over all HSPs.
    Returns {query: [subject records sorted by bitscore desc]}."""
    hsps = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(BLAST_COLS):
                continue
            r = dict(zip(BLAST_COLS, f))
            for k in ("pident", "evalue", "bitscore", "qcovs", "qcovhsp"):
                r[k] = float(r[k])
            for k in ("length", "qstart", "qend", "sstart", "send", "qlen", "slen"):
                r[k] = int(r[k])
            hsps[(r["qseqid"], r["sseqid"])].append(r)
    by_q = defaultdict(list)
    for (q, s), hs in hsps.items():
        best = sorted(hs, key=lambda h: (-h["bitscore"], h["evalue"]))[0]
        rec = dict(best)
        rec["n_hsp"] = len(hs)
        rec["scov"] = merged_len([(h["sstart"], h["send"]) for h in hs]) / best["slen"]
        rec["qcov_merged"] = merged_len([(h["qstart"], h["qend"]) for h in hs]) / best["qlen"]
        by_q[q].append(rec)
    for q in by_q:
        by_q[q].sort(key=lambda r: (-r["bitscore"], r["evalue"], r["sseqid"]))
    return by_q


def collapse_loci(set_name, hits, annot):
    """Collapse subject entries of one query to loci. For Arabidopsis, entries sharing any AGI are merged
    (union-find), so an entry listing two AGIs joins the entries of either. Representative entry of a locus:
    highest bitscore, then reviewed, then accession. Returns locus records sorted by bitscore desc."""
    keys = {}
    for h in hits:
        row = annot.get(h["sseqid"], {"Entry": h["sseqid"]})
        keys[h["sseqid"]] = species_locus(set_name, row)
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for acc, key in keys.items():
        find(acc)
        if set_name.startswith("at_") and re.match(r"AT[1-5CM]G\d{5}", key):
            for a in key.split("/"):
                union(acc, "AGI:" + a)
        else:
            union(acc, "KEY:" + key)
    groups = defaultdict(list)
    for h in hits:
        groups[find(h["sseqid"])].append(h)
    loci = []
    for members in groups.values():
        members.sort(key=lambda h: (-h["bitscore"], 0 if annot.get(h["sseqid"], {}).get("Reviewed") == "reviewed"
                                    else 1, h["sseqid"]))
        rep = members[0]
        row = annot.get(rep["sseqid"], {})
        if set_name.startswith("at_"):
            ag = sorted({a for m in members for a in keys[m["sseqid"]].split("/")
                         if re.match(r"AT[1-5CM]G\d{5}$", a)})
            locus = "/".join(ag) if ag else keys[rep["sseqid"]]
        else:
            ag = []
            locus = keys[rep["sseqid"]]
        loci.append({"locus": locus, "agis": ag, "rep": rep, "row": row,
                     "accessions": [m["sseqid"] for m in members]})
    loci.sort(key=lambda L: (-L["rep"]["bitscore"], L["rep"]["evalue"], L["locus"]))
    return loci


def load_annot(path):
    return {r["Entry"]: r for r in read_tsv(path)}
