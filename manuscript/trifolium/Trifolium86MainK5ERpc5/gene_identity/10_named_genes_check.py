#!/usr/bin/env python3
"""10_named_genes_check.py -- curated identities of the named genes vs the re-run evidence, with a verdict.

Run with `python3 -I` (stdlib only; gi_common.py is loaded from this script's directory).

Usage:
  python3 -I 10_named_genes_check.py <data_dir> <out_dir> <named_genes_spec.tsv> \
      <region_gene_evidence.tsv> <evidence_rgt12.tsv> [tie_frac=0.90] [flank_bp=10000]

The spec gives one row per sub-claim: the level, the test and the expectation. The family, class and member rows
were written before the BLAST results of the named genes were summarised; context rows were added later and only
report. Expected loci are AGI codes or "symbol:<name>", resolved
through the UniProt gene names of UP000006548 (the script stops if a symbol maps to more or fewer than one AGI).

Tests (status in brackets):
  cdd_specific    the CD-Search specific hits include the expected domain accession
                  [reproduced | contradicted: other specific domain only | not reproduced: no specific hit]
  at_family       the best Arabidopsis locus shares a PANTHER family with the reference locus
                  [reproduced | contradicted | not reproduced: no PANTHER family on one side]
  at_member       expected AGI is the best locus [reproduced]; or lies in the near-tie set (bitscore >= tie_frac x
                  best) and is the synteny choice of 08 [reproduced]; or is in the near-tie set and the claim allows
                  a tie (tie_ok) [reproduced]; else in the near-tie set [near-tie]; else [contradicted]
  at_best         expected AGI is the best locus [reproduced]; in the near-tie set with tie_ok
                  [reproduced, order within near-tie differs]; in the near-tie set [near-tie]; else [contradicted]
  at_tie_set      every expected AGI is in the near-tie set [reproduced] else [contradicted]
  at_syntenic     expected AGI is the synteny choice of 08 [reproduced]; else [contradicted] (syntenic but another
                  candidate has more support / no collinear support / not a candidate)
  mt_best_member  among the published Medicago family members (06) hit by the gene, the expected member has the
                  highest bitscore [reproduced]; within tie_frac of it [near-tie]; else [contradicted];
                  member without a UniProt entry [not reproduced]
  mt_class        the class of the re-run's best member equals the expected class [reproduced] or not
                  [contradicted]; the class is taken from the source text, else from the figure reading of 06
                  (status suffixed "(figure read by eye)"); no class in either [not reproduced]
  mt_class_not    as mt_class, inverted (claim: the best member is NOT in that class)
  at_pident       context: % identity of the Trifolium protein to the expected locus vs the curated number
  at_report       context: best locus, near-ties and synteny choice
Verdict per gene (context rows never change it):
  disagrees                   a family- or class-level sub-claim is contradicted
  not reproduced              a family- or class-level sub-claim could not be tested, and none is contradicted
  agrees                      every family-, class- and member-level sub-claim is reproduced
  agrees at family level only family (and class) hold; a member-level sub-claim is contradicted, a near-tie or
                              not reproduced
TBLASTN (09) results are summarised for the genes that triggered it.
Writes <out_dir>/named_genes_check.tsv and <data_dir>/work/named_genes_subclaims.tsv
"""
import os
import re
import sys
from collections import defaultdict

sys.dont_write_bytecode = True  # keep the scripts directory free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gi_common as gc  # noqa: E402

def is_ok(status):
    return status.startswith("reproduced")


def is_contra(status):
    return status.startswith("contradicted")


def f1(x):
    return f"{float(x):.1f}"


def short_name(r):
    return r["gene_symbol"] or r["protein_name"][:30]


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    data, out_dir, spec_path, ev_path, rgt12_path = sys.argv[1:6]
    tie_frac = float(sys.argv[6]) if len(sys.argv) > 6 else 0.90
    flank = int(sys.argv[7]) if len(sys.argv) > 7 else 10000
    work = os.path.join(data, "work")

    spec = gc.read_tsv(spec_path)
    curated = {r["gene_id"]: r["identity"] for r in gc.read_tsv(ev_path)}
    for r in gc.read_tsv(rgt12_path):
        if r["identity"]:
            curated[r["transcript"] or r["gene_id"] + ".t1"] = r["identity"]
    genes = {r["transcript_id"]: r for r in gc.read_tsv(os.path.join(work, "genes.tsv"))}
    annot = gc.load_annot(os.path.join(data, "uniprot", "at_proteome.annot.tsv"))

    # symbol -> AGI and AGI -> PANTHER families, from the UniProt annotation of UP000006548
    sym2agi, agi_panther = defaultdict(set), defaultdict(set)
    for r in annot.values():
        a = gc.agis(r.get("Gene Names (ordered locus)")) or gc.agis(r.get("Gene Names"))
        for tok in (r.get("Gene Names") or "").split():
            sym2agi[tok].update(a)
        for x in a:
            agi_panther[x] |= set(gc.panther(r.get("PANTHER"))[0])

    def resolve(token):
        if token.startswith("symbol:"):
            s = sym2agi.get(token[7:], set())
            if len(s) != 1:
                sys.exit(f"[10] ERROR: {token} resolves to {sorted(s)}")
            return next(iter(s))
        return token.upper()

    at = defaultdict(list)
    for r in gc.read_tsv(os.path.join(work, "arabidopsis_loci.tsv")):
        r["agiset"] = {a for a in r["agis"].split(";") if a}
        at[r["transcript_id"]].append(r)
    syn = defaultdict(list)
    for r in gc.read_tsv(os.path.join(work, "microsynteny_candidates.tsv")):
        syn[r["transcript_id"]].append(r)
    cds = defaultdict(list)
    cds_path = os.path.join(data, "cdsearch", "cdsearch_hits.tsv")
    cds_ran = os.path.exists(cds_path)
    if cds_ran:
        for r in gc.read_tsv(cds_path):
            cds[r["query_id"]].append(r)
    fam_rows = gc.read_tsv(os.path.join(work, "medicago_family_names.tsv"))
    acc2member = {a: r["published_name"] for r in fam_rows for a in r["uniprot_accessions"].split(";") if a}
    member_info = {r["published_name"]: r for r in fam_rows}
    mt_blast = gc.read_blast(os.path.join(data, "blast", "mt_proteome.blastp.tsv"))
    mt_top = {}
    for r in gc.read_tsv(os.path.join(work, "blast_top_hits.tsv")):
        if r["database"] == "mt_proteome" and r["rank"] == "1":
            mt_top[r["transcript_id"]] = r
    seeds = {r["transcript_id"]: r for r in gc.read_tsv(os.path.join(work, "eggnog_seed_medicago.tsv"))}

    def members_ranked(q, prefix):
        best = {}
        for h in mt_blast.get(q, []):
            m = acc2member.get(h["sseqid"])
            if m and m.startswith(prefix) and (m not in best or h["bitscore"] > best[m]["bitscore"]):
                best[m] = h
        return sorted(best.items(), key=lambda kv: -kv[1]["bitscore"])

    # TBLASTN summaries (09)
    tbl = {}
    trig_path = os.path.join(data, "tblastn", "triggers.tsv")
    if os.path.exists(trig_path):
        spans, tid_re = {}, re.compile(r'transcript_id "([^"]+)"')
        with open(os.path.join(data, "lab", "Trifolium_genes.gtf")) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) > 8 and f[2] == "CDS":
                    t = tid_re.search(f[8]).group(1)
                    c, s, e = spans.get(t, (f[0], int(f[3]), int(f[4])))
                    spans[t] = (f[0], min(s, int(f[3])), max(e, int(f[4])))
        for t in gc.read_tsv(trig_path):
            if t["trigger"] != "yes":
                continue
            q = t["gene"]
            chrom, rng = t["window"].split(":")
            w0 = int(rng.split("-")[0])
            hsps = []
            with open(os.path.join(data, "tblastn", f"{q}.tblastn.tsv")) as fh:
                for line in fh:
                    x = line.rstrip("\n").split("\t")
                    s, e = sorted((w0 + int(x[8]) - 1, w0 + int(x[9]) - 1))
                    over = sorted(tt for tt, (c, a, b) in spans.items() if c == chrom and a <= e and b >= s)
                    hsps.append({"qs": int(x[6]), "qe": int(x[7]), "s": s, "e": e, "pid": float(x[2]),
                                 "genes": over, "qlen": int(x[12])})
            _, gs, ge = spans[q]
            qlen = hsps[0]["qlen"] if hsps else 0
            own = [h for h in hsps if q in h["genes"]]
            free = [h for h in hsps if not h["genes"]]
            near_free = [h for h in free if h["e"] >= gs - flank and h["s"] <= ge + flank]
            own_txt = ", ".join(f"aa {h['qs']}-{h['qe']}" for h in own) or "nothing"
            own_pct = 100 * gc.merged_len([(h["qs"], h["qe"]) for h in own]) / qlen if qlen else 0
            first = min(hsps, key=lambda h: h["qs"]) if hsps else None
            tbl[q] = (f"{t['accession']} ({qlen} aa) vs {t['window']}: {len(hsps)} HSPs, {len(free)} outside any "
                      f"gene model; the {q.split('.')[0]} model matches {own_txt} ({own_pct:.0f}%); "
                      f"{len(near_free)} HSP within {flank // 1000} kb of it outside a model "
                      f"({'possible missed exon' if near_free else 'no sign of a truncated model'}); "
                      + (f"the most N-terminal match in the window starts at aa {first['qs']} "
                         f"({', '.join(g.split('.')[0] for g in first['genes']) or 'no model'})" if first else ""))

    sub_rows, out_rows = [], []
    by_gene = defaultdict(list)
    for r in spec:
        by_gene[r["gene"]].append(r)
    for q, claims in by_gene.items():
        loci = at.get(q, [])
        best = loci[0] if loci else None
        top = float(best["bitscore"]) if best else 0.0
        ties = [r for r in loci if float(r["bitscore"]) >= tie_frac * top] if best else []
        syn_rows = syn.get(q, [])
        choice = next((r for r in syn_rows if r["synteny_choice"] == "yes"), None)
        results = []
        for c in claims:
            test, exp, tie_ok = c["test"], c["expected"], c["tie_ok"] == "yes"
            status, note, short = "not reproduced", "", ""
            if test == "cdd_specific":
                hits = [h for h in cds.get(q, []) if h["hit_type"] == "specific"]
                mine = [h for h in hits if h["accession"] == exp]
                if not cds_ran:
                    note = short = "CD-Search output absent"
                elif mine:
                    status = "reproduced"
                    note = f"CD-Search specific {exp} {mine[0]['short_name']} aa {mine[0]['from']}-{mine[0]['to']}"
                    short = f"CD-Search {exp} {mine[0]['short_name']}"
                elif hits:
                    status = "contradicted"
                    note = short = "CD-Search specific: " + ", ".join(f"{h['accession']} {h['short_name']}" for h in hits)
                else:
                    note = short = "no CD-Search specific hit"
            elif test.startswith("at_") and best is None:
                note = short = "no Arabidopsis hit"
            elif test == "at_family":
                ref = resolve(exp)
                ref_fam = agi_panther.get(ref, set())
                best_fam = set().union(*(agi_panther.get(a, set()) for a in best["agiset"])) if best["agiset"] else set()
                if not ref_fam or not best_fam:
                    note = short = f"PANTHER family missing (best {sorted(best_fam)}, reference {sorted(ref_fam)})"
                elif ref_fam & best_fam:
                    status = "reproduced"
                    fam = "/".join(sorted(ref_fam & best_fam))
                    note = f"best {best['locus']} {short_name(best)} in {fam}, the PANTHER family of {ref}"
                    short = f"best Arabidopsis hit {best['locus']} {short_name(best)} ({fam})"
                else:
                    status = "contradicted"
                    note = short = (f"best {best['locus']} {short_name(best)} in {'/'.join(sorted(best_fam))}; "
                                    f"{ref} is in {'/'.join(sorted(ref_fam))}")
            elif test in ("at_member", "at_best"):
                ref = resolve(exp)
                hit = next((r for r in loci if ref in r["agiset"]), None)
                in_tie = hit is not None and float(hit["bitscore"]) >= tie_frac * top
                lab = f"{hit['locus']} {short_name(hit)}" if hit else ref
                if hit is best:
                    status = "reproduced"
                    runner = f"; next {ties[1]['locus']} {short_name(ties[1])} at {ties[1]['frac_of_best_bitscore']}" \
                        if len(ties) > 1 else ""
                    note = f"best locus {lab} ({f1(best['pident'])}%, {best['bitscore']} bits){runner}"
                    short = f"{lab} best ({f1(best['pident'])}%{runner})"
                elif in_tie and test == "at_member" and choice and ref in choice["candidate_locus"].split("/"):
                    status = "reproduced"
                    note = (f"{lab} at {hit['frac_of_best_bitscore']} of best {best['locus']} {short_name(best)}; "
                            f"near-tie resolved by synteny")
                    short = f"{lab} near-tie with {short_name(best)}, chosen by synteny"
                elif in_tie and tie_ok:
                    status = "reproduced" if test == "at_member" else "reproduced, order within near-tie differs"
                    note = short = f"{lab} at {hit['frac_of_best_bitscore']} of best {best['locus']} {short_name(best)}"
                elif in_tie:
                    status = "near-tie"
                    note = short = f"{lab} at {hit['frac_of_best_bitscore']} of best {best['locus']} {short_name(best)}"
                else:
                    status = "contradicted"
                    note = short = (f"{lab} only at {hit['frac_of_best_bitscore']} of best {best['locus']} "
                                    f"{short_name(best)}" if hit else f"no hit to {ref}; best {best['locus']} {short_name(best)}")
            elif test == "at_tie_set":
                refs = [resolve(x) for x in exp.split(";")]
                parts, missing = [], []
                for ref in refs:
                    hit = next((r for r in loci if ref in r["agiset"]), None)
                    if hit is None or float(hit["bitscore"]) < tie_frac * top:
                        missing.append(ref)
                    if hit:
                        parts.append(f"{short_name(hit)} {hit['bitscore']} bits ({hit['frac_of_best_bitscore']})")
                status = "contradicted" if missing else "reproduced"
                note = short = " vs ".join(parts) + (f"; outside the near-tie set: {', '.join(missing)}" if missing else "")
            elif test == "at_syntenic":
                ref = resolve(exp)
                row = next((r for r in syn_rows if ref in r["candidate_locus"].split("/")), None)
                others = ", ".join(f"{r['gene_symbol'] or r['candidate_locus']} {r['support_n']}"
                                   for r in syn_rows if r is not row and r.get("candidate_locus"))
                if row is None:
                    status = "contradicted"
                    note = short = f"{ref} not a synteny candidate (below {tie_frac} of the best bitscore)"
                elif row["synteny_choice"] == "yes":
                    status = "reproduced"
                    note = (f"synteny choice {row['candidate_locus']} {row['gene_symbol']} (support "
                            f"{row['support_n']} genes / {row['support_loci_n']} Arabidopsis loci: "
                            f"{row['supporting_genes']})"
                            + (f"; other candidates (support genes): {others}" if others else ""))
                    short = (f"synteny picks {row['gene_symbol'] or row['candidate_locus']} ({row['support_n']} genes"
                             f" / {row['support_loci_n']} loci" + (f"; others {others}" if others else "") + ")")
                elif row["syntenic"] == "yes":
                    status = "contradicted"
                    note = short = (f"{ref} support {row['support_n']}, {choice['candidate_locus']} has more "
                                    f"({choice['support_n']})")
                else:
                    status = "contradicted"
                    note = short = f"{ref} has no collinear support ({row['support_n']} genes)"
            elif test == "mt_best_member":
                prefix = re.match(r"(Mt[A-Z]+)", exp).group(1)
                ranked = members_ranked(q, prefix)
                if exp not in member_info or not member_info[exp]["uniprot_accessions"]:
                    note = short = f"{exp} has no UniProt entry in UP000002051"
                elif not ranked:
                    note = short = f"no published {prefix} member among the Medicago hits"
                else:
                    m1, h1 = ranked[0]
                    mine = dict(ranked).get(exp)
                    if m1 == exp:
                        status = "reproduced"
                    elif mine and mine["bitscore"] >= tie_frac * h1["bitscore"]:
                        status = "near-tie"
                    else:
                        status = "contradicted"
                    note = "members ranked: " + "; ".join(f"{m} {h['pident']:.1f}% {h['bitscore']:.0f} bits"
                                                          for m, h in ranked[:3])
                    short = (f"best member {m1} ({h1['pident']:.1f}%, {h1['bitscore']:.0f} bits"
                             + (f"; next {ranked[1][0]} {ranked[1][1]['bitscore']:.0f}" if len(ranked) > 1 else "") + ")")
            elif test in ("mt_class", "mt_class_not"):
                prefix = "MtYUC" if exp.startswith("clade") else "MtBURP"
                ranked = members_ranked(q, prefix)
                if not ranked:
                    note = short = "no published member hit"
                else:
                    m1 = ranked[0][0]
                    cls = member_info[m1]["class_from_text"]
                    figc = member_info[m1].get("class_from_figure_read_by_eye", "")
                    if cls != "not stated in text":
                        src, tag = "text", ""
                    elif figc:
                        cls, src, tag = figc.split(" (")[0], "figure, read by eye", " (figure read by eye)"
                    else:
                        src = ""
                    if not src:
                        note = short = f"{m1}: no class stated in the source text"
                    else:
                        has = exp in cls
                        if test == "mt_class":
                            status = ("reproduced" if has else "contradicted") + tag
                        else:
                            status = ("contradicted" if has else "reproduced") + tag
                        note = (f"{m1}: {cls} ({member_info[m1]['class_source'][:100]})" if src == "text"
                                else f"{m1}: {figc}; the text states no class for {m1}")
                        short = f"{m1} {cls} in the source {'text' if src == 'text' else 'figures (read by eye)'}"
            elif test == "at_pident":
                ref, val = exp.split(":")
                hit = next((r for r in loci if ref in r["agiset"]), None)
                if hit:
                    d = float(hit["pident"]) - float(val)
                    status = "reported" + (" (differs by > 5 points)" if abs(d) > 5 else "")
                    note = short = (f"Trifolium protein vs {hit['locus']}: {f1(hit['pident'])}% identity, query "
                                    f"cover {hit['qcovs']}%, subject cover {hit['scov']}% (curated ~{val}%, "
                                    f"query: {c['curated_query']})")
                else:
                    note = short = f"no hit to {ref}"
            elif test == "at_report":
                tie_txt = ", ".join(f"{r['locus']} {short_name(r)} {r['frac_of_best_bitscore']}" for r in ties[1:]) or "none"
                syn_txt = (f"synteny choice {choice['candidate_locus']} {choice['gene_symbol']} (support "
                           f"{choice['support_n']} genes / {choice['support_loci_n']} loci: "
                           f"{choice['supporting_genes']})") if choice else "no synteny call"
                status = "reported"
                note = short = f"best {best['locus']} {short_name(best)}; near-ties {tie_txt}; {syn_txt}"
            results.append((c, status, note, short))
            sub_rows.append({"gene": q, "claim_id": c["claim_id"], "level": c["level"], "test": test,
                             "expected": exp, "manuscript_claim": c["manuscript_claim"], "status": status,
                             "evidence": note})

        def st(level):
            return [s for c, s, _, _ in results if c["level"] == level]
        fam_s, cls_s, mem_s = st("family"), st("class"), st("member")
        if any(is_contra(s) for s in fam_s + cls_s):
            verdict = "disagrees"
        elif any(s.startswith("not reproduced") for s in fam_s + cls_s):
            verdict = "not reproduced"
        elif all(is_ok(s) for s in fam_s + cls_s + mem_s):
            verdict = "agrees"
        else:
            verdict = "agrees at family level only"
        failed = [(c, s, sh) for c, s, _, sh in results if c["level"] != "context" and not is_ok(s)]
        held = [sh for c, s, _, sh in results if c["level"] != "context" and is_ok(s)]
        reason = "; ".join(held)
        if failed:
            reason = ("; ".join(f"{c['manuscript_claim']}: {s} ({sh})" for c, s, sh in failed)
                      + (f". Reproduced: {reason}" if reason else ""))

        if syn_rows and not choice:
            syn_txt = "no synteny call; candidates (support genes/loci): " + ", ".join(
                f"{r['candidate_locus']} {r['gene_symbol']} {r['support_n']}/{r['support_loci_n']}"
                + (f" [{r['supporting_genes']}]" if r.get("supporting_genes") else "") for r in syn_rows)
        elif choice:
            syn_txt = (f"{choice['candidate_locus']} {choice['gene_symbol']} (support {choice['support_n']} genes / "
                       f"{choice['support_loci_n']} Arabidopsis loci); candidates (genes/loci): "
                       + ", ".join(f"{r['gene_symbol'] or r['candidate_locus']} {r['support_n']}/{r['support_loci_n']}"
                                   for r in syn_rows))
        else:
            syn_txt = "not run (region not in 08)"
        mt = mt_top.get(q)
        mt_rank = members_ranked(q, "MtBURP") or members_ranked(q, "MtYUC")
        mt_txt = (f"top {mt['accession']} {mt['locus']}" + (f" = {mt['published_family_name']}" if mt["published_family_name"] else "")
                  + f" {mt['pident']}% {mt['bitscore']} bits") if mt else "no hit"
        if mt_rank:
            mt_txt += "; published members: " + ", ".join(f"{m} {h['bitscore']:.0f}" for m, h in mt_rank[:3])
        seed = seeds.get(q)
        basis = ""
        if seed:
            basis = (f"eggNOG seed {seed['em_target']} = UniProt {seed['uniprot_accession']} {seed['uniprot_locus']}"
                     + (f" = {seed['published_family_name']}" if seed["published_family_name"] else ""))
        else:
            basis = f"eggNOG seed {genes[q]['em_target']}"
        out_rows.append({
            "region_id": genes[q]["region_id"], "gene": q,
            "manuscript_claim": "; ".join(dict.fromkeys(c["manuscript_claim"] for c, *_ in results
                                                         if c["level"] != "context")),
            "curated_identity": curated.get(q, ""), "curated_basis": f"{claims[0]['curated_query']} [{basis}]",
            "rerun_arabidopsis_best": (f"{best['locus']} {short_name(best)} ({best['accession']}) {f1(best['pident'])}% id, "
                                       f"qcov {best['qcovs']}%, scov {best['scov']}%, {best['bitscore']} bits, "
                                       f"PANTHER {best['panther']}") if best else "none",
            "rerun_arabidopsis_near_ties": ", ".join(f"{r['locus']} {short_name(r)} {r['bitscore']} bits "
                                                     f"({r['frac_of_best_bitscore']})" for r in ties[1:]) or "none",
            "rerun_synteny": syn_txt, "rerun_medicago": mt_txt,
            "rerun_cdd_specific": (", ".join(f"{h['accession']} {h['short_name']} {h['from']}-{h['to']}"
                                             for h in cds.get(q, []) if h["hit_type"] == "specific") or "none")
            if cds_ran else "CD-Search not run",
            "rerun_tblastn": tbl.get(q, "not triggered"),
            "verdict": verdict, "reason": reason,
            "context_notes": " | ".join(n for c, s, n, _ in results if c["level"] == "context")})

    cols = ["region_id", "gene", "manuscript_claim", "curated_identity", "curated_basis", "rerun_arabidopsis_best",
            "rerun_arabidopsis_near_ties", "rerun_synteny", "rerun_medicago", "rerun_cdd_specific", "rerun_tblastn",
            "verdict", "reason", "context_notes"]
    gc.write_tsv(os.path.join(out_dir, "named_genes_check.tsv"), cols, out_rows)
    gc.write_tsv(os.path.join(work, "named_genes_subclaims.tsv"),
                 ["gene", "claim_id", "level", "test", "expected", "manuscript_claim", "status", "evidence"], sub_rows)
    for r in out_rows:
        print(f"[10] {r['gene']}\t{r['verdict']}")


if __name__ == "__main__":
    main()
