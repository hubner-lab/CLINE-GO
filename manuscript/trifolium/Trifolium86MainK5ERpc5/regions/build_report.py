#!/usr/bin/env python3
"""Region-QC report page v2 for the 15 Trifolium86MainK5ERpc5 (5 PCs) regions (plan ~/.claude/plans/peaceful-kindling-micali.md,
Phase 7, 2026-10-16). COPY of the PC3 run's regionsER/build_report.py (plan wild-cooking-kite.md).

Reads (nothing is recomputed):
  report/region_genes.tsv            from prep_regionplots.py
  figures/input/regionplot_calls.tsv from prep_regionplots.py
  figures/regionplot_<order>_<region>.png  from figures/run_regionplots.sh
  tables/region_sig_snps.tsv, tables/sig_snp_patterns.tsv
  REGIONS_v2.md                      the overview table rows (support verdict) and section headlines
Writes
  report/regions_qc.html  (page body for the Artifact tool; images referenced as img/regionplot_<order>.png)
  report/img/regionplot_<order>.png  copies of the plots
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/build_report.py   (from the repo root)
"""
import csv
import html
import os
import re
import shutil
import sys

B = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(B, "report")
T = os.path.join(B, "tables")
THR = {"EMMAX": "2.82e-7", "RDA": "1.03e-5"}


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def e(x):
    return html.escape(str(x))


calls = rows(os.path.join(B, "figures", "input", "regionplot_calls.tsv"))
genes = rows(os.path.join(R, "region_genes.tsv"))
sigs = rows(os.path.join(T, "region_sig_snps.tsv"))
pats = {(r["region_id"], r["snp_id"]): r for r in rows(os.path.join(T, "sig_snp_patterns.tsv"))}
md = open(os.path.join(B, "REGIONS_v2.md")).read().split("\n")

# overview rows and headlines from REGIONS_v2.md
over, head = {}, {}
for line in md:
    m = re.match(r"\| (\d+) \| (\S+) \| (\d+) \| (\d+) \| (\d+) \| ([^|]+) \| ([^|]+) \| ([^|]+) \|$", line)
    if m:
        over[m.group(2)] = [x.strip() for x in m.groups()]
    m = re.match(r"## (\d+)\. (\S+) — (.+)$", line)
    if m:
        head[m.group(2)] = m.group(3)
assert len(over) == 15 and len(head) == 15, (len(over), len(head))

# checks: every in-region gene is a key or host gene of its region, every outside gene is a neighbourhood gene
nbh = {}
for r in rows(os.path.join(T, "neighbourhood.tsv")):
    for rid in r["region_ids"].split(","):
        nbh.setdefault(rid, set()).add(r["gene_id"].split(".")[0])
ev = {(r["region_id"], r["gene_id"].split(".")[0]) for r in rows(os.path.join(T, "region_gene_evidence.tsv"))}
for g in genes:
    k = (g["region_id"], g["gene_id"])
    ok = (k in ev) if g["where"] != "outside (guess)" else (g["gene_id"] in nbh.get(g["region_id"], set()))
    if not ok:
        sys.exit(f"ERROR: {g['gene_id']} ({g['where']}) does not belong to {g['region_id']}")


def support_class(s):
    s = s.lower()
    if s.startswith("strong"):
        return "s-strong"
    if s.startswith("moderate"):
        return "s-mod"
    if s.startswith("weak to"):
        return "s-weakmod"
    if s.startswith("weak"):
        return "s-weak"
    return "s-none"


def support_word(s):
    m = re.match(r"(strong class|moderate|weak to moderate|weak|none)", s)
    return m.group(1) if m else s


def pmid_links(cell):
    def rep(m):
        p = m.group(0)
        return f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/">{p}</a>'
    return re.sub(r"\b\d{7,8}\b", rep, e(cell))


def fmt_p(x):
    if not x or x in ("NA", "nan"):
        return "—"
    v = float(x)
    return f"{v:.1e}" if v < 1e-3 else f"{v:.3f}".rstrip("0")


os.makedirs(os.path.join(R, "img"), exist_ok=True)
MANH = os.path.join(os.getcwd(), "Trifolium86MainK5ERpc5_results/GEA/plots/manhattan/combined/manhattan_combined_K5.png")
if not os.path.getsize(MANH):
    sys.exit(f"ERROR: empty {MANH}")
shutil.copyfile(MANH, os.path.join(R, "img", "manhattan_combined.png"))
sections, toc = [], []
for c in calls:
    o, rid = int(c["order"]), c["region_id"]
    png = os.path.join(B, "figures", f"regionplot_{o}_{rid}.png")
    if not os.path.getsize(png):
        sys.exit(f"ERROR: empty plot {png}")
    shutil.copyfile(png, os.path.join(R, "img", f"regionplot_{o}.png"))
    ov = over[rid]
    sup = ov[7]
    # signal line per significant SNP
    sl = []
    for s in [x for x in sigs if x["region_id"] == rid]:
        p = pats[(rid, s["snp_id"])]
        em = f"EMMAX {fmt_p(s['EMMAX_minp'])} ({s['EMMAX_pred']})" if "EMMAX" in s["sig_methods"] else \
            f"EMMAX {fmt_p(s['EMMAX_minp'])}, n.s."
        rd = f"RDA {fmt_p(s['RDA_p'])}" + ("" if "RDA" in s["sig_methods"] else ", n.s.")
        uniq = "unique vector" if p["locus_unique"] == "True" else f"vector shared by {p['pattern_n']} SNPs"
        sl.append(f"<li><span class=mono>{e(s['snp_id'])}</span> {e(s['sig_methods'])} &middot; {e(em)} &middot; "
                  f"{e(rd)} &middot; carriers <span class=mono>{e(p['carrier_sites'])}</span> &middot; {e(uniq)}"
                  f" &middot; in {e(s['gene_hit'] or '—')}</li>")
    gr = []
    for g in [x for x in genes if x["region_id"] == rid]:
        cls = {"key": "w-key", "sig-SNP host": "w-host", "outside (guess)": "w-out"}[g["where"]]
        tag = {"key": "key", "sig-SNP host": "SNP host", "outside (guess)": "outside"}[g["where"]]
        body_prom = "—" if g["where"] == "outside (guess)" else f"{g['snps_body']} / {g['snps_promoter']}"
        pe = "—" if g["where"] == "outside (guess)" else \
            f"{fmt_p(g['best_p_EMMAX_body'])} / {fmt_p(g['best_p_EMMAX_promoter'])}"
        lab = e(g["label"].rstrip("*")) if g["label"] else '<span class=muted>not on plot</span>'
        gr.append(f"<tr class={cls}><td><span class='tag {cls}'>{tag}</span></td><td class=mono>{e(g['gene_id'])}</td>"
                  f"<td>{lab}</td><td>{e(g['identity'])}</td><td class=num>{e(g['distance'])}</td>"
                  f"<td class=num>{body_prom}</td><td class=num>{pe}</td><td>{pmid_links(g['pmids']) or '—'}</td></tr>")
    notes = []
    if rid.startswith("2_176097944"):
        notes.append("The six EMMAX SNPs (bio_15) share one coastal site-indicator vector (28955 + 28960) with 135 "
                     "significant SNPs on all six chromosomes; only the RDA SNP 2:176317977 is locus-specific.")
    if rid.startswith("4_100931231"):
        notes.append("Same highland haplotype as region 9 (4:101438316, r2 0.77): read the two together.")
    if rid.startswith("4_101307505"):
        notes.append("New at 5 PCs. Same carriers as region 8's SNP 4:101062042, 376 kb left (r2 0.77): one haplotype "
                     "tagged twice.")
    if rid.startswith("6_116698620"):
        notes.append("6:116959134 and 6:116969913 share their vectors with 6:109080530-588, 7.9 Mb away (PC3 region 11, "
                     "dropped at 5 PCs because its EMMAX p rose above the threshold).")
    if rid in ("2_145944878-146206511", "5_8217606-8479228"):
        notes.append("New at 5 PCs.")
    note_html = "".join(f"<p class=note>{e(n)}</p>" for n in notes)
    toc.append(f"<tr><td class=num>{o}</td><td><a href='#r{o}' class=mono>{e(rid)}</a></td><td class=num>{ov[2]}</td>"
               f"<td class=num>{ov[3]}</td><td class=num>{ov[4]}</td><td class=num>{e(ov[5])}</td>"
               f"<td>{e(ov[6])}</td><td><span class='pill {support_class(sup)}'>{e(support_word(sup))}</span> "
               f"<span class=muted>{e(sup[len(support_word(sup)):].strip(' ()'))}</span></td></tr>")
    sections.append(f"""
<section id="r{o}" class=region>
  <header class=rhead>
    <span class=rnum>{o}</span>
    <div class=rtitle><h2 class=mono>{e(rid)}</h2><p>{e(head[rid])}</p></div>
    <span class='pill {support_class(sup)}'>{e(support_word(sup))}</span>
  </header>
  <ul class=signal>{''.join(sl)}</ul>
  {note_html}
  <figure><img src="img/regionplot_{o}.png" alt="topr region plot of {e(rid)}, region plus 100 kb each side" loading="lazy" width="4200" height="1500"></figure>
  <div class=tablewrap><table class=genes>
    <thead><tr><th></th><th>Gene</th><th>Plot label</th><th>Identity / guess</th><th>Position</th>
    <th>SNPs body / promoter</th><th>Best EMMAX p body / promoter</th><th>Literature (PMID)</th></tr></thead>
    <tbody>{''.join(gr)}</tbody></table></div>
</section>""")

page = f"""<title>Trifolium PC5 Region Check</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Sans+Condensed:wght@600&display=swap">
<style>
/* Layout: one narrow reading column; overview table, then one block per region: header, SNP list, plot, gene table. */
:root {{
  --bg: #f6f8f7; --surface: #ffffff; --fg: #17201d; --muted: #5d6b66; --line: #d7dfdc;
  --accent: #0b775e; --emmax: #8b0000; --rda: #36648b;
  --strong: #0b775e; --mod: #3b7f8f; --weak: #a8690c; --none: #8a5a5a;
  --tag-key: #e2f1ec; --tag-host: #e8eef6; --tag-out: #f5efe2;
  --f-body: "IBM Plex Sans", system-ui, sans-serif; --f-head: "IBM Plex Sans Condensed", "IBM Plex Sans", system-ui, sans-serif;
  --f-mono: "IBM Plex Mono", ui-monospace, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #111816; --surface: #18211e; --fg: #e3ebe8; --muted: #9aaba5; --line: #2c3934;
  --accent: #4fc3a1; --emmax: #e0605a; --rda: #79a8d6;
  --strong: #4fc3a1; --mod: #6fb6c6; --weak: #e0a94a; --none: #c28f8f;
  --tag-key: #1d3a31; --tag-host: #1f2c3b; --tag-out: #3a3222; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #111816; --surface: #18211e; --fg: #e3ebe8; --muted: #9aaba5; --line: #2c3934;
  --accent: #4fc3a1; --emmax: #e0605a; --rda: #79a8d6;
  --strong: #4fc3a1; --mod: #6fb6c6; --weak: #e0a94a; --none: #c28f8f;
  --tag-key: #1d3a31; --tag-host: #1f2c3b; --tag-out: #3a3222; color-scheme: dark; }}
body {{ background: var(--bg); color: var(--fg); font: 15px/1.5 var(--f-body); }}
.wrap {{ max-width: 1180px; margin: 0 auto; padding-inline: 20px; padding-block: 28px 64px; display: grid; gap: 28px; }}
h1 {{ font: 600 1.9rem/1.15 var(--f-head); margin: 0; text-wrap: balance; }}
h2 {{ font-size: 1.05rem; margin: 0; font-weight: 500; }}
a {{ color: var(--accent); }}
.mono {{ font-family: var(--f-mono); }}
.muted {{ color: var(--muted); }}
.num {{ font-variant-numeric: tabular-nums; white-space: nowrap; }}
.wrap > *, .region > * {{ min-width: 0; }}
.signal li, .rtitle h2, .rtitle p {{ overflow-wrap: anywhere; }}
.intro {{ display: grid; gap: 10px; max-width: 75ch; }}
.intro p {{ margin: 0; }}
.legend {{ display: flex; flex-wrap: wrap; gap: 8px 22px; font-size: .9rem; color: var(--muted); }}
.legend span {{ display: inline-flex; align-items: center; gap: 8px; }}
.dash {{ display: inline-block; width: 26px; border-top: 2px dashed currentColor; }}
.dot {{ display: inline-block; width: 26px; border-top: 2px dotted var(--muted); }}
.l-emmax {{ color: var(--emmax); }} .l-rda {{ color: var(--rda); }}
.tablewrap {{ overflow-x: auto; min-width: 0; border: 1px solid var(--line); border-radius: 6px; background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: .86rem; }}
th, td {{ text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--line); vertical-align: top; }}
th {{ font-weight: 600; color: var(--muted); font-size: .76rem; letter-spacing: .03em; text-transform: uppercase; }}
tbody tr:last-child td {{ border-bottom: 0; }}
.pill {{ display: inline-block; padding: 1px 9px; border-radius: 999px; font-size: .78rem; font-weight: 600;
  border: 1px solid currentColor; white-space: nowrap; }}
.s-strong {{ color: var(--strong); }} .s-mod {{ color: var(--mod); }} .s-weakmod {{ color: var(--weak); }}
.s-weak {{ color: var(--weak); }} .s-none {{ color: var(--none); }}
.region {{ display: grid; gap: 12px; padding-top: 20px; border-top: 2px solid var(--line); }}
.rhead {{ display: flex; align-items: center; gap: 14px; flex-wrap: wrap; }}
.rnum {{ font: 600 1.6rem/1 var(--f-head); color: var(--accent); min-width: 1.6em; }}
.rtitle {{ flex: 1 1 300px; min-width: 0; }}
.rtitle p {{ margin: 2px 0 0; color: var(--muted); }}
.signal {{ margin: 0; padding-left: 20px; font-size: .86rem; display: grid; gap: 2px; }}
.note {{ margin: 0; font-size: .86rem; padding: 6px 10px; border-left: 3px solid var(--weak); background: var(--surface); }}
figure {{ margin: 0; background: #ffffff; border: 1px solid var(--line); border-radius: 6px; overflow: hidden; }}
figure img {{ display: block; width: 100%; height: auto; }}
.tag {{ font-size: .72rem; padding: 1px 6px; border-radius: 4px; white-space: nowrap; }}
.tag.w-key {{ background: var(--tag-key); }} .tag.w-host {{ background: var(--tag-host); }} .tag.w-out {{ background: var(--tag-out); }}
tr.w-out td {{ color: var(--muted); }}
@media (max-width: 600px) {{ .wrap {{ padding-inline: 16px; }} h1 {{ font-size: 1.5rem; }} }}
</style>
<div class=wrap>
<header class=intro>
  <h1>Trifolium EMMAX + RDA regions at 5 PCs: plots and main genes</h1>
  <p>The 15 regions of Trifolium86MainK5ERpc5: EMMAX (5 PCs, Bonferroni 0.05) and RDA (5 conditioning PCs, q 0.05)
  combined by the pipeline's Overlap rule (2 of 2), sNMF K 5, imputed MAF 0.1. Each plot shows the region plus 100 kb
  on each side. Every series is a predictor and method with a significant SNP in that window.</p>
  <div class=legend>
    <span class=l-emmax><i class=dash></i>EMMAX Bonferroni, p {THR['EMMAX']}</span>
    <span class=l-rda><i class=dash></i>RDA q 0.05, p {THR['RDA']}</span>
    <span><i class=dot></i>region edges</span>
  </div>
  <p class=muted>Gene rows: <b>key</b> = key gene from the literature report (REGIONS_v2.md); <b>SNP host</b> = gene with a
  significant SNP in its body or 5 kb promoter; <b>outside</b> = a priority-A guess beyond the region edge, from the
  annotation only, with no literature check. Outside genes within 100 kb are labelled on the plot with a trailing *. SNP counts and best p
  are over all tested SNPs in the gene body / promoter. Most signals are carried by site 29072, the site with the widest
  temperature range and the coldest winter.</p>
</header>
<div class=tablewrap><table>
  <thead><tr><th>#</th><th>Region</th><th>kb</th><th>Genes</th><th>Sig. SNPs</th><th>Best p EMMAX / RDA</th>
  <th>Predictors</th><th>Candidate support</th></tr></thead>
  <tbody>{''.join(toc)}</tbody></table></div>
<section class=region id=manhattan>
  <header class=rhead><div class=rtitle><h2>Combined Manhattan, EMMAX + RDA</h2>
  <p>The pipeline's own plot (<span class=mono>GEA/plots/manhattan/combined/manhattan_combined_K5.png</span>), unchanged.
  Colour = predictor, shape = method, rectangles = the 15 regions. It draws one threshold line only, the strictest
  (EMMAX Bonferroni, p {THR['EMMAX']}); RDA hits pass at p {THR['RDA']}, so they sit below that line.</p></div></header>
  <figure><img src="img/manhattan_combined.png" alt="Combined Manhattan plot of EMMAX and RDA, all predictors, K 5" loading="lazy" width="4200" height="1500"></figure>
</section>
{''.join(sections)}
</div>
"""
# images are embedded as data URIs, so the page is one self-contained file (Artifact file upload was unavailable)
import base64
def _inline(m):
    with open(os.path.join(R, m.group(1)), "rb") as fh:
        return 'src="data:image/png;base64,' + base64.b64encode(fh.read()).decode() + '"'
page = re.sub(r'src="(img/[^"]+\.png)"', _inline, page)
assert len(page) < 15_000_000, len(page)
open(os.path.join(R, "regions_qc.html"), "w").write(page)
print(f"INFO: wrote report/regions_qc.html ({len(page)//1024} KB, images embedded), {len(sections)} regions, {len(genes)} gene rows")
