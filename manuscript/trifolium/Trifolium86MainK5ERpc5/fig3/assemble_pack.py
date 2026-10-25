#!/usr/bin/env python3
"""Assemble the merged Trifolium figure pack (plan idempotent-conjuring-noodle.md, Phase 5).

Numbering (user 2026-10-17, provisional layout): Fig 3 = PCA 1-2 + Miami, Fig 4 = rgt_01 + rgt_12, Fig 5 = offset;
supplement S9-S19 in order of first mention, one S number per figure, the region set merged into S13 (three pages).
Composed figures come from fig3/compose_figs.py (already in manuscript_fig3/figures/); every other item is copied
from its draft / offset_minou source under its final name (offset S15-S19 = the *_ms final-size renders, SUPP_MS=1). Writes manuscript_fig3/supplement_numbering.tsv (the
single source for captions and the zip), checks every SVG for fonts and raw site ids, and writes manuscript_fig3.zip.
Run: python3 -I fig3/assemble_pack.py <analysis dir>
"""
import csv
import re
import shutil
import sys
import zipfile
from pathlib import Path

AN = Path(sys.argv[1])
MF = AN / "manuscript_fig3"
DR = MF / "draft"
OP = AN / "fig3/offset_minou/paper_pack"
OS = AN / "fig3/offset_minou/site_maps"
FIG = MF / "figures"
COMPOSED = None  # marker: file already written by compose_figs.py
FONT = "font-family: Arial, Helvetica, sans-serif;"

# number, final stem, source stem (None = composed), item ids, first-mention anchor
ITEMS = [
    ("Fig. 3", "Fig3_structure_association", COMPOSED, "st_pcapie12 + as_miami_v2", "Results, structure and association"),
    ("Fig. 4", "Fig4_regions", COMPOSED, "rgt_01 + rgt_12", "Results, candidate regions"),
    ("Fig. 5", "Fig5_offset", OP / "fig_main_ms", "off_f1 (fig_main_ms)", "Results, genomic offset"),
    ("Fig. S9", "FigS9_pca34", DR / "st_pcapie34", "st_pcapie34", "Results, structure"),
    ("Fig. S10", "FigS10_ancestry_by_site", DR / "st_bars_split", "st_bars_split", "Results, structure"),
    ("Fig. S11", "FigS11_piemap_bio3", DR / "piemap/piemap_K5_bio_3", "st_piemap_bio3", "Results, structure"),
    ("Fig. S12", "FigS12_emmax_per_predictor", DR / "as_pred14_v2", "as_pred14_v2", "Results, association"),
    ("Fig. S13", "FigS13_regions_p1", COMPOSED, "rgt_02, rgt_06, rgt_07", "Results, candidate regions"),
    ("Fig. S13", "FigS13_regions_p2", COMPOSED, "rgt_10, rgt_11, rgt_14", "Results, candidate regions"),
    ("Fig. S13", "FigS13_regions_p3", COMPOSED, "rgt_15", "Results, candidate regions"),
    ("Fig. S14", "FigS14_lfmm_k5", DR / "lf_k5_manhattan", "lf_k5_manhattan", "Results, association (stripes, LFMM)"),
    ("Fig. S15", "FigS15_site_offsets_by_method", OS / "fig_null_compare_27_ms", "off_s1", "Results, genomic offset"),
    ("Fig. S16", "FigS16_random_snps", OP / "fig_s_random_ms", "off_s2", "Results, genomic offset"),
    ("Fig. S17", "FigS17_predictor_change_vs_rank", OP / "fig_s_delta_heat_ms", "off_s3", "Results, genomic offset"),
    ("Fig. S18", "FigS18_offset_maps", OP / "fig_s_offset_maps_ms", "off_s5", "Results, genomic offset"),
    ("Fig. S19", "FigS19_k_choice", OP / "fig_s_k_choice_ms", "off_s6a", "Methods, choice of K"),
]

for num, stem, src, _, _ in ITEMS:
    for ext in ("svg", "png"):
        dst = FIG / f"{stem}.{ext}"
        if src is not None:
            shutil.copy2(f"{src}.{ext}", dst)
            if ext == "svg":   # ms_save's rewrite (benchmarks/ms_style.R:77) for figures saved by plain ggsave
                dst.write_text(dst.read_text().replace('font-family: "Nimbus Sans";', FONT))
        assert dst.exists() and dst.stat().st_size > 0, f"{dst} missing"

with open(MF / "supplement_numbering.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["number", "file", "items", "source", "first_mention"])
    for num, stem, src, items, anchor in ITEMS:
        w.writerow([num, f"figures/{stem}.svg", items,
                    "fig3/compose_figs.py" if src is None else str(Path(f"{src}.svg").relative_to(AN)), anchor])

# checks: fonts per SVG (reported), raw site ids in text elements (must be none)
ids = [r["site_id"] for r in csv.DictReader(open(AN / "fig3/site_names.tsv"), delimiter="\t")]
id_re = re.compile(r"\b(" + "|".join(ids) + r")\b")
bad = []
for num, stem, *_ in ITEMS:
    s = (FIG / f"{stem}.svg").read_text()
    fams = sorted(set(re.findall(r"font-family:[^;]*;", s)))
    texts = " ".join(re.findall(r"<text[^>]*>([^<]*)</text>", s))
    hit = sorted(set(id_re.findall(texts)))
    if hit:
        bad.append((stem, hit))
    print(f"{num:9s} {stem:34s} fonts={fams} raw_ids={hit or '-'}")
    assert fams == [FONT], f"{stem}: fonts {fams}"
assert not bad, f"raw site ids left in figure text: {bad}"

# zip: figures, tables, legends + offset text, numbering; drafts and review stay out
zp = AN / "manuscript_fig3.zip"
with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(FIG.iterdir()):
        z.write(p, f"manuscript_fig3/figures/{p.name}")
    for p in sorted((MF / "tables").iterdir()):
        z.write(p, f"manuscript_fig3/tables/{p.name}")
    for name in ("fig3.md", "supplement_numbering.tsv", "pipeline_deviations.md"):
        z.write(MF / name, f"manuscript_fig3/{name}")
    z.write(AN / "manuscript_offset_v2/manuscript_offset.ai.md", "manuscript_fig3/manuscript_offset.ai.md")
print(f"wrote {zp.name}: {len(zipfile.ZipFile(zp).namelist())} files")
