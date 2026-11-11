#!/usr/bin/env python3
"""[final redraw 2026-11-11] Copy of manuscript/trifolium/Trifolium86MainK5ERpc5/fig3/compose_figs.py cut to Fig 3
(the Fig 4 and S13 blocks and their REG lookup removed; panel(), letter(), compose() and the Fig 3 call unchanged).
Run: python3 -I compose_fig3.py <analysis dir with manuscript_fig3/draft/{st_pcapie12.svg,piemap/piemap_K5_bio_3.svg,as_miami_v2.svg}>

Compose the manuscript figures from panel SVGs drawn at final size (plan idempotent-conjuring-noodle.md, Phase 5).

Layout (user 2026-10-17, provisional): Fig 3 = a PCA 1-2 + b Miami; Fig 4 = a rgt_01 + b rgt_12; Fig 5 = offset
Fig 1 redrawn at 180 mm (fig_main_ms, letters a-c already drawn by patchwork); S13 = the supplementary region set,
three regions per page, letters a-g running across pages.

Panels are placed unscaled (every one was saved at its final mm size by ms_save) as nested <svg> elements; svglite
clip-path ids are prefixed per panel so they stay unique in the composite. Panel letters: Arial bold 8 pt.
Writes manuscript_fig3/figures/<name>.svg; PNGs are rendered at 600 dpi by run_compose.sh (rsvg-convert).
2026-10-18 (plan let-s-implement-this-plan-mighty-bubble.md): Fig 3 = a PCA + b BIO3 piemap over c Miami; a panel may
be cropped to a vertical band (empty top/bottom margins trimmed, content unscaled); Fig 5 is no longer composed
(build_paper_pack.R MAIN6_MS=1 writes it at final size; copied into figures/ as it is).
Run: python3 -I fig3/compose_figs.py <analysis dir>
"""
import re
import sys
from pathlib import Path

PT = 72 / 25.4                      # pt per mm
AN = Path(sys.argv[1])
DR = AN / "manuscript_fig3/draft"
OM = AN / "fig3/offset_minou/paper_pack"
OUT = AN / "manuscript_fig3/figures"
OUT.mkdir(parents=True, exist_ok=True)
FONT = "font-family: Arial, Helvetica, sans-serif;"



def panel(path: Path, x_mm: float, y_mm: float, tag: str, crop: tuple[float, float] | None = None) -> tuple[str, float, float]:
    """Return (nested svg markup, width mm, height mm) of one panel at (x, y) mm, ids prefixed with tag.
    crop = (top mm, bottom mm): show only that band of the panel, unscaled (trims empty margins)."""
    s = path.read_text()
    s = re.sub(r"<\?xml[^>]*\?>\s*", "", s)
    m = re.match(r"\s*<svg[^>]*width='([\d.]+)pt' height='([\d.]+)pt' viewBox='([^']+)'[^>]*>", s)
    assert m, f"{path.name}: no svglite root"
    w, h, vb = float(m.group(1)), float(m.group(2)), m.group(3)
    body = s[m.end():s.rindex("</svg>")]
    ids = set(re.findall(r"id='([^']+)'", body))
    for i in ids:
        body = body.replace(f"id='{i}'", f"id='{tag}_{i}'").replace(f"url(#{i})", f"url(#{tag}_{i})")
    assert not re.search(r"href='#", body), f"{path.name}: internal href not handled"
    if crop:
        x0, y0, vw, vh = map(float, vb.split())
        assert (x0, y0, vw, vh) == (0.0, 0.0, w, h), f"{path.name}: viewBox not in pt"
        assert 0 <= crop[0] < crop[1] <= h / PT + 0.01, f"{path.name}: crop outside the panel"
        h = (crop[1] - crop[0]) * PT
        vb = f"0 {crop[0] * PT:.2f} {w:.2f} {h:.2f}"
    nested = (f"<svg x='{x_mm * PT:.2f}' y='{y_mm * PT:.2f}' width='{w:.2f}' height='{h:.2f}' "
              f"viewBox='{vb}' overflow='hidden'>{body}</svg>")
    return nested, w / PT, h / PT


def letter(ch: str, x_mm: float, y_mm: float) -> str:
    return (f"<text x='{(x_mm + 1.0) * PT:.2f}' y='{(y_mm + 3.2) * PT:.2f}' "
            f"style='font-size: 8.00px; font-weight: bold; fill: #222222; {FONT}'>{ch}</text>")


def compose(name: str, items: list[tuple], w_mm: float, h_mm: float) -> None:
    parts = []
    for k, (p, x, y, ch, *crop) in enumerate(items):
        nested, pw, ph = panel(p, x, y, f"p{k}", crop[0] if crop else None)
        assert x + pw <= w_mm + 0.05 and y + ph <= h_mm + 0.05, f"{name}: {p.name} overflows the page"
        parts.append(nested)
        if ch:
            parts.append(letter(ch, x, y))
    svg = ("<?xml version='1.0' encoding='UTF-8' ?>\n"
           "<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' "
           f"width='{w_mm * PT:.2f}pt' height='{h_mm * PT:.2f}pt' viewBox='0 0 {w_mm * PT:.2f} {h_mm * PT:.2f}'>\n"
           "<rect width='100%' height='100%' style='stroke: none; fill: #FFFFFF;'/>\n"
           + "\n".join(parts) + "\n</svg>\n")
    fams = set(re.findall(r"font-family:[^;]*;", svg))
    assert fams == {FONT}, f"{name}: fonts {fams}"
    assert len(re.findall(r"id='([^']+)'", svg)) == len(set(re.findall(r"id='([^']+)'", svg))), f"{name}: duplicate ids"
    (OUT / f"{name}.svg").write_text(svg)
    print(f"wrote {name}.svg  {w_mm:.0f} x {h_mm:.1f} mm, {len(items)} panels")


# Fig 3 (2026-10-18): a PCA 1-2 (80 x 112) | b BIO3 piemap (100 x 112, = FigS11 render), both cropped to the band
# 7.5-104.5 mm (their drawn content spans 9.2-102.7 and 13.2-98.8 mm), over c Miami (180 x 90)
ROW1 = (7.5, 104.5)
H1 = ROW1[1] - ROW1[0]
compose("Fig3_structure_association",
        [(DR / "st_pcapie12.svg", 0, 0, "a", ROW1), (DR / "piemap/piemap_K5_bio_3.svg", 80, 0, "b", ROW1),
         (DR / "as_miami_v2.svg", 0, H1, "c")], 180, H1 + 90)
