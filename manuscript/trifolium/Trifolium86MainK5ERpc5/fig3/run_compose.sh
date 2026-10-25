#!/usr/bin/env bash
# Compose the Trifolium main figures and render their PNGs at 600 dpi (named in compose_figs.py's docstring; written
# 2026-10-18, plan let-s-implement-this-plan-mighty-bubble.md). Host, from anywhere:
#   bash ANALYSIS/Trifolium86MainK5ERpc5/fig3/run_compose.sh [name ...]   (default: Fig3_structure_association Fig4_regions)
# 1. fig3/compose_figs.py writes manuscript_fig3/figures/<name>.svg (all composed figures, FigS13 pages included);
# 2. rsvg-convert (nix librsvg; Arial resolves to the metric-compatible Liberation Sans) renders <name>.png at 600 dpi;
# 3. rsvg-convert writes no dpi tag: a pHYs chunk (23622 px/m = 600 dpi, unit metre) is inserted before the first
#    IDAT (stdlib struct + zlib.crc32), then every PNG is parsed back and its pHYs asserted.
# Fig 5 is not composed: the ms_save() PNG of fig5_offset_ms already carries the tag.
set -euo pipefail
cd /mnt/data/eugene/ADAPTOGENE
AN=ANALYSIS/Trifolium86MainK5ERpc5
FIG=$AN/manuscript_fig3/figures
NAMES=("$@"); [ ${#NAMES[@]} -eq 0 ] && NAMES=(Fig3_structure_association Fig4_regions)
python3 -I $AN/fig3/compose_figs.py $AN
for n in "${NAMES[@]}"; do
    nix shell nixpkgs#librsvg -c rsvg-convert --unlimited -d 600 -p 600 -b white -f png -o "$FIG/$n.png" "$FIG/$n.svg"
done
python3 -I - "${NAMES[@]/#/$FIG/}" <<'EOF'
import struct, sys, zlib
PPM = 23622                                    # 600 dpi in pixels per metre
def chunks(b):
    i = 8
    while i < len(b):
        n, = struct.unpack(">I", b[i:i + 4]); t = b[i + 4:i + 8]
        yield i, t, b[i + 8:i + 8 + n]; i += 12 + n
for stem in sys.argv[1:]:
    f = stem + ".png"; b = open(f, "rb").read()
    assert b[:8] == b"\x89PNG\r\n\x1a\n", f
    if not any(t == b"pHYs" for _, t, _ in chunks(b)):
        idat = next(i for i, t, _ in chunks(b) if t == b"IDAT")
        body = struct.pack(">IIB", PPM, PPM, 1)
        ch = struct.pack(">I", len(body)) + b"pHYs" + body + struct.pack(">I", zlib.crc32(b"pHYs" + body) & 0xFFFFFFFF)
        b = b[:idat] + ch + b[idat:]
        open(f, "wb").write(b)
    ihdr = next(d for _, t, d in chunks(b) if t == b"IHDR"); w, h = struct.unpack(">II", ihdr[:8])
    x, y, u = struct.unpack(">IIB", next(d for _, t, d in chunks(b) if t == b"pHYs"))
    assert (x, y, u) == (PPM, PPM, 1), (f, x, y, u)
    print(f"{f}: {w} x {h} px = {w / 600 * 25.4:.1f} x {h / 600 * 25.4:.1f} mm, pHYs {x}/{y} px/m ({x * 0.0254:.0f} dpi)")
EOF
