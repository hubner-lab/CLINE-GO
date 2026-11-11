#!/usr/bin/env bash
# compose_fig3.sh <repro|new> -- Fig3_structure_association from panels a, b (the published 150 kb draft panels, read
# from work/ld_window_lab/scratch_audit3b/fig4_an, whose composite equals the published b918f2) and panel c from
# scratch/<variant>/j5/as_miami_v2.svg; then the 600 dpi PNG with its pHYs tag, as run_compose.sh / run_compose_3b.sh.
set -euo pipefail
cd /mnt/data/eugene/ADAPTOGENE
V=${1:?repro|new}; W=work/final_redraw_2026-11
SRC=work/ld_window_lab/scratch_audit3b/fig4_an/manuscript_fig3/draft
AN=$W/scratch/fig3_an_$V; D=$AN/manuscript_fig3/draft
mkdir -p $D/piemap
cp -p $SRC/st_pcapie12.svg $D/; cp -p $SRC/piemap/piemap_K5_bio_3.svg $D/piemap/; cp -p $W/scratch/$V/j5/as_miami_v2.svg $D/
python3 -I $W/scripts/compose_fig3.py $AN
FIG=$AN/manuscript_fig3/figures/Fig3_structure_association
nix shell --offline nixpkgs#librsvg -c rsvg-convert --unlimited -d 600 -p 600 -b white -f png -o "$FIG.png" "$FIG.svg"
python3 -I - "$FIG.png" <<'PY'
import struct, sys, zlib
PPM = 23622
f = sys.argv[1]; b = open(f, "rb").read()
def chunks(b):
    i = 8
    while i < len(b):
        n, = struct.unpack(">I", b[i:i + 4]); t = b[i + 4:i + 8]; yield i, t, b[i + 8:i + 8 + n]; i += 12 + n
if not any(t == b"pHYs" for _, t, _ in chunks(b)):
    idat = next(i for i, t, _ in chunks(b) if t == b"IDAT"); body = struct.pack(">IIB", PPM, PPM, 1)
    b = b[:idat] + struct.pack(">I", len(body)) + b"pHYs" + body + struct.pack(">I", zlib.crc32(b"pHYs" + body) & 0xFFFFFFFF) + b[idat:]
    open(f, "wb").write(b)
w, h = struct.unpack(">II", next(d for _, t, d in chunks(b) if t == b"IHDR")[:8])
x, y, u = struct.unpack(">IIB", next(d for _, t, d in chunks(b) if t == b"pHYs")); assert (x, y, u) == (PPM, PPM, 1)
print(f"{f}: {w} x {h} px, pHYs {x} px/m")
PY
