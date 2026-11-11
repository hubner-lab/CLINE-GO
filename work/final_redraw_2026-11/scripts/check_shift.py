#!/usr/bin/env python3
"""check_shift.py <ref svg> <new svg> <shift mm> <out tsv> -- Jobs 3 and 4: panel a moved, not resized.
Pairs the clip rectangles (svglite clipPath <rect>, x/y/width/height in pt) and the <text> nodes of the two SVGs in
document order: every pair must have the same width, height and y, and an x that is either unchanged (panel b, the
page) or moved right by exactly <shift mm> (panel a and its letter). Reports the counts and the largest residual."""
import re, sys
ref, new, shift_mm, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
PT = 72 / 25.4; sh = shift_mm * PT
def items(p):
    s = open(p).read()
    rects = [tuple(map(float, m)) for m in re.findall(r"<clipPath[^>]*>\s*<rect x='([-\d.]+)' y='([-\d.]+)' width='([-\d.]+)' height='([-\d.]+)'", s)]
    texts = [(float(x), float(y), t) for x, y, t in re.findall(r"<text[^>]*?x='([-\d.]+)' y='([-\d.]+)'[^>]*>([^<]*)</text>", s)]
    return rects, texts
(r1, t1), (r2, t2) = items(ref), items(new)
assert len(r1) == len(r2) and len(t1) == len(t2), (len(r1), len(r2), len(t1), len(t2))
rows, worst, moved, kept = [], 0.0, 0, 0
for kind, A, B in (("clip", r1, r2), ("text", t1, t2)):
    for a, b in zip(A, B):
        dx = b[0] - a[0]
        same = [abs(b[i] - a[i]) for i in range(1, 4 if kind == "clip" else 2)]
        if kind == "text":
            assert a[2] == b[2], (a, b)
        res = min(abs(dx), abs(dx - sh))
        cls = "moved" if abs(dx - sh) < abs(dx) else "kept"
        moved += cls == "moved"; kept += cls == "kept"
        worst = max(worst, res, *same)
        rows.append((kind, cls, f"{a[0]:.2f}", f"{b[0]:.2f}", f"{dx:.2f}", f"{max(same):.3f}", a[2] if kind == "text" else ""))
with open(out, "w") as f:
    f.write("kind\tclass\tx_ref_pt\tx_new_pt\tdx_pt\tmax_other_diff_pt\ttext\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
print(f"{out}: {len(rows)} items ({len(r1)} clips, {len(t1)} texts); moved by {shift_mm} mm: {moved}, unchanged: {kept}; "
      f"largest residual {worst:.3f} pt")
