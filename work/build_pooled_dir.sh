#!/usr/bin/env bash
# Stage I -- build benchmarks/mvp_eval/offset12_ssclines_pooled/ from the five per-block
# snapshots, so the 600-replicate report reads one directory.
#
# THE TABLES ARE NOT ALL THE SAME KIND. Verified on the b1/b2 snapshots:
#
#   garden_performance.tsv        120 / 120 seeds   PER-BLOCK  -> row-bind the five
#   source_performance.tsv        120 / 120         PER-BLOCK  -> row-bind
#   phase1_seed_medians_solo.tsv  120 / 120         PER-BLOCK  -> row-bind
#   panel_pr_recomputed.tsv       332 / 452         CUMULATIVE -> take the LAST block's copy
#   panel_offset_exclusions.tsv   160 / 167         CUMULATIVE -> take the LAST block's copy
#   panel_coverage.tsv            (no seed column)             -> take the LAST block's copy
#
# mvp_panel_tables.R rebuilds the cumulative ones for EVERY seed it can find, so b5's copy
# already contains all 692 manifest seeds. Row-binding those five times would duplicate every
# legacy row and silently inflate every denominator computed from them.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

EVAL=benchmarks/mvp_eval
OUT=$EVAL/offset12_ssclines_pooled
LAST=${LAST:-ssclines_b5}
BLOCKS=(ssclines_b1 ssclines_b2 ssclines_b3 ssclines_b4 ssclines_b5)

for b in "${BLOCKS[@]}"; do
    [[ -f "$EVAL/offset12_$b/garden_performance.tsv" ]] || { echo "FATAL: $b not snapshotted"; exit 1; }
done
mkdir -p "$OUT"

OUT=$OUT EVAL=$EVAL LAST=$LAST python3 - <<'PY' || exit 1
import csv, os, sys
out, ev, last = os.environ['OUT'], os.environ['EVAL'], os.environ['LAST']
blocks = ['ssclines_b1','ssclines_b2','ssclines_b3','ssclines_b4','ssclines_b5']

PER_BLOCK = ['garden_performance.tsv','source_performance.tsv','phase1_seed_medians_solo.tsv']
CUMULATIVE = ['panel_pr_recomputed.tsv','panel_offset_exclusions.tsv','panel_coverage.tsv']

for t in PER_BLOCK:
    header, rows, seeds = None, [], set()
    for b in blocks:
        p = f"{ev}/offset12_{b}/{t}"
        with open(p) as f:
            r = csv.reader(f, delimiter='\t')
            h = next(r)
            if header is None: header = h
            elif h != header:
                print(f"FATAL: {t} header differs in {b}"); sys.exit(1)
            si = h.index('seed') if 'seed' in h else None
            for row in r:
                rows.append(row)
                if si is not None: seeds.add(row[si])
    with open(f"{out}/{t}", 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(header); w.writerows(rows)
    print(f"  bound   {t:32s} {len(rows):7d} rows, {len(seeds) or '-'} distinct seeds")
    if seeds and len(seeds) != 600:
        print(f"  FATAL: {t} has {len(seeds)} distinct seeds, expected 600"); sys.exit(1)

import shutil
for t in CUMULATIVE:
    src = f"{ev}/offset12_{last}/{t}"
    if not os.path.exists(src):
        print(f"  (absent in {last}: {t})"); continue
    shutil.copy2(src, f"{out}/{t}")
    with open(src) as f:
        r = csv.reader(f, delimiter='\t'); h = next(r)
        n = sum(1 for _ in r)
    print(f"  copied  {t:32s} {n:7d} rows (cumulative, from {last})")
print("pooled dir built:", out)
PY

echo "=== contents ==="; ls -1 "$OUT"
