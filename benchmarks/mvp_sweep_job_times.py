#!/usr/bin/env python3
"""Per-job wall time of a garden sweep, parsed from its Snakemake run logs.

Why: pricing a garden sweep needs to know WHICH rule the hours go to, and Snakemake
writes no benchmark table for mode=maladaptation. The per-seed sweep logs
(mvp_garden_sweep.sh LOGDIR/MVP{seed}.log) do carry a timestamp before every
`rule X:` block and before every `Finished jobid: N`, so the duration of each job is
recoverable post hoc. Durations are WALL time under whatever concurrency that block
ran at (contention included), not CPU time.

Usage:
  python3 benchmarks/mvp_sweep_job_times.py <offset_dir> <out.tsv>
    reads <offset_dir>/sweeplogs_*/MVP*.log, writes one row per finished job:
    block  seed  rule  panel  sec
"""
import glob
import os
import re
import sys
from datetime import datetime

if len(sys.argv) != 3:
    sys.exit(__doc__)
root, out = sys.argv[1], sys.argv[2]

TS = re.compile(r'^\[(\w{3} \w{3} +\d+ [\d:]+ \d{4})\]')
RULE = re.compile(r'^(?:local)?rule (\w+):')
JOB = re.compile(r'\s+jobid: (\d+)')
WC = re.compile(r'\s+wildcards: (.*)')
FIN = re.compile(r'^Finished jobid: (\d+)')

rows = []
dirs = sorted(glob.glob(os.path.join(root, 'sweeplogs_*')))
if not dirs:
    sys.exit('no sweeplogs_* directories under ' + root)
for d in dirs:
    block = os.path.basename(d).replace('sweeplogs_', '')
    for f in sorted(glob.glob(os.path.join(d, 'MVP*.log'))):
        seed = os.path.basename(f)[3:-4]
        start, ts, cur = {}, None, None
        for line in open(f, errors='replace'):
            m = TS.match(line)
            if m:
                ts = datetime.strptime(re.sub(' +', ' ', m.group(1)), '%a %b %d %H:%M:%S %Y')
                continue
            m = RULE.match(line)
            if m:
                cur = {'rule': m.group(1), 't0': ts, 'panel': ''}
                continue
            if cur is not None:
                m = JOB.match(line)
                if m:
                    start[m.group(1)] = cur
                    continue
                m = WC.match(line)
                if m:
                    w = dict(x.split('=', 1) for x in m.group(1).split(', ') if '=' in x)
                    cur['panel'] = w.get('run_label', '')
                    continue
            m = FIN.match(line)
            if m and m.group(1) in start:
                j = start[m.group(1)]
                rows.append((block, seed, j['rule'], j['panel'], int((ts - j['t0']).total_seconds())))

with open(out, 'w') as o:
    o.write('block\tseed\trule\tpanel\tsec\n')
    for r in rows:
        o.write('\t'.join(map(str, r)) + '\n')
print('%d jobs from %d log dirs -> %s' % (len(rows), len(dirs), out))
