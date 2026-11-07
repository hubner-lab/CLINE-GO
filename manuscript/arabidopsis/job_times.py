#!/usr/bin/env python3
"""job_times.py: wall time of every job in one Snakemake run log.

    python3 -I manuscript/arabidopsis/job_times.py <run log> <out.tsv>

The pipeline records no per-job time (no benchmark: directives), so this reads the run log that
manuscript/arabidopsis/run_mode.sh keeps: each job opens with a "[timestamp]" line followed by
"localrule <rule>:" and its "jobid: N", and closes with "[timestamp]" then
"Finished jobid: N (Rule: <rule>)". Output columns: jobid, rule, wildcards, start, end, seconds;
a job that never finished has an empty end and seconds. Timestamps are the log's local time.
"""
import csv
import re
import sys
from datetime import datetime

STAMP = re.compile(r"^\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]$")
START = re.compile(r"^(?:local)?rule (\S+):$")
JOBID = re.compile(r"^\s+jobid: (\d+)$")
WILDCARDS = re.compile(r"^\s+wildcards: (.*)$")
FINISH = re.compile(r"^Finished jobid: (\d+) \(Rule: (\S+)\)$")


def parse(path):
    jobs = {}
    stamp = None
    current = None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            m = STAMP.match(line)
            if m:
                stamp = datetime.strptime(" ".join(m.group(1).split()), "%a %b %d %H:%M:%S %Y")
                current = None
                continue
            m = START.match(line)
            if m:
                current = {"rule": m.group(1), "start": stamp, "wildcards": "", "end": None}
                continue
            if current is not None:
                m = JOBID.match(line)
                if m:
                    jobs[int(m.group(1))] = current
                    continue
                m = WILDCARDS.match(line)
                if m:
                    current["wildcards"] = m.group(1)
                    continue
            m = FINISH.match(line)
            if m:
                job = jobs.get(int(m.group(1)))
                if job is None or job["rule"] != m.group(2):
                    sys.exit(f"job_times.py: finish of job {m.group(1)} ({m.group(2)}) has no matching start")
                job["end"] = stamp
    return jobs


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    jobs = parse(sys.argv[1])
    if not jobs:
        sys.exit(f"job_times.py: no job found in {sys.argv[1]}")
    with open(sys.argv[2], "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["jobid", "rule", "wildcards", "start", "end", "seconds"])
        for jobid in sorted(jobs):
            job = jobs[jobid]
            end = job["end"]
            seconds = int((end - job["start"]).total_seconds()) if end and job["start"] else ""
            writer.writerow([jobid, job["rule"], job["wildcards"],
                             job["start"].isoformat() if job["start"] else "",
                             end.isoformat() if end else "", seconds])


if __name__ == "__main__":
    main()
