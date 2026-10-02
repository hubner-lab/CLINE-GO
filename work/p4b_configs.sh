#!/usr/bin/env bash
# p4b_configs.sh -- Phase 4b step 8: one sweep config per SS-Clines replicate, per-seed roster.
#
# The roster comes from offset13/snp_sets_summary.tsv (mvp_build_snp_sets.R), never from a
# fixed list: a panel is swept only if usable_for_offset is TRUE (n >= 3 -- scripts/rda_offset.R
# FATALs below that and kills the whole seed's run). Swept panels: the 6 GEA candidates + truth,
# plus the size-curve panels where they exist (one block). all / neutral_all are never swept.
# Seeds sharing one roster go through ONE mvp_write_sweep_config.R call.
#
#   work/p4b_configs.sh            # all 600
#   work/p4b_configs.sh 1231118    # one seed (pilot)
set -euo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT"
SUM=benchmarks/mvp_eval/offset13/snp_sets_summary.tsv
ONLY="${1:-}"
DK=(nix shell nixpkgs#docker-client -c docker)

SUM="$SUM" ONLY="$ONLY" python3 - > work/p4b_rosters.tsv <<'PY'
import csv, os, collections
SWEPT = ["truth", "union", "best", "intersect3", "solo_lfmm", "solo_rda", "solo_emmax",
         "union_top0.1pct", "union_top0.5pct", "union_top1pct", "union_top2pct"]
rows = list(csv.DictReader(open(os.environ["SUM"]), delimiter="\t"))
only = os.environ["ONLY"]
by_seed = collections.defaultdict(list)
for r in rows:
    if r["set"] in SWEPT and r["usable_for_offset"] == "TRUE":
        by_seed[r["seed"]].append(r["set"])
seeds = sorted({r["seed"] for r in rows})
assert len(seeds) == 600, len(seeds)
if only: seeds = [only]
groups = collections.defaultdict(list)
for s in seeds:
    roster = [p for p in SWEPT if p in by_seed[s]]
    assert {"truth", "union", "best", "solo_lfmm", "solo_rda", "solo_emmax"} <= set(roster), (s, roster)
    groups[",".join(roster)].append(s)
for roster, ss in groups.items():
    print(f"{roster}\t{','.join(ss)}\t{len(ss)}")
PY
cut -f1,3 work/p4b_rosters.tsv
while IFS=$'\t' read -r roster seeds n; do
  "${DK[@]}" run --rm --name "mvp-gcfg-p4b-$RANDOM" --user "$(id -u):$(id -g)" -e USER=adaptogene \
    -e PIPELINE_ROOT=/pipeline -v "$PWD":/pipeline cline-go:latest \
    Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds="$seeds" --sets="$roster" | tail -1
done < work/p4b_rosters.tsv

# gate: every config lists exactly its roster, never all / neutral_all / rand_* / solo
bad=0
while IFS=$'\t' read -r roster seeds n; do
  for s in ${seeds//,/ }; do
    got=$(awk '/^  snp_sets:/{f=1;next} f && /^    - /{print $2; next} f{exit}' "config_MVP${s}_sweep.yaml" | paste -sd,)
    [[ "$got" == "$roster" ]] || { echo "ROSTER MISMATCH MVP$s: $got vs $roster"; bad=1; }
  done
done < work/p4b_rosters.tsv
(( bad == 0 )) && echo "config gate passed: $(awk -F'\t' '{n+=$3} END {print n}' work/p4b_rosters.tsv) configs"
exit $bad
