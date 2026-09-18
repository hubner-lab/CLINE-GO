#!/usr/bin/env bash
# Regenerate b4's garden sweep configs after the sensitivity_check fix in
# benchmarks/mvp_write_sweep_config.R. Step 8 only -- steps 6 and 7 already passed.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1
COHORT=ssclines_b4
DK=(nix shell nixpkgs#docker-client -c docker)
UIDGID="$(id -u):$(id -g)"
SEEDS=$(paste -sd, work/seeds_b4.txt)
ALLSETS=truth,union,best,intersect3,rand_best1,solo_lfmm,solo_rda,solo_emmax

echo "--- step 8 regenerate all 120  $(date -Is)"
"${DK[@]}" run --rm --name mvp-gcfg-${COHORT}-redo --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds="$SEEDS" \
    --sets=$ALLSETS > logs_gcfg_${COHORT}.log 2>&1 || { echo "GATE FAIL: step 8"; exit 1; }
tail -1 logs_gcfg_${COHORT}.log

echo "--- step 8 per-seed rewrite for the 35 underfilled seeds  $(date -Is)"
for s in $(cut -f1 panel_underfilled_${COHORT}.tsv | sort -u | sed 's/^MVP//'); do
  keep=$(COHORT=$COHORT SEED=$s python3 -c "
import os
ALL='truth,union,best,intersect3,rand_best1,solo_lfmm,solo_rda,solo_emmax'.split(',')
c=os.environ['COHORT']; s='MVP'+os.environ['SEED']
drop={l.split(chr(9))[1] for l in open('panel_underfilled_%s.tsv'%c) if l.split(chr(9))[0]==s}
print(','.join(p for p in ALL if p not in drop))")
  "${DK[@]}" run --rm --name mvp-gcfgfix-$s --user "$UIDGID" -e USER=adaptogene \
    -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g \
    -v "$PWD":/pipeline cline-go:latest \
    Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds=$s --sets="$keep" \
    > work/gcfg_fix_b4_${s}.log 2>&1 || echo "FAILED gcfg fix $s"
done

echo "--- gate: sensitivity_check present in all 120, no all/neutral_all"
miss=0; bad=0
for s in $(cat work/seeds_b4.txt); do
  grep -q "sensitivity_check: false" "config_MVP${s}_sweep.yaml" || { echo "  MISSING in MVP${s}"; miss=1; }
  grep -qE '^\s*-\s*(all|neutral_all)\s*$' "config_MVP${s}_sweep.yaml" && { echo "  all/neutral_all in MVP${s}"; bad=1; }
done
(( miss == 1 || bad == 1 )) && { echo "GATE FAIL: step 8 redo"; exit 1; }
echo "=== ssclines_b4 prep done (step 8 regenerated)  $(date -Is)"
