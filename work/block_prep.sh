#!/usr/bin/env bash
# Steps 6, 7, 8 for one SS-Clines block. Generalised from work/b4_prep.sh.
#
#   work/block_prep.sh <cohort> <cohort_whose_garden_qc_must_be_snapshotted_first>
#   work/block_prep.sh ssclines_b5 ssclines_b4
#
# Step 7's mvp_garden_fitness.R writes the fixed-name offset12/garden_qc.tsv, so it must not
# start until the previous block's copy has been snapshotted out. The wait is on that file,
# not a timer.
#
# Step 6's n>=3 pre-screen is not optional: a panel with <3 SNPs FATALs in
# scripts/rda_offset.R:174 and aborts the WHOLE seed's garden run, not just that panel.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

COHORT="${1:?cohort}"
PREV="${2:?previous cohort}"
SHORT="${COHORT#ssclines_}"
OFFSET_DIR=benchmarks/mvp_eval/offset12
DK=(nix shell nixpkgs#docker-client -c docker)
UIDGID="$(id -u):$(id -g)"
SEEDS_FILE=work/seeds_${SHORT}.txt
SEEDS=$(paste -sd, "$SEEDS_FILE")
NSEED=$(wc -l < "$SEEDS_FILE")
ALLSETS=truth,union,best,intersect3,rand_best1,solo_lfmm,solo_rda,solo_emmax

die() { echo "GATE FAIL: $*  $(date -Is)"; exit 1; }

echo "=== $COHORT prep: waiting for $PREV snapshot (garden_qc.tsv)  $(date -Is)"
for i in $(seq 1 2880); do
    [[ -f benchmarks/mvp_eval/offset12_${PREV}/garden_qc.tsv ]] && break
    sleep 30
done
[[ -f benchmarks/mvp_eval/offset12_${PREV}/garden_qc.tsv ]] || die "$PREV snapshot never appeared"
echo "=== $PREV snapshot present, starting $COHORT prep  $(date -Is)"

# ---------------------------------------------------------------- step 6
echo "--- step 6 build_snp_sets  $(date -Is)"
"${DK[@]}" run --rm --name mvp-snpsets-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_build_snp_sets.R --seeds="$SEEDS" \
  > logs_snpsets_${COHORT}.log 2>&1 || die "step 6 rc=$?"
tail -2 logs_snpsets_${COHORT}.log

echo "--- step 6 pre-screen: panels with n < 3 SNPs"
: > panel_underfilled_${COHORT}.tsv
for s in $(cat "$SEEDS_FILE"); do
  for p in truth union best intersect3 rand_best1 solo_lfmm solo_rda solo_emmax; do
    f="MVP${s}_results/_intermediate/snp_sets/${p}/selected_snps.tsv"
    if [[ -f "$f" ]]; then n=$(( $(wc -l < "$f") - 1 )); else n=0; fi
    (( n < 3 )) && printf 'MVP%s\t%s\t%s\n' "$s" "$p" "$n" >> panel_underfilled_${COHORT}.tsv
  done
done
echo "underfilled rows: $(wc -l < panel_underfilled_${COHORT}.tsv)  seeds affected: $(cut -f1 panel_underfilled_${COHORT}.tsv | sort -u | wc -l)"

# ---------------------------------------------------------------- step 7
# ORDER MATTERS AND IS NOT THE RUNBOOK'S. mvp_garden_fitness.R:59 REQUIRES a fitness-identity
# row for every seed it builds, and dies `no fitness-identity row` without one (:84). Its
# --identity default is benchmarks/mvp_eval/offset08/fitness_identity.tsv -- a STALE 32-row
# August file holding none of the SS-Clines seeds. The live table is the ACCUMULATED
# offset12/fitness_identity.tsv: 120 rows after b1, 240 after b2, 360 after b3, gates all TRUE.
# So: identity FIRST, merge into the accumulated table, THEN garden_fitness with --identity.
echo "--- step 7a fitness identity QC gate  $(date -Is)"
mkdir -p work/identity_${SHORT}
"${DK[@]}" run --rm --name mvp-identity-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=32g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_fitness_identity.R --seeds="$SEEDS" \
    --outdir=/pipeline/work/identity_${SHORT} > logs_identity_${COHORT}.log 2>&1 || die "step 7a rc=$?"
tail -2 logs_identity_${COHORT}.log
grep -q "All ${NSEED} seeds pass the identity gate" logs_identity_${COHORT}.log \
  || die "identity gate did not pass for all $NSEED seeds -- see logs_identity_${COHORT}.log"

echo "--- step 7b merge identity into the accumulated table  $(date -Is)"
IDENT=$OFFSET_DIR/fitness_identity.tsv
cp -p "$IDENT" "${IDENT}.pre_${COHORT}"
SHORT=$SHORT IDENT=$IDENT python3 - <<'PY' || exit 1
import csv, os, sys
ident = os.environ['IDENT']; short = os.environ['SHORT']
new_f = 'work/identity_%s/fitness_identity.tsv' % short
with open(ident) as f: old = list(csv.DictReader(f, delimiter='\t'))
with open(new_f) as f: new = list(csv.DictReader(f, delimiter='\t'))
have = {r['seed'] for r in old}
add  = [r for r in new if r['seed'] not in have]
dup  = [r['seed'] for r in new if r['seed'] in have]
bad  = [r['seed'] for r in new if r['gate'] != 'TRUE']
if bad:
    print("  FATAL: %d seed(s) with gate != TRUE: %s" % (len(bad), bad[:5])); sys.exit(1)
if dup:
    print("  NOTE: %d seed(s) already present, not re-added: %s" % (len(dup), dup[:5]))
rows = old + add
with open(ident, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(old[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print("  accumulated identity table: %d -> %d rows (+%d)" % (len(old), len(rows), len(add)))
PY
echo "--- step 7c garden_fitness  $(date -Is)"
"${DK[@]}" run --rm --name mvp-fitness-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 --cpus=16 --memory=96g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_garden_fitness.R --seeds="$SEEDS" \
    --identity=/pipeline/$OFFSET_DIR/fitness_identity.tsv \
    --outdir=/pipeline/$OFFSET_DIR > logs_fitness_${COHORT}.log 2>&1 || die "step 7c rc=$?"
tail -2 logs_fitness_${COHORT}.log

echo "--- step 7c garden env tables (12 concurrent)  $(date -Is)"
# --gardens is MANDATORY (defaults to offset09); --outdir must NOT be passed.
n=0
for s in $(cat "$SEEDS_FILE"); do
  (
    "${DK[@]}" run --rm --name mvp-genv-$s --user "$UIDGID" -e USER=adaptogene \
      -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g \
      -v "$PWD":/pipeline cline-go:latest \
      Rscript /pipeline/benchmarks/mvp_write_garden_env.R --seed=$s \
        --gardens=/pipeline/$OFFSET_DIR/gardens_${s}.tsv \
      > work/genv_${SHORT}_${s}.log 2>&1 || echo "FAILED genv $s"
  ) &
  n=$((n+1)); (( n % 12 == 0 )) && wait
done
wait
echo "genv done  $(date -Is)"

echo "--- step 7 gate"
short=0
for s in $(cat "$SEEDS_FILE"); do
  [[ -s "$OFFSET_DIR/gardens_${s}.tsv" ]]        || { echo "  missing gardens_${s}.tsv"; short=1; }
  [[ -s "$OFFSET_DIR/garden_fitness_${s}.tsv" ]] || { echo "  missing garden_fitness_${s}.tsv"; short=1; }
  c=$(ls data/mvp/MVP${s}/gardens/ 2>/dev/null | wc -l)
  (( c != 112 )) && { echo "  MVP${s} env tables=$c (want 112)"; short=1; }
done
(( short == 1 )) && die "step 7 gate"
echo "step 7 gate passed: $NSEED seeds x 112 env tables"

# ---------------------------------------------------------------- step 8
echo "--- step 8 garden sweep configs  $(date -Is)"
"${DK[@]}" run --rm --name mvp-gcfg-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds="$SEEDS" \
    --sets=$ALLSETS > logs_gcfg_${COHORT}.log 2>&1 || die "step 8 rc=$?"
tail -1 logs_gcfg_${COHORT}.log

echo "--- step 8 per-seed rewrite for underfilled seeds"
for s in $(cut -f1 panel_underfilled_${COHORT}.tsv | sort -u | sed 's/^MVP//'); do
  keep=$(COHORT=$COHORT SEED=$s python3 -c "
import os
ALL='truth,union,best,intersect3,rand_best1,solo_lfmm,solo_rda,solo_emmax'.split(',')
c=os.environ['COHORT']; s='MVP'+os.environ['SEED']
drop={l.split(chr(9))[1] for l in open('panel_underfilled_%s.tsv'%c) if l.split(chr(9))[0]==s}
print(','.join(p for p in ALL if p not in drop))")
  echo "  MVP$s -> $keep"
  "${DK[@]}" run --rm --name mvp-gcfgfix-$s --user "$UIDGID" -e USER=adaptogene \
    -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g \
    -v "$PWD":/pipeline cline-go:latest \
    Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds=$s --sets="$keep" \
    > work/gcfg_fix_${SHORT}_${s}.log 2>&1 || echo "FAILED gcfg fix $s"
done

echo "--- step 8 gate: no neutral_all / bare 'all' in $COHORT sweep configs"
bad=0
for s in $(cat "$SEEDS_FILE"); do
  grep -qE '^\s*-\s*(all|neutral_all)\s*$' "config_MVP${s}_sweep.yaml" 2>/dev/null \
    && { echo "  FAIL config_MVP${s}_sweep.yaml contains all/neutral_all"; bad=1; }
done
(( bad == 1 )) && die "step 8 gate"
one=$(head -1 "$SEEDS_FILE"); grep -A2 emit_plots "config_MVP${one}_sweep.yaml" | head -3
echo "=== $COHORT prep done  $(date -Is)"
