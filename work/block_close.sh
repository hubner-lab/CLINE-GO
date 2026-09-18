#!/usr/bin/env bash
# Steps 10 + 11 + snapshot for one SS-Clines block. Generalised from work/b3_close.sh.
#
#   work/block_close.sh <cohort> <prev_cohort_for_recheck>
#   work/block_close.sh ssclines_b4 ssclines_b3
#
# THE SERIALIZED TRIPLE. eval_offset_lind.R and mvp_panel_tables.R write fixed filenames into
# the shared offset12/ (garden_performance, source_performance, scoring_skipped, panel_*), so
# no other block's step 10/11 may run while this does. flock on work/.offset12_score.lock is
# what enforces that between the b3/b4/b5 chains.
#
# The snapshot copies MORE than the runbook's six files: garden_qc.tsv and scoring_skipped.tsv
# are the same class of fixed-name overwrite but are absent from the runbook's list.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

COHORT="${1:?cohort, e.g. ssclines_b4}"
PREV="${2:?previous cohort for the medians re-check, e.g. ssclines_b3}"
SHORT="${COHORT#ssclines_}"                      # b4
OFFSET_DIR=benchmarks/mvp_eval/offset12
SNAP=benchmarks/mvp_eval/offset12_${COHORT}
DK=(nix shell nixpkgs#docker-client -c docker)
UIDGID="$(id -u):$(id -g)"
SEEDS_FILE=work/seeds_${SHORT}.txt
SEEDS=$(paste -sd, "$SEEDS_FILE")

exec 9>work/.offset12_score.lock
flock 9 || { echo "FATAL: cannot take offset12 scoring lock"; exit 1; }
echo "=== $COHORT close-out, scoring lock held  $(date -Is)"

# ---------------------------------------------------------------- step 10
echo "--- step 10 score  $(date -Is)"
"${DK[@]}" run --rm --name mvp-score-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=8 --cpus=16 --memory=96g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/eval_offset_lind.R --seeds="$SEEDS" \
    --outdir=/pipeline/$OFFSET_DIR > logs_score_${COHORT}.log 2>&1
rc=$?
tail -3 logs_score_${COHORT}.log
[[ $rc -ne 0 ]] && { echo "GATE FAIL: step 10 rc=$rc  $(date -Is)"; exit 1; }

EXPECT=$(COHORT=$COHORT SHORT=$SHORT python3 - <<'PY'
import collections, os
c = os.environ['COHORT']; s = os.environ['SHORT']
drop = collections.Counter()
try:
    for line in open('panel_underfilled_%s.tsv' % c):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2: drop[f[0]] += 1
except FileNotFoundError:
    pass
seeds = [l.strip() for l in open('work/seeds_%s.txt' % s) if l.strip()]
print(112 * 4 * sum(8 - drop.get('MVP'+x, 0) for x in seeds))
PY
)
GOT=$(( $(wc -l < $OFFSET_DIR/garden_performance.tsv) - 1 ))
echo "step 10 gate: garden_performance rows got=$GOT expect=$EXPECT"
[[ "$GOT" == "$EXPECT" ]] || { echo "GATE FAIL: row count mismatch  $(date -Is)"; exit 1; }
if [[ -s "$OFFSET_DIR/scoring_skipped.tsv" ]]; then
    echo "--- scoring_skipped.tsv (cross-check against panel_underfilled_${COHORT}.tsv):"
    head -20 "$OFFSET_DIR/scoring_skipped.tsv"
fi

# ---------------------------------------------------------------- step 11
echo "--- step 11 tables, --check=offset11  $(date -Is)"
"${DK[@]}" run --rm --name mvp-tables-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_panel_tables.R \
    --outdir=/pipeline/$OFFSET_DIR --check=/pipeline/benchmarks/mvp_eval/offset11 \
  > logs_tables_${COHORT}.log 2>&1
echo "step 11 --check=offset11 exit=$?  (2 = NOT APPLICABLE, expected, NOT a pass)"
tail -5 logs_tables_${COHORT}.log

# --- direction-aware re-check against the previous block's snapshot.
#
# `mvp_panel_tables.R --check` CANNOT be used for this. panel_pr_recomputed.tsv is CUMULATIVE
# -- rebuilt for every seed with a truth table and snp_sets on disk -- so a later block's
# panels appear in it as `n = 0` rows (the empty-set branch) BEFORE that block's step 6 runs,
# and as real values after. Run on b3 it reported `max |diff| = 3.16e+04 ... GATE FAILED`
# purely from 1427 rows going 0 -> populated. Every future block would fail it the same way.
#
# What actually matters is the DIRECTION of each change:
#   empty -> populated  = this block's panels being built. Benign.
#   populated -> changed / populated -> empty = a real regression. Stop.
echo "--- step 11 direction-aware re-check vs offset12_${PREV}  $(date -Is)"
PREV=$PREV OFFSET_DIR=$OFFSET_DIR python3 - <<'PY'
import csv, os, sys, collections
prev = 'benchmarks/mvp_eval/offset12_%s/panel_pr_recomputed.tsv' % os.environ['PREV']
new  = os.path.join(os.environ['OFFSET_DIR'], 'panel_pr_recomputed.tsv')
def load(p):
    d = {}
    try:
        with open(p) as f:
            for row in csv.DictReader(f, delimiter='\t'): d[(row['seed'], row['set'])] = row
    except FileNotFoundError:
        print("  (no %s -- re-check skipped)" % p); sys.exit(0)
    return d
A, B = load(prev), load(new)
same = grew = changed = lost = 0
bad = []
for k in set(A) & set(B):
    a, b = A[k], B[k]
    if a == b: same += 1
    elif a['n'] == '0' and b['n'] != '0': grew += 1
    elif a['n'] != '0' and b['n'] == '0': lost += 1; bad.append((k, a['n'], b['n']))
    else: changed += 1; bad.append((k, a['n'], b['n']))
print("  common keys=%d  identical=%d  empty->populated=%d  changed=%d  populated->empty=%d"
      % (len(set(A) & set(B)), same, grew, changed, lost))
for x in bad[:5]: print("   REGRESSION:", x)
print("  VERDICT:", "CLEAN (only new panels appeared)" if not bad else "REAL REGRESSION -- diagnose")
PY

# ---------------------------------------------------------------- snapshot
echo "--- snapshot -> $SNAP  $(date -Is)"
mkdir -p "$SNAP"
for f in garden_performance.tsv source_performance.tsv phase1_seed_medians_solo.tsv \
         panel_pr_recomputed.tsv panel_offset_exclusions.tsv panel_coverage.tsv \
         garden_qc.tsv scoring_skipped.tsv; do
    if [[ -f "$OFFSET_DIR/$f" ]]; then cp -p "$OFFSET_DIR/$f" "$SNAP/$f"; else echo "  (absent: $f)"; fi
done
cp -p "panel_underfilled_${COHORT}.tsv" "$SNAP/" 2>/dev/null
# Snapshot BOTH identity tables: this block's own rows, and the accumulated table that
# mvp_garden_fitness.R --identity actually reads (120 rows after b1, 240 after b2, ...).
# b1/b2 preserved only the accumulated one, as fitness_identity.b1only/.b1b2.tsv.
[[ -f work/identity_${SHORT}/fitness_identity.tsv ]] && \
    cp -p work/identity_${SHORT}/fitness_identity.tsv "$SNAP/fitness_identity.${SHORT}.tsv"
[[ -f $OFFSET_DIR/fitness_identity.tsv ]] && \
    cp -p $OFFSET_DIR/fitness_identity.tsv "$SNAP/fitness_identity.accumulated.tsv"

COHORT=$COHORT SHORT=$SHORT python3 - <<'PY'
import collections, os
c = os.environ['COHORT']; s = os.environ['SHORT']
ALL = ['truth','union','best','intersect3','rand_best1','solo_lfmm','solo_rda','solo_emmax']
drop = collections.defaultdict(set)
try:
    for line in open('panel_underfilled_%s.tsv' % c):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2: drop[f[0]].add(f[1])
except FileNotFoundError:
    pass
with open('sets_per_seed_%s.tsv' % c, 'w') as out:
    for seed in (l.strip() for l in open('work/seeds_%s.txt' % s) if l.strip()):
        d = drop.get('MVP'+seed, set())
        out.write("%s\t%s\t%s\n" % (seed, ','.join(p for p in ALL if p not in d),
                                    '|'.join(sorted(d)) or 'none'))
print("sets_per_seed_%s.tsv written" % c)
PY
cp -p sets_per_seed_${COHORT}.tsv "$SNAP/"

echo "=== $COHORT snapshot contents:"; ls -1 "$SNAP"
echo "=== $COHORT close-out done  $(date -Is)"
