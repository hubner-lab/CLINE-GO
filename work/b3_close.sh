#!/usr/bin/env bash
# Stage B -- close out b3 (ssclines_nequal_mconst): step 10 score, step 11 tables, snapshot.
#
# This is the SERIALIZED TRIPLE. eval_offset_lind.R and mvp_panel_tables.R both write fixed
# filenames into the shared offset12/ (garden_performance, source_performance,
# scoring_skipped, panel_*), so no other block's step 10/11 may run while this does.
# A lock file enforces that against the b4/b5 chains.
#
# The snapshot copies MORE than the runbook's six files: garden_qc.tsv and scoring_skipped.tsv
# are the same class of fixed-name overwrite but are absent from the runbook list, and b3's
# garden_qc.tsv is destroyed the moment b4's step 7 runs. Copy before Stage C starts.
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

COHORT=ssclines_b3
TAG=ssclines_nequal_mconst
OFFSET_DIR=benchmarks/mvp_eval/offset12
SNAP=benchmarks/mvp_eval/offset12_${COHORT}
LOCK=work/.offset12_score.lock
DK=(nix shell nixpkgs#docker-client -c docker)
UIDGID="$(id -u):$(id -g)"
SEEDS=$(paste -sd, work/seeds_b3.txt)

exec 9>"$LOCK"
flock 9 || { echo "FATAL: cannot take offset12 scoring lock"; exit 1; }
echo "=== b3 close-out, scoring lock held  $(date -Is)"

# ---------------------------------------------------------------- step 10: score
echo "--- step 10 score  $(date -Is)"
"${DK[@]}" run --rm --name mvp-score-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=8 --cpus=16 --memory=96g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/eval_offset_lind.R --seeds="$SEEDS" \
    --outdir=/pipeline/$OFFSET_DIR > logs_score_${COHORT}.log 2>&1
rc=$?
tail -3 logs_score_${COHORT}.log
[[ $rc -ne 0 ]] && { echo "GATE FAIL: step 10 rc=$rc"; exit 1; }

# Gate: row count must equal 112 gardens x 4 methods x (panels kept per seed), where the
# panels dropped at step 6 are exactly those in panel_underfilled_ssclines_b3.tsv.
EXPECT=$(python3 - <<'PY'
import collections
drop = collections.Counter()
for line in open('panel_underfilled_ssclines_b3.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 2: drop[f[0]] += 1
seeds = [l.strip() for l in open('work/seeds_b3.txt') if l.strip()]
panels = sum(8 - drop.get('MVP'+s, 0) for s in seeds)
print(112 * 4 * panels)
PY
)
GOT=$(( $(wc -l < $OFFSET_DIR/garden_performance.tsv) - 1 ))
echo "step 10 gate: garden_performance rows got=$GOT expect=$EXPECT"
[[ "$GOT" == "$EXPECT" ]] || { echo "GATE FAIL: row count mismatch"; exit 1; }
NA=$(awk -F'\t' 'NR>1 && /\tNA\t|\tNA$/ {n++} END {print n+0}' $OFFSET_DIR/garden_performance.tsv)
echo "step 10 gate: rows containing NA = $NA"

# scoring_skipped must name nothing outside panel_underfilled
if [[ -s "$OFFSET_DIR/scoring_skipped.tsv" ]]; then
    echo "--- scoring_skipped.tsv contents (cross-check against panel_underfilled):"
    head -20 "$OFFSET_DIR/scoring_skipped.tsv"
fi

# ---------------------------------------------------------------- step 11: tables
echo "--- step 11 tables, regression gate vs offset11  $(date -Is)"
"${DK[@]}" run --rm --name mvp-tables-${COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_panel_tables.R \
    --outdir=/pipeline/$OFFSET_DIR --check=/pipeline/benchmarks/mvp_eval/offset11 \
  > logs_tables_${COHORT}.log 2>&1
echo "step 11 --check=offset11 exit=$?  (2 = NOT APPLICABLE, expected, NOT a pass)"
tail -5 logs_tables_${COHORT}.log

# The medians table gets no regression check against offset11 (no shared seeds), so re-check
# against the previous block's own snapshot instead.
echo "--- step 11 re-check vs offset12_ssclines_b2  $(date -Is)"
"${DK[@]}" run --rm --name mvp-tables-${COHORT}-recheck --user "$UIDGID" -e USER=adaptogene \
  -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_panel_tables.R \
    --outdir=/pipeline/$OFFSET_DIR --check=/pipeline/benchmarks/mvp_eval/offset12_ssclines_b2 \
  > logs_tables_${COHORT}_recheck.log 2>&1
echo "step 11 --check=b2 exit=$?"
tail -5 logs_tables_${COHORT}_recheck.log

# ---------------------------------------------------------------- snapshot
echo "--- snapshot -> $SNAP  $(date -Is)"
mkdir -p "$SNAP"
for f in garden_performance.tsv source_performance.tsv phase1_seed_medians_solo.tsv \
         panel_pr_recomputed.tsv panel_offset_exclusions.tsv panel_coverage.tsv \
         garden_qc.tsv scoring_skipped.tsv; do
    if [[ -f "$OFFSET_DIR/$f" ]]; then cp -p "$OFFSET_DIR/$f" "$SNAP/$f"; else echo "  (absent: $f)"; fi
done
cp -p "panel_underfilled_${COHORT}.tsv" "$SNAP/" 2>/dev/null
[[ -f work/identity_b3/fitness_identity.tsv ]] && cp -p work/identity_b3/fitness_identity.tsv "$SNAP/fitness_identity.b3.tsv"

# sets_per_seed: seed <tab> sets actually used <tab> dropped|none  (b2's 3-column shape)
python3 - <<'PY'
import collections
ALL = ['truth','union','best','intersect3','rand_best1','solo_lfmm','solo_rda','solo_emmax']
drop = collections.defaultdict(set)
for line in open('panel_underfilled_ssclines_b3.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 2: drop[f[0]].add(f[1])
with open('sets_per_seed_ssclines_b3.tsv','w') as out:
    for s in (l.strip() for l in open('work/seeds_b3.txt') if l.strip()):
        d = drop.get('MVP'+s, set())
        kept = [p for p in ALL if p not in d]
        out.write("%s\t%s\t%s\n" % (s, ','.join(kept), '|'.join(sorted(d)) or 'none'))
print("sets_per_seed_ssclines_b3.tsv written")
PY
cp -p sets_per_seed_ssclines_b3.tsv "$SNAP/"

echo "=== b3 snapshot contents:"; ls -1 "$SNAP"
echo "=== b3 close-out done  $(date -Is)"
