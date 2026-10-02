#!/usr/bin/env bash
# p4b_sweep.sh -- Phase 4b: garden sweep of all 600 SS-Clines replicates into offset13, then score.
#
# Prerequisites (Phase 4a + 4b prep, all done 2026-10-01): panels in MVP*_results/_intermediate/
# snp_sets (mvp_build_snp_sets.R), sweep configs with per-seed rosters (work/p4b_configs.sh, gate
# passed), fitness tables copied into offset13, all 600 dry runs clean (offset13/dryrun/summary.tsv),
# scripts/geometric_offset.R single-fit (equivalence: offset13/pilot/equiv/).
#
# Resources (user, 2026-10-01): up to 175 CPU / 800 GB. 40 seeds x snakemake -c4 = 160 cores;
# pilot peak 6.8 GiB for one seed at -c12, so 40 seeds stay far below the 800 GB host ceiling
# that mvp_garden_sweep.sh's watchdog enforces on TOTAL host use. Block b1 (11 panels) runs first.
#
#   nohup work/p4b_sweep.sh > benchmarks/mvp_eval/offset13/p4b_sweep.log 2>&1 &
set -uo pipefail
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1
O=benchmarks/mvp_eval/offset13
DK=(nix shell nixpkgs#docker-client -c docker)
UIDGID="$(id -u):$(id -g)"
SEED_JOBS="${SEED_JOBS:-40}"; SNAKE_CORES="${SNAKE_CORES:-4}"; MEM_PER_SEED="${MEM_PER_SEED:-24g}"
die() { echo "GATE FAIL: $*  $(date -Is)"; exit 1; }

# ---- preconditions
[[ $(wc -l < $O/dryrun/summary.tsv) -eq 600 ]] || die "dry-run summary is not 600 rows"
awk -F'\t' '$2!=0 || $3!="none"' $O/dryrun/summary.tsv | grep -q . && die "a dry run failed or wants upstream rules"
# b1 (11-panel rosters) first, then the rest
SEEDS=$(awk -F'\t' '{print gsub(/,/,",",$1)+1 "\t" $2}' work/p4b_rosters.tsv | sort -k1,1nr | cut -f2 | paste -sd,)
[[ $(tr ',' '\n' <<< "$SEEDS" | sort -u | wc -l) -eq 600 ]] || die "seed list is not 600 unique seeds"

# ---- sweep, with aggregate memory sampling
SP=$O/memtrace_sweep.tsv
( echo -e "time\tn_containers\tmem_gib\tcpu_pct\thost_used_gib"
  while true; do
    s=$("${DK[@]}" stats --no-stream --format '{{.Name}}\t{{.MemUsage}}\t{{.CPUPerc}}' 2>/dev/null | awk -F'\t' 'index($1,"gsweep_")==1 {split($2,a," "); m=a[1]; v=m+0; if (m ~ /GiB/) g=v; else if (m ~ /MiB/) g=v/1024; else if (m ~ /KiB/) g=v/1048576; else g=v/1073741824; n++; G+=g; C+=$3+0} END {printf "%d\t%.2f\t%.1f", n, G, C}')
    echo -e "$(date -Is)\t$s\t$(free -g | awk '/^Mem:/ {print $3}')"; sleep 30
  done ) > "$SP" 2>/dev/null & MS=$!
trap 'kill $MS 2>/dev/null' EXIT

echo "=== sweep start $(date -Is): 600 seeds, $SEED_JOBS x -c$SNAKE_CORES, mem/seed $MEM_PER_SEED"
OUT="$PWD/$O" LOGDIR="$PWD/$O/sweeplogs" CONFIG_SUFFIX=_sweep BLAS_THREADS=2 HOST_RAM_CEILING_GB=800 \
  nix shell nixpkgs#docker-client -c bash benchmarks/mvp_garden_sweep.sh "$SEEDS" "$SEED_JOBS" "$SNAKE_CORES" "$MEM_PER_SEED"
echo "=== sweep end $(date -Is)"
grep -c "FAILED" $O/sweeplogs/*.log 2>/dev/null | awk -F: '$2>0' | head

# ---- completeness gate: 112 gardens x 4 offset methods x roster size, per seed
short=0
while IFS=$'\t' read -r roster seeds n; do
  np=$(awk -F, '{print NF}' <<< "$roster"); want=$((112 * 4 * np))
  for s in ${seeds//,/ }; do
    got=$(find "$O/gardens/$s" -name '*.tsv' 2>/dev/null | wc -l)
    (( got == want )) || { echo "SHORT MVP$s: $got of $want"; short=$((short + 1)); }
  done
done < work/p4b_rosters.tsv
(( short == 0 )) || die "$short seed(s) short"
echo "completeness gate passed: $(find $O/gardens -name '*.tsv' | wc -l) harvested offset tables"

# ---- score (Lind & Lotterhos tau per garden / per source) and tables
"${DK[@]}" run --rm --name mvp-score-p4b --user "$UIDGID" -e USER=adaptogene -e PIPELINE_ROOT=/pipeline \
  -e OPENBLAS_NUM_THREADS=8 --cpus=16 --memory=96g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/eval_offset_lind.R --seeds="$SEEDS" --outdir=/pipeline/$O \
  > $O/logs_score.log 2>&1 || die "scoring rc=$?"
want=$(awk -F'\t' '{s += 112 * 4 * (gsub(/,/,",",$1)+1) * $3} END {print s}' work/p4b_rosters.tsv)
got=$(( $(wc -l < $O/garden_performance.tsv) - 1 ))
echo "scoring: garden_performance rows $got (harvested tables $want; scoring_skipped below)"
[[ -s $O/scoring_skipped.tsv ]] && { echo "scoring_skipped rows: $(( $(wc -l < $O/scoring_skipped.tsv) - 1 ))"; head -5 $O/scoring_skipped.tsv; }

"${DK[@]}" run --rm --name mvp-tables-p4b --user "$UIDGID" -e USER=adaptogene -e PIPELINE_ROOT=/pipeline \
  -e MVP_PANELS=offset13 -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_panel_tables.R --seeds="$SEEDS" --outdir=/pipeline/$O \
  > $O/logs_tables.log 2>&1 || die "tables rc=$?"
tail -3 $O/logs_tables.log
echo "=== Phase 4b sweep + scoring done $(date -Is)"
