#!/usr/bin/env bash
# =============================================================================
# sweep_rda.sh — RDA Condition()-PC ladder with per-SNP p-values retained.
#
# Why this exists: mode=pregea's RDA block persists only ladder DIAGNOSTICS
# (candidate counts, gif, R2adj, VIF) — not per-SNP p-values — so its rungs
# cannot be scored against the causal-locus truth set. This driver calls the
# pipeline's own scripts/rda.R once per rung, writing a full p-value table per
# condition_pcs value into its own directory. No pipeline code is modified.
#
# Every input path is resolved by glob at runtime, because the _work/ directory
# names embed the filter and LD parameters (maf*/ld*). A missing or ambiguous
# match aborts immediately rather than handing a bad path to rda.R.
#
# Usage:
#   benchmarks/sweep_rda.sh [PROJECT] [COND_PCS_CSV] [JOBS]
# Defaults: LARUSON1K, 0..10, 3 concurrent
#
# Run inside the cline-go container with /pipeline mounted.
# =============================================================================
set -euo pipefail

PROJECT="${1:-LARUSON1K}"
COND_PCS_CSV="${2:-0,1,2,3,4,5,6,7,8,9,10}"
JOBS="${3:-3}"

ROOT="${PIPELINE_ROOT:-/pipeline}"
RES="${ROOT}/${PROJECT}_results"
OUT_ROOT="${ROOT}/benchmarks/laruson_sweep/${PROJECT}/rda"

# Tuning knobs held FIXED across the ladder — only condition_pcs varies.
PREDICTORS="${PREDICTORS:-bio_1,bio_2}"
AXES="${AXES:-2}"          # rank ceiling with 2 predictors; fixing it skips the by-axis anova.cca
AXIS_ALPHA="${AXIS_ALPHA:-0.05}"
PERMS="${PERMS:-999}"
FIT_MODE="${FIT_MODE:-full}"
MAX_SNPS="${MAX_SNPS:-10000000}"
MAX_GB="${MAX_GB:-100}"
SEED="${SEED:-42}"
CPU_PER_JOB="${CPU_PER_JOB:-8}"
CAND_ADJUST="${CAND_ADJUST:-qval}"
CAND_VALUE="${CAND_VALUE:-0.1}"

die() { echo "ERROR: $*" >&2; exit 1; }

# Resolve exactly one path from a glob, or abort with the reason.
resolve1() {
    local desc="$1" pattern="$2"
    local matches=()
    # shellcheck disable=SC2206
    matches=( $(compgen -G "$pattern" || true) )
    (( ${#matches[@]} == 1 )) || die "$desc: expected exactly 1 match for '$pattern', got ${#matches[@]}${matches[*]+ (${matches[*]})}"
    [[ -s "${matches[0]}" ]] || die "$desc: resolved to an empty file: ${matches[0]}"
    echo "${matches[0]}"
}

[[ -d "$RES" ]] || die "No results directory: $RES (run mode=processing/structure/climate first)"

LFMM_FULL=$(resolve1   "lfmm_imp_full_climate"  "${RES}/_intermediate/climate_subset/lfmm_imp_full_climate.lfmm")
LFMM_PRUNED=$(resolve1 "lfmm_imp_climate"       "${RES}/_intermediate/climate_subset/lfmm_imp_climate.lfmm")
CLIMATE=$(resolve1     "climate_site_scaled"    "${RES}/climate/tables/present/climate_present_site_scaled.tsv")
VCFSNP_FULL=$(resolve1 "vcfsnp full"            "${RES}/_work/maf*/*.vcfsnp")
VCFSNP_PRUNED=$(resolve1 "vcfsnp pruned"        "${RES}/_work/maf*/ld*/*.vcfsnp")
PCA=$(resolve1         "LEA projections"        "${RES}/_work/maf*/ld*/*.pca/*.projections")
SAMPLES_ORDER=$(resolve1 "samples_order"        "${RES}/_intermediate/samples/samples_order.list")
CLIMATE_VALID=$(resolve1 "climate-valid list"   "${RES}/_intermediate/samples/coord_valid_samples.list")

# K_BEST is used by rda.R for filenames/titles only, never as an axis count.
K_BEST=$(basename "$(resolve1 "K-imputed lfmm" "${RES}/_work/maf*/*_K*imp.lfmm")" | sed -E 's/.*_K([0-9]+)imp\.lfmm/\1/')
[[ "$K_BEST" =~ ^[0-9]+$ ]] || die "Could not infer K_BEST from the imputed-matrix filename"

echo "INFO: project=$PROJECT  K_BEST=$K_BEST  predictors=$PREDICTORS  axes=$AXES  perms=$PERMS"
echo "INFO: genotypes  = $LFMM_FULL"
echo "INFO: vcfsnp     = $VCFSNP_FULL"
echo "INFO: pca        = $PCA"
echo "INFO: output     = $OUT_ROOT"

run_rung() {
    local n="$1"
    local d="${OUT_ROOT}/pc${n}"
    mkdir -p "${d}/plots"
    echo "INFO: [pc${n}] starting"
    if Rscript "${ROOT}/scripts/rda.R" \
        "$LFMM_FULL" "$LFMM_PRUNED" "$CLIMATE" \
        "$VCFSNP_FULL" "$VCFSNP_PRUNED" \
        "$PCA" "$SAMPLES_ORDER" "$CLIMATE_VALID" \
        "$PREDICTORS" "$n" "$AXES" "$AXIS_ALPHA" \
        "$PERMS" "$FIT_MODE" "$MAX_SNPS" "$MAX_GB" \
        "$SEED" "$CPU_PER_JOB" "$K_BEST" \
        "${d}/RDA_pvalues.tsv" "${d}/RDA_candidates.tsv" \
        "${d}/RDA_diagnostics.tsv" "${d}/RDA_anova.tsv" \
        "${d}/plots/" "$CAND_ADJUST" "$CAND_VALUE" > "${d}/rda.log" 2>&1
    then
        echo "INFO: [pc${n}] done — $(grep -c . "${d}/RDA_pvalues.tsv" 2>/dev/null || echo 0) lines; $(tail -1 "${d}/rda.log" | cut -c1-140)"
    else
        echo "ERROR: [pc${n}] FAILED — see ${d}/rda.log" >&2
        tail -5 "${d}/rda.log" >&2
        return 1
    fi
}

mkdir -p "$OUT_ROOT"
IFS=',' read -r -a RUNGS <<< "$COND_PCS_CSV"

# Bounded concurrency: each worker holds its own copy of the genotype matrix.
failed=0
for n in "${RUNGS[@]}"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || failed=1; done
    run_rung "$n" &
done
while (( $(jobs -rp | wc -l) > 0 )); do wait -n || failed=1; done

if (( failed )); then
    echo "WARNING: at least one rung failed — check the per-rung rda.log files" >&2
    exit 1
fi
echo "INFO: RDA ladder complete: $OUT_ROOT"
