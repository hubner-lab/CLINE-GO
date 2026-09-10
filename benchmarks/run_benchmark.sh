#!/usr/bin/env bash
# run_benchmark.sh — Headless end-to-end CLINE-GO pipeline run for benchmarking.
#
# Chains Snakemake modes sequentially (never parallel — all modes share one results dir/lock).
# Runs: processing → prestructure → structure → gea → [promote_snp_set] → maladaptation
#
# Usage:
#   bash benchmarks/run_benchmark.sh <configfile> [<set_name>] [<docker_image>] [<cpu>]
#
#   configfile    Path to config YAML (e.g. config_SIMDATA_custom.yaml)
#   set_name      SNP set name for maladaptation (default: "gea_sig")
#                 Must match Maladaptation.snp_sets in the configfile.
#   docker_image  Docker image tag (default: cline-go:latest)
#   cpu           Number of cores (default: 4)
#
# DEFERRED (Phase 4/5 — requires Gate G1 dataset deep-research):
#   - Reciprocal-transplant garden loop (multiple future scenarios per garden)
#   - Two-axis evaluation: eval_offset.R (offset vs simulated fitness),
#     eval_detection.R (precision/recall vs known causal loci)
#   These will be added once the Láruson archive format is confirmed.

set -euo pipefail

CONFIGFILE="${1:?Usage: $0 <configfile> [set_name] [docker_image] [cpu]}"
SET_NAME="${2:-gea_sig}"
DOCKER_IMAGE="${3:-cline-go:latest}"
CPU="${4:-4}"

# ── Derive project name from config ──────────────────────────────────────────
PROJECT=$(python3 -c "
import yaml, sys
with open('$CONFIGFILE') as f:
    c = yaml.safe_load(f)
print(c.get('project_name', c.get('Input', {}).get('project_name', 'BENCHMARK')))
" 2>/dev/null || grep -E '^project_name:' "$CONFIGFILE" | head -1 | awk '{print $2}' | tr -d '"'"'')

if [[ -z "$PROJECT" ]]; then
    echo "ERROR: Could not extract project_name from $CONFIGFILE" >&2
    exit 1
fi

RESULTS_DIR="${PROJECT}_results"
INTER_DIR="${RESULTS_DIR}/_intermediate"
SNP_SETS_DIR="${INTER_DIR}/snp_sets"

echo "========================================================"
echo "  CLINE-GO headless benchmark"
echo "  Config:   $CONFIGFILE"
echo "  Project:  $PROJECT"
echo "  SNP set:  $SET_NAME"
echo "  Image:    $DOCKER_IMAGE"
echo "  CPUs:     $CPU"
echo "  Started:  $(date)"
echo "========================================================"

# ── Docker run helper ─────────────────────────────────────────────────────────
run_mode() {
    local mode="$1"
    echo ""
    echo "── mode=$mode ─────────────────────────────────────────"
    docker run \
        --user "$(id -u):$(id -g)" \
        --rm \
        --memory=20g \
        -v "$PWD:/pipeline" \
        "$DOCKER_IMAGE" \
        snakemake -c"$CPU" -s Snakefile \
            --config mode="$mode" \
            --configfile "$CONFIGFILE" \
            --scheduler greedy
    echo "── mode=$mode complete ✓"
}

# ── Run all modes ─────────────────────────────────────────────────────────────
run_mode processing
run_mode prestructure
run_mode structure
run_mode gea

# ── Promote GEA SNPs → maladaptation SNP set ─────────────────────────────────
GEA_SNPS="${RESULTS_DIR}/GEA/tables/selected_snps.tsv"
if [[ ! -f "$GEA_SNPS" ]]; then
    echo "WARNING: GEA selected_snps.tsv not found at $GEA_SNPS" >&2
    echo "         Check that mode=gea produced significant SNPs." >&2
    echo "         Maladaptation mode will only run if snp_sets exist." >&2
else
    echo ""
    echo "── promote_snp_set: $GEA_SNPS → $SET_NAME ────────────"
    docker run \
        --user "$(id -u):$(id -g)" \
        --rm \
        --memory=4g \
        -v "$PWD:/pipeline" \
        "$DOCKER_IMAGE" \
        Rscript /pipeline/scripts/promote_snp_set.R \
            "/pipeline/${GEA_SNPS}" \
            "$SET_NAME" \
            "/pipeline/${SNP_SETS_DIR}" \
            "GEA"
    echo "── promote_snp_set complete ✓"
fi

run_mode maladaptation

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo "  Benchmark complete: $(date)"
echo "  Results in: $RESULTS_DIR/"
echo ""
echo "  Key outputs:"
echo "    GEA sig SNPs:     $GEA_SNPS"
echo "    SNP set:          $SNP_SETS_DIR/$SET_NAME/selected_snps.tsv"
for method in gradient_forest geometric_offset; do
    site_tsv=$(find "${RESULTS_DIR}/Maladaptation/tables/${method}" \
                    -name "genetic_offset_site.tsv" 2>/dev/null | head -1 || true)
    [[ -n "$site_tsv" ]] && echo "    ${method} offsets: $site_tsv"
done
echo ""
echo "  DEFERRED (Phase 4/5 — needs Gate G1 dataset deep-research):"
echo "    - Reciprocal-transplant garden loop"
echo "    - eval_offset.R  (offset vs simulated fitness)"
echo "    - eval_detection.R (precision/recall vs causal loci)"
echo "========================================================"
