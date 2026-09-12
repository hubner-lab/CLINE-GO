#!/usr/bin/env bash
# One command for the whole test suite.
#
#   tests/run_all.sh                          # the merge gate: all unit suites
#   tests/run_all.sh --heavy                  # + the tool-dependent wrapper tests
#   tests/run_all.sh --invariants SIMDATA_results     # path relative to the repo root
#   tests/run_all.sh --image cline-go:latest
#
# Runs every suite to completion and reports a one-line summary per suite, then
# exits non-zero if any suite failed. It deliberately does NOT stop at the first
# failure: the point is to learn everything that broke in one pass.
#
# THREE GATES, THREE DIFFERENT CONTRACTS. Do not confuse them:
#
#   default    MUST be green. Hermetic, no genomics tool, seconds to minutes.
#   --heavy    EXPECTED green. Needs plink / LEA / EMMAX / vcftools; separated
#              from the default gate by RUNTIME and tool dependency, NOT by
#              correctness. Run it before publishing a pipeline version. A missing
#              tool FAILS rather than skips — inside the image they are all
#              present, so absence means the image is wrong.
#   --invariants  EXPECTED RED. See the note below.
#
# ON --invariants: the validator is EXPECTED to be red on the current
# SIMDATA_results tree. Every violation it reports there corresponds to a defect
# already filed in docs/pipeline_improvement_requests.md (stale GEAxGWAS
# pairwise table, overlap_traits/overlap_snps, duplicated collapsed-gene
# columns, n_climate_variables). That is why it is opt-in rather than part of
# the default gate — the default must stay meaningfully green.

set -uo pipefail

IMAGE="cline-go:latest"
RUN_HEAVY=0
INVARIANTS_DIR=""
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --heavy)      RUN_HEAVY=1;             shift 1 ;;
        --invariants) INVARIANTS_DIR="${2:-}"; shift 2 ;;   # path RELATIVE to the repo root
        --image)      IMAGE="${2:-}";          shift 2 ;;
        -h|--help)    sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# -v is required for every suite: the image ships no pipeline code.
# --name so a killed docker CLIENT leaves a container that can be stopped by
# name rather than hunted for in docker ps (see CLAUDE.md).
docker_run() {
    local name="$1"; shift
    docker run --rm --name "$name" \
        --user "$(id -u):$(id -g)" -e USER=pipeline \
        -v "${REPO_ROOT}:/pipeline" "$IMAGE" "$@"
}

FAILED=()
SUMMARY=()

run_suite() {
    local label="$1" name="$2"; shift 2
    echo
    echo "══ ${label} ═════════════════════════════════════════════"
    if docker_run "$name" "$@"; then
        SUMMARY+=("PASS  ${label}")
    else
        SUMMARY+=("FAIL  ${label}")
        FAILED+=("$label")
    fi
}

# 1. Pipeline libraries: scripts/R/lib + scripts/R/utils.
#    NOT via scripts/clinego.app/dev.R and NOT with .R_libs_dev on .libPaths() —
#    that carries testthat 3.3.2 while the image pins 3.2.3, and the two
#    disagree about test_dir()'s return shape (tests/run_tests.R:15-17).
run_suite "pipeline libs (tests/)" clinego_tests_libs \
    Rscript /pipeline/tests/run_tests.R

# 2. The Shiny app package.
run_suite "shiny app (clinego.app)" clinego_tests_app \
    Rscript -e 'setwd("/pipeline/scripts/clinego.app/tests"); source("testthat.R")'

# 3. Python: scripts/*.py + workflow/methods/.
#    stdlib unittest, deliberately not pytest: the image has python3.12 + numpy
#    + pandas + scipy already, while pytest is absent and PEP 668
#    (/usr/lib/python3.12/EXTERNALLY-MANAGED) plus cleaned apt lists make adding
#    it a Dockerfile layer for no gain — nothing under test needs fixtures or
#    parametrize. -t puts tests/python on sys.path so `import _support` resolves
#    without an __init__.py; -B keeps __pycache__ out of the mounted repo.
run_suite "python (tests/python/)" clinego_tests_python \
    python3 -B -m unittest discover -s /pipeline/tests/python -t /pipeline/tests/python

# 4. Heavy tier: wrapper tests needing plink / LEA / EMMAX / vcftools. Opt-in for
#    runtime, not for correctness — expected green. Its own testthat root, so the
#    default gate neither sources nor skips it (tests/run_heavy.R explains why).
if [[ "$RUN_HEAVY" -eq 1 ]]; then
    run_suite "heavy wrappers (tests/heavy/)" clinego_tests_heavy \
        Rscript /pipeline/tests/run_heavy.R
fi

# 5. Invariants over a real results tree — opt-in, see the header note.
if [[ -n "$INVARIANTS_DIR" ]]; then
    # The container sees the repo at /pipeline, so the argument must be a path
    # relative to the repo root. Accept an absolute host path too by stripping
    # the repo prefix — "/pipeline/${abs}" would otherwise be nonsense.
    rel="${INVARIANTS_DIR#./}"
    rel="${rel#"${REPO_ROOT}/"}"
    if [[ "$rel" = /* ]]; then
        echo "--invariants must be a path inside the repo (got: ${INVARIANTS_DIR})" >&2
        exit 2
    fi
    if [[ ! -d "${REPO_ROOT}/${rel}" ]]; then
        echo "--invariants: no such directory: ${REPO_ROOT}/${rel}" >&2
        exit 2
    fi
    run_suite "invariants (${rel})" clinego_tests_invariants \
        Rscript /pipeline/scripts/check_invariants.R "/pipeline/${rel}"
fi

echo
echo "══ summary ══════════════════════════════════════════════"
printf '%s\n' "${SUMMARY[@]}"

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo
    echo "FAILED: ${FAILED[*]}"
    exit 1
fi
echo
echo "All suites passed."
