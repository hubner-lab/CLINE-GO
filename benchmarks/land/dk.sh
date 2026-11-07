#!/usr/bin/env bash
# dk.sh <run> <job> [docker-run extra flags ...] -- <command ...>
# work/ld_window_lab/dk.sh pinned to the SS-Clines LD-window code (plan nifty-prancing-grove):
#   - image cline-go:paper-d150, its id asserted (93bd6025, the image of offset14);
#   - the repo at /pipeline (results, configs, data, benchmarks at main) with the pipeline code
#     overlaid read-only from the c80710e export: Snakefile, workflow/, scripts/;
#   - named container <job>, --rm, TEMP /tmp under /mnt/data/eugene/TEMP/mvp_land/<run>/<job>
#     (removed after exit), thread caps; pass --cpus/--memory as extra flags.
# Exit status = the container's.
set -uo pipefail
RUN="${1:?run}"; JOB="${2:?job}"; shift 2
EXTRA=()
while [[ $# -gt 0 && "$1" != "--" ]]; do EXTRA+=("$1"); shift; done
[[ "${1:-}" == "--" ]] || { echo "dk.sh: missing -- before the command" >&2; exit 2; }
shift
[[ $# -gt 0 ]] || { echo "dk.sh: empty command" >&2; exit 2; }
REPO=/mnt/data/eugene/ADAPTOGENE
CODE=/mnt/data/eugene/TEMP/mvp_land/code_c80710e
IMG=cline-go:paper-d150
WANT=sha256:93bd6025e77b18e9b97ae1bea4a37e9c43a5b64c8acc067819632e09d550546f
DOCKER=(nix shell nixpkgs#docker-client -c docker)
[[ "$(cat "$CODE/COMMIT" 2>/dev/null)" == c80710e ]] || { echo "dk.sh: pinned code export missing at $CODE" >&2; exit 2; }
id=$("${DOCKER[@]}" image inspect "$IMG" --format '{{.Id}}' 2>/dev/null)
[[ "$id" == "$WANT" ]] || { echo "dk.sh: $IMG is $id, expected $WANT" >&2; exit 2; }
T="/mnt/data/eugene/TEMP/mvp_land/$RUN/$JOB"
mkdir -p "$T"
"${DOCKER[@]}" run --rm --name "$JOB" -u "$(id -u):$(id -g)" -e USER=pipeline \
  -v "$T":/tmp -e TMPDIR=/tmp -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e R_DATATABLE_NUM_THREADS=2 \
  "${EXTRA[@]}" \
  -v "$REPO":/pipeline \
  -v "$CODE/Snakefile":/pipeline/Snakefile:ro \
  -v "$CODE/workflow":/pipeline/workflow:ro \
  -v "$CODE/scripts":/pipeline/scripts:ro \
  -w /pipeline "$IMG" "$@"
rc=$?
rm -rf "$T"
exit $rc
