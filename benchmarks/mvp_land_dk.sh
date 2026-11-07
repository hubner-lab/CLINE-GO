#!/usr/bin/env bash
# mvp_land_dk.sh <job> [docker-run extra flags ...] -- <command ...>
# One compliant container for the landscapes pack: a copy of work/mvp_land/dk.sh with the TEMP root
# of the pack.
#   - image cline-go:paper-d150, its id asserted (93bd6025, the image of offset14 and offset15);
#   - the repo at /pipeline, the pipeline code overlaid read-only from the c80710e export (Snakefile,
#     workflow/, scripts/), so scripts that source /pipeline/scripts/... read the code of the runs;
#   - named container lpack_<job>, --rm, TEMP /tmp under /mnt/data/eugene/TEMP/landscapes/<job>
#     (removed after exit), thread caps; pass --cpus/--memory as extra flags;
#   - DK_WORKDIR overrides the working directory (default /pipeline).
# Exit status = the container's. A later -e for the same variable overrides the default cap.
set -uo pipefail
JOB="${1:?job}"; shift
[[ "$JOB" != reference ]] || { echo "mvp_land_dk.sh: 'reference' holds the published tables" >&2; exit 2; }
EXTRA=()
while [[ $# -gt 0 && "$1" != "--" ]]; do EXTRA+=("$1"); shift; done
[[ "${1:-}" == "--" ]] || { echo "mvp_land_dk.sh: missing -- before the command" >&2; exit 2; }
shift
[[ $# -gt 0 ]] || { echo "mvp_land_dk.sh: empty command" >&2; exit 2; }
REPO=/mnt/data/eugene/ADAPTOGENE
CODE=/mnt/data/eugene/TEMP/mvp_land/code_c80710e
IMG=cline-go:paper-d150
WANT=sha256:93bd6025e77b18e9b97ae1bea4a37e9c43a5b64c8acc067819632e09d550546f
DOCKER=(nix shell nixpkgs#docker-client -c docker)
[[ "$(cat "$CODE/COMMIT" 2>/dev/null)" == c80710e ]] || { echo "mvp_land_dk.sh: pinned code export missing at $CODE" >&2; exit 2; }
id=$("${DOCKER[@]}" images --no-trunc --format '{{.ID}}' "$IMG" 2>/dev/null | head -1)
[[ "$id" == "$WANT" ]] || { echo "mvp_land_dk.sh: $IMG is $id, expected $WANT" >&2; exit 2; }
T="/mnt/data/eugene/TEMP/landscapes/$JOB"
mkdir -p "$T"
"${DOCKER[@]}" run --rm --name "lpack_$JOB" -u "$(id -u):$(id -g)" -e USER=pipeline \
  -v "$T":/tmp -e TMPDIR=/tmp -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e R_DATATABLE_NUM_THREADS=2 \
  "${EXTRA[@]}" \
  -v "$REPO":/pipeline \
  -v "$CODE/Snakefile":/pipeline/Snakefile:ro \
  -v "$CODE/workflow":/pipeline/workflow:ro \
  -v "$CODE/scripts":/pipeline/scripts:ro \
  -w "${DK_WORKDIR:-/pipeline}" "$IMG" "$@"
rc=$?
rm -rf "$T"
exit $rc
