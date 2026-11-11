#!/usr/bin/env bash
# frdk.sh <paper|latest> <job> [docker-run extra flags ...] -- <command ...>
# Final redraw (lab, 2026-11-11): one compliant container. A copy of benchmarks/mvp_land_dk.sh (paper) and of
# work/closeout_lab/results/scripts/tools/dkc.sh (latest), with the TEMP root of this job:
#   paper  = image cline-go:paper-d150 (id 93bd6025 asserted), pipeline code overlaid read-only from the c80710e export,
#            as the Data S1 jobs B and F ran (benchmarks/mvp_land_*);
#   latest = image cline-go:latest (id 005aa13b asserted), as the closeout copies ran (build_paper_pack_d150.R,
#            plot_random_climate_ms_d150.R, build_fig3.R only=as_pred14_v2) and the audit 3b s_qc run.
# Named fr_<job>, --rm, -u, USER=pipeline, /tmp = /mnt/data/eugene/TEMP/final_redraw/<job>_<timestamp> (new per run; it is
# NOT removed here: the permission mode refuses rm -rf, the user removes the TEMP root at the end), TMPDIR=HOME=/tmp,
# OMP/OPENBLAS threads 1, nix client --offline, working directory $FR_WORKDIR (default /pipeline).
set -uo pipefail
VAR="${1:?paper|latest}"; shift
JOB="${1:?job}"; shift
EXTRA=()
while [[ $# -gt 0 && "$1" != "--" ]]; do EXTRA+=("$1"); shift; done
[[ "${1:-}" == "--" ]] || { echo "frdk.sh: missing -- before the command" >&2; exit 2; }
shift
[[ $# -gt 0 ]] || { echo "frdk.sh: empty command" >&2; exit 2; }
REPO=/mnt/data/eugene/ADAPTOGENE
DOCKER=(nix shell --offline nixpkgs#docker-client -c docker)
case $VAR in
  paper)  IMG=cline-go:paper-d150; WANT=sha256:93bd6025e77b18e9b97ae1bea4a37e9c43a5b64c8acc067819632e09d550546f
          CODE=/mnt/data/eugene/TEMP/mvp_land/code_c80710e
          [[ "$(cat "$CODE/COMMIT" 2>/dev/null)" == c80710e ]] || { echo "frdk.sh: code export missing at $CODE" >&2; exit 2; }
          OVERLAY=(-v "$CODE/Snakefile":/pipeline/Snakefile:ro -v "$CODE/workflow":/pipeline/workflow:ro
                   -v "$CODE/scripts":/pipeline/scripts:ro) ;;
  latest) IMG=cline-go:latest; WANT=sha256:005aa13b65f0877e310509c79b17bb52565032b21d19168e81bcfc81c21baa58
          OVERLAY=() ;;
  *) echo "frdk.sh: variant must be paper or latest" >&2; exit 2 ;;
esac
id=$("${DOCKER[@]}" images --no-trunc --format '{{.ID}}' "$IMG" 2>/dev/null | head -1)
[[ "$id" == "$WANT" ]] || { echo "frdk.sh: $IMG is $id, expected $WANT" >&2; exit 2; }
T="/mnt/data/eugene/TEMP/final_redraw/${JOB}_$(date +%Y%m%dT%H%M%S)"   # one per run: reruns of a step get their own
[[ ! -e "$T" ]] || { echo "frdk.sh: $T exists" >&2; exit 2; }
mkdir -p "$T"
"${DOCKER[@]}" run --rm --name "fr_$JOB" -u "$(id -u):$(id -g)" -e USER=pipeline \
  -v "$T":/tmp -e TMPDIR=/tmp -e HOME=/tmp \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e R_DATATABLE_NUM_THREADS=2 \
  "${EXTRA[@]}" \
  -v "$REPO":/pipeline "${OVERLAY[@]}" \
  -w "${FR_WORKDIR:-/pipeline}" "$IMG" "$@"
