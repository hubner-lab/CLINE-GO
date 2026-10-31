#!/usr/bin/env bash
# Build the documentation site locally, exactly as Read the Docs does (warnings are errors).
#
#   nix shell nixpkgs#python312 -c documentation/build_docs.sh <outdir> <venv dir>
#
# <venv dir> is created with the python3 on PATH when it does not exist yet, then reused;
# the pinned requirements.txt is installed into it on every run (a no-op when satisfied).
# Nothing is written inside the repository: both paths must lie outside it, and
# PYTHONDONTWRITEBYTECODE keeps conf.py's imports from leaving __pycache__/ behind.
# -E rebuilds from scratch every time, as Read the Docs does: with the doctrees reused, a warning
# raised while reading an unchanged page is not raised again, and a red build reruns green.
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "usage: $0 <outdir> <venv dir>" >&2
  exit 2
fi

DOCS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
REPO_DIR=$(cd "$DOCS_DIR/.." && pwd -P)
OUT_DIR=$(realpath -m "$1")
VENV_DIR=$(realpath -m "$2")

for p in "$OUT_DIR" "$VENV_DIR"; do
  case "$p" in
    "$REPO_DIR" | "$REPO_DIR"/*)
      echo "build_docs.sh: $p is inside the repository ($REPO_DIR); choose a path outside it" >&2
      exit 2
      ;;
  esac
done

command -v python3 >/dev/null || { echo "build_docs.sh: no python3 on PATH" >&2; exit 1; }

export PYTHONDONTWRITEBYTECODE=1

if [ ! -x "$VENV_DIR/bin/python" ]; then
  echo "build_docs.sh: creating venv $VENV_DIR with $(python3 --version)"
  python3 -m venv "$VENV_DIR"
fi
"$VENV_DIR/bin/python" -m pip install --quiet --disable-pip-version-check \
  -r "$DOCS_DIR/requirements.txt"

mkdir -p "$OUT_DIR"
"$VENV_DIR/bin/python" -m sphinx -E -W --keep-going -b html "$DOCS_DIR" "$OUT_DIR"
echo "build_docs.sh: site written to $OUT_DIR/index.html"
