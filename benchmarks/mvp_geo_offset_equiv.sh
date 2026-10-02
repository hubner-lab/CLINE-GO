#!/usr/bin/env bash
# mvp_geo_offset_equiv.sh SEED PANEL RUN SCRIPT -- replay one seed's geometric_offset job with SCRIPT
#
# Why: proving that the single-fit scripts/geometric_offset.R (2026-10-01) writes the same
# offsets as the per-scenario LEA::genetic.gap() version needs both scripts run on the SAME
# argv. The argv is read back from Snakemake's own record of the job
# (MVP{seed}_results/.snakemake/metadata/<urlsafe b64 of an output path>, key `shellcmd`),
# so nothing is reconstructed by hand. Every OUTPUT argument (by position: 14-17, 19) is
# redirected under $EQ_DIR/RUN/ (same relative path), the script path is replaced by SCRIPT,
# and the results tree is only read.
#
#   SEED    replicate id            PANEL  snp set (run_label), e.g. truth
#   RUN     label of this replay    SCRIPT path inside the container (/pipeline/...)
# Env: PIPELINE_ROOT (host repo), EQ_DIR (host dir under the repo), BLAS_THREADS (default 2)
set -euo pipefail
SEED="$1"; PANEL="$2"; RUN="$3"; SCRIPT="$4"
ROOT="${PIPELINE_ROOT:?export PIPELINE_ROOT}"
EQ_DIR="${EQ_DIR:?export EQ_DIR}"
BLAS_THREADS="${BLAS_THREADS:-2}"
case "$EQ_DIR" in "$ROOT"/*) ;; *) echo "FATAL: EQ_DIR must be under $ROOT" >&2; exit 1;; esac
EQ_IN="/pipeline/${EQ_DIR#"$ROOT"/}/$RUN"
mkdir -p "$EQ_DIR/$RUN"

ARGV_FILE="$EQ_DIR/$RUN/argv.txt"
ROOT="$ROOT" SEED="$SEED" PANEL="$PANEL" EQ_IN="$EQ_IN" SCRIPT="$SCRIPT" EQ_HOST="$EQ_DIR/$RUN" \
python3 - > "$ARGV_FILE" <<'PY'
import os, sys, json, base64, shlex
root, seed, panel = os.environ["ROOT"], os.environ["SEED"], os.environ["PANEL"]
eq_in, script, eq_host = os.environ["EQ_IN"], os.environ["SCRIPT"], os.environ["EQ_HOST"]
res = f"/pipeline/MVP{seed}_results/"
diag = f"{res}Maladaptation/tables/geometric_offset/{panel}_nospatial/geometric_offset_diagnostics.tsv"
meta = os.path.join(root, f"MVP{seed}_results/.snakemake/metadata",
                    base64.urlsafe_b64encode(diag.encode()).decode())
if not os.path.isfile(meta): sys.exit(f"FATAL: no Snakemake metadata for {diag}")
tok = shlex.split(json.load(open(meta))["shellcmd"])
if tok[:2] != ["Rscript", "/pipeline/scripts/geometric_offset.R"]: sys.exit(f"FATAL: command head {tok[:2]}")
i = tok.index(">")
a = tok[2:i]
if len(a) != 19: sys.exit(f"FATAL: {len(a)} positional args, expected 19")
# OUTPUT arguments by position (scripts/geometric_offset.R header): 14 site TSVs, 15 map TSVs,
# 16 offset rasters (these live under _intermediate/, NOT Maladaptation/ -- a prefix test on
# Maladaptation/ misses them and the replay overwrites the real tree's rasters), 17 importance
# PNG, 19 diagnostics. Every one of them moves under the replay dir; inputs stay untouched.
OUT_POS = {13, 14, 15, 16, 18}          # 0-based into a
out = []
for j, x in enumerate(a):
    parts = x.split(",")
    if j in OUT_POS:
        if not all(p.startswith(res) for p in parts): sys.exit(f"FATAL: output arg {j+1} outside {res}: {x[:80]}")
        parts = [eq_in + "/" + p[len(res):] for p in parts]
        for p in parts:   # create the host-side directory of every redirected output
            os.makedirs(os.path.dirname(eq_host + "/" + p[len(eq_in) + 1:]), exist_ok=True)
    out.append(",".join(parts))
for j in OUT_POS:
    if any(p.startswith(res) for p in out[j].split(",")): sys.exit(f"FATAL: output arg {j+1} still in the results tree")
print("\0".join([script] + out), end="")
PY

mapfile -d '' ARGS < "$ARGV_FILE"
echo "[$RUN] seed $SEED panel $PANEL script ${ARGS[0]}  start $(date -Is)"
t0=$(date +%s)
nix shell nixpkgs#docker-client -c docker run --rm --name "geo-equiv-${SEED}-${RUN}" \
  --user "$(id -u):$(id -g)" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" -e OMP_NUM_THREADS="$BLAS_THREADS" \
  --cpus=4 --memory=32g -v "$ROOT":/pipeline cline-go:latest \
  Rscript "${ARGS[@]}" > "$EQ_DIR/$RUN/run.log" 2>&1
echo "[$RUN] done rc=$? in $(( $(date +%s) - t0 )) s"
