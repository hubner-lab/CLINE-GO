#!/usr/bin/env bash
# =============================================================================
# mvp_rda_uncorrected_run.sh -- re-run the GEA RDA scan UNCORRECTED (condition_pcs = 0) on
# SS-Clines replicates, beside the corrected one-fit RDA they already carry.
#
# WHY. Journal 18: one-fit (corrected, Condition(PC1..PCk)) RDA ranks causal loci worse than the
# frozen pmax arm did (AUC-PR median 0.072 -> 0.040). pmax was max(p_partial, p_unconstrained), so
# the question is whether the uncorrected fit alone carries that ranking. Phase 3b plan:
# ~/.claude/plans/glittery-fluttering-seahorse.md. A check only -- Phase 4 panels keep corrected RDA.
#
# HOW. scripts/rda.R is run DIRECTLY, not through Snakemake: a snakemake re-run with condition_pcs 0
# would overwrite GEA/tables/methods/RDA/*_K{k}.tsv in MVP{seed}_results/ and rebuild selected_snps,
# regions and WZA from them. Instead, per seed, the exact command Snakemake ran is read back from its
# metadata (MVP{seed}_results/.snakemake/metadata/<urlsafe b64 of the output path>, key `shellcmd`),
# asserted to have the expected shape, and re-issued with exactly five substitutions:
#   arg 10  condition_pcs  k_best -> 0
#   args 20-23  pvalues / candidates / diagnostics / anova -> params_rdaunc/MVP{seed}/c1/RDA_*.tsv
#   arg 24  PLOT_DIR -> params_rdaunc/MVP{seed}/c1/plots/   (rda.R names its plots rda_*_K{k} itself)
# Every input, CPU (2), SEED (42), permutations (99), fit mode and candidate rule stay verbatim; the
# results trees are only READ. Image cline-go:latest (93bd6025) -- the image of both other arms.
#
# A seed is complete when params_rdaunc/MVP{seed}/c1/DONE exists; it is written only after rda.R
# exits 0 (a killed run can leave a partial p-table, so the table itself is not the sentinel).
# Longest seeds first; the first wave is staggered STAGGER_S apart so the ~25-29 GB rda() peaks do
# not coincide. Memory of every mvp-rdaunc-* container is sampled every 15 s into memtrace.tsv.
#
# Run under the nix docker client (the snap client cannot run in this session):
#   setsid nohup nix shell nixpkgs#docker-client -c \
#     benchmarks/mvp_rda_uncorrected_run.sh <seeds|all> > benchmarks/mvp_eval/runlogs_rdaunc/driver.log 2>&1 &
#   <seeds>  comma-separated seed ids, or `all` = the 600 of remeasure600/onefit/covariates.tsv
#   JOBS (30)  STAGGER_S (20)  IMAGE (cline-go:latest)  MEM (40g)
# =============================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EVAL="$ROOT/benchmarks/mvp_eval"
OUT="$EVAL/params_rdaunc"
LOGD="$EVAL/runlogs_rdaunc"
COV="$EVAL/remeasure600/onefit/covariates.tsv"
JOBS="${JOBS:-30}"; STAGGER_S="${STAGGER_S:-20}"; IMAGE="${IMAGE:-cline-go:latest}"; MEM="${MEM:-40g}"
SEL="${1:?seeds (comma list) or all}"
export ROOT EVAL OUT LOGD COV IMAGE MEM STAGGER_S JOBS

command -v docker > /dev/null && docker version > /dev/null 2>&1 \
    || { echo "FATAL: no working docker client -- run under 'nix shell nixpkgs#docker-client -c'"; exit 1; }
[[ -s "$COV" ]] || { echo "FATAL: missing $COV"; exit 1; }
mkdir -p "$OUT" "$LOGD"

# seed -> k_best (the K{k} of the harvested file names and rda.R's K_BEST arg)
kbest_of() { awk -F'\t' -v s="$1" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "k_best") c = i; next }
                                    $1 == s { print $c }' "$COV"; }
export -f kbest_of

# The re-issued argv, NUL-separated, or a FATAL line on stderr and exit 1.
build_argv() {
    python3 - "$1" "$2" <<'PY'
import sys, os, json, base64, shlex
seed, k = sys.argv[1], sys.argv[2]
root, out = os.environ["ROOT"], os.environ["OUT"]
pv = f"/pipeline/MVP{seed}_results/GEA/tables/methods/RDA/RDA_pvalues_K{k}.tsv"
meta = os.path.join(root, f"MVP{seed}_results/.snakemake/metadata",
                    base64.urlsafe_b64encode(pv.encode()).decode())
def fatal(msg):
    sys.stderr.write(f"FATAL seed {seed}: {msg}\n"); sys.exit(1)
if not os.path.isfile(meta): fatal(f"no Snakemake metadata for {pv}")
tok = shlex.split(json.load(open(meta))["shellcmd"])
if tok[:2] != ["Rscript", "/pipeline/scripts/rda.R"]: fatal(f"unexpected command head {tok[:2]}")
if ">" not in tok: fatal("no log redirect in shellcmd")
i = tok.index(">")
a = tok[2:i]
if len(a) != 26: fatal(f"{len(a)} positional args, expected 26")
if tok[i + 1:] != [f"/pipeline/MVP{seed}_logs/gea/assoc_rda.log", "2>&1"]: fatal(f"unexpected redirect {tok[i + 1:]}")
if a[9] != k:  fatal(f"arg 10 condition_pcs = {a[9]}, expected k_best {k}")
if a[18] != k: fatal(f"arg 19 K_BEST = {a[18]}, expected {k}")
if a[19] != pv: fatal(f"arg 20 = {a[19]}, expected {pv}")
d = f"/pipeline/benchmarks/mvp_eval/params_rdaunc/MVP{seed}/c1"
a[9] = "0"
a[19:23] = [f"{d}/RDA_pvalues.tsv", f"{d}/RDA_candidates.tsv", f"{d}/RDA_diagnostics.tsv", f"{d}/RDA_anova.tsv"]
a[23] = f"{d}/plots/"
sys.stdout.write("\0".join(["Rscript", "/pipeline/scripts/rda.R"] + a))
PY
}
export -f build_argv

run_one() {   # $1 = queue index, $2 = seed
    local idx="$1" seed="$2" k d log t0 rc
    d="$OUT/MVP$seed/c1"; log="$LOGD/MVP$seed.log"
    [[ -f "$d/DONE" ]] && { echo "[$seed] SKIP (DONE)"; return 0; }
    k="$(kbest_of "$seed")"
    [[ -n "$k" ]] || { echo "[$seed] FATAL: no k_best in $COV"; return 1; }
    local -a argv
    readarray -d '' argv < <(build_argv "$seed" "$k") || true
    (( ${#argv[@]} == 28 )) || { echo "[$seed] FATAL: argv not built (see stderr)"; return 1; }
    (( idx < JOBS )) && sleep $(( idx * STAGGER_S ))
    mkdir -p "$d/plots"
    t0=$(date +%s); echo "[$seed] START $(date -Is) k_best=$k"
    docker run --rm --name "mvp-rdaunc-$seed" --user "$(id -u):$(id -g)" -e USER=adaptogene \
        -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory="$MEM" -v "$ROOT":/pipeline "$IMAGE" \
        "${argv[@]}" > "$log" 2>&1
    rc=$?
    if (( rc == 0 )); then
        printf "seed\t%s\nk_best\t%s\ncondition_pcs\t0\nseconds\t%s\nfinished\t%s\n" \
            "$seed" "$k" "$(( $(date +%s) - t0 ))" "$(date -Is)" > "$d/DONE"
        echo "[$seed] DONE $(date -Is) $(( $(date +%s) - t0 )) s"
    else
        echo "[$seed] FAILED exit $rc after $(( $(date +%s) - t0 )) s (log $log)"
    fi
    return $rc
}
export -f run_one

memtrace() {
    local f="$LOGD/memtrace.tsv"
    [[ -s "$f" ]] || printf "ts\tcontainer\tmem_gib\n" > "$f"
    while :; do
        docker stats --no-stream --format '{{.Name}}\t{{.MemUsage}}' 2>/dev/null | grep '^mvp-rdaunc-' |
            awk -F'\t' -v t="$(date -Is)" '
                function gib(x,  v, u) { v = x + 0; u = x; gsub(/[0-9.]/, "", u)
                    if (u ~ /^Ki/) return v / 1048576; if (u ~ /^Mi/) return v / 1024
                    if (u ~ /^Gi/) return v; if (u ~ /^Ti/) return v * 1024; return v / 1073741824 }
                { split($2, m, " / "); printf "%s\t%s\t%.3f\n", t, $1, gib(m[1]) }' >> "$f"
        sleep 15
    done
}

# ---- seed list, longest first (p-table rows of the corrected arm)
if [[ "$SEL" == "all" ]]; then
    mapfile -t SEEDS < <(awk -F'\t' 'NR > 1 { print $1 }' "$COV")
else
    IFS=',' read -ra SEEDS <<< "$SEL"
fi
QUEUE="$LOGD/queue.tsv"
for s in "${SEEDS[@]}"; do
    f="$EVAL/params_onefit/MVP$s/c1/RDA_pvalues.tsv"
    [[ -s "$f" ]] || { echo "FATAL: no corrected RDA table for $s ($f)"; exit 1; }
    printf "%s\t%s\n" "$s" "$(( $(wc -l < "$f") - 1 ))"
done | sort -t$'\t' -k2,2nr > "$QUEUE"
echo "rdaunc driver: ${#SEEDS[@]} seed(s), JOBS=$JOBS STAGGER_S=$STAGGER_S IMAGE=$IMAGE MEM=$MEM start $(date -Is)"

memtrace & MT=$!
trap 'kill $MT 2>/dev/null' EXIT
awk -F'\t' '{ print NR - 1, $1 }' "$QUEUE" | xargs -P "$JOBS" -n 2 bash -c 'run_one "$0" "$1"'
rc=$?
n_done=$(for s in "${SEEDS[@]}"; do [[ -f "$OUT/MVP$s/c1/DONE" ]] && echo x; done | wc -l)
echo "rdaunc driver: $n_done of ${#SEEDS[@]} DONE, xargs exit $rc, end $(date -Is)"
(( n_done == ${#SEEDS[@]} ))
