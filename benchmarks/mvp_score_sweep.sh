#!/usr/bin/env bash
# =============================================================================
# mvp_score_sweep.sh -- score the harvested parameter cells against the truth tables.
#
#   benchmarks/mvp_score_sweep.sh [SEEDS_CSV|all] [SEED_JOBS] [CELLS]
#
# Consumes benchmarks/mvp_eval/params/MVP{seed}/c{i}/{METHOD}_pvalues.tsv (written by
# mvp_run_sweep.sh) and writes into benchmarks/mvp_eval/sweep/MVP{seed}/.
#
# THE SCORING CONTRACT (benchmarks/MVP_README.md:41-62) -- getting this wrong silently
# converts real detections into false positives, so it is encoded here once:
#
#   EMMAX / LFMM / BLINK   truth_temp + bio_1     the partially confounded axis
#                          truth_sal  + bio_2     the orthogonal control axis
#                          truth_any  + all       union
#   RDA                    truth_any  + all       ONLY -- rda.R writes a single
#                                                 climate_multivariate column, there is no
#                                                 per-axis p-value to score
#
# Three bundles per cell:
#   any   4 methods on the full tables, truth_any. The only bundle RDA can join, and the
#         headline multi-method comparison.
#   temp  3 univariate methods on bio_1-only tables, truth_temp.
#   sal   3 univariate methods on bio_2-only tables, truth_sal.
# The per-axis bundles exist because sweep_thresholds.R hardcodes load_pvalues(..., "all")
# and has no --traits switch: restricting the axis means restricting the table.
#
# The four 1-trait seeds (1231858, 1232908, 1233133, 1233208) have ZERO causal loci on the
# salinity axis. That is a feature -- bio_2 is a pure false-positive control there. Rank mode
# is skipped for them (eval_detection.R:74 hard-errors with no causal loci in the tested set);
# sweep_thresholds.R has no such guard and is run anyway, where n_called / fp_background are
# the meaningful columns and precision / recall are correctly NaN.
# =============================================================================
set -uo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
SEEDS_ARG="${1:-all}"
SEED_JOBS="${2:-7}"
CELLS="${3:-5}"

MANIFEST="$PIPELINE_ROOT/benchmarks/mvp_seeds.tsv"
PARAMS_DIR="${PARAMS_DIR:-$PIPELINE_ROOT/benchmarks/mvp_eval/params}"
SWEEP_DIR="${SWEEP_DIR:-$PIPELINE_ROOT/benchmarks/mvp_eval/sweep}"

# Method lists are env-overridable so journal 06's 11-method run reuses this driver rather
# than forking it -- the scoring contract encoded in the header is the thing that must not
# be duplicated. Defaults reproduce journal 05 exactly.
#   SWEEP_UNIVARIATE  methods with one p-value column PER PREDICTOR: EMMAX, LFMM and every
#                     GAPIT model. These get the per-axis bundles (bio_1/temp, bio_2/sal).
#   SWEEP_ALL_METHODS the full set including RDA, which is any-axis only.
# COMBINE_WINDOW_KB / TOP_MAX are passed straight through to sweep_thresholds.R.
read -r -a UNIVARIATE  <<< "${SWEEP_UNIVARIATE:-EMMAX LFMM BLINK}"
read -r -a ALL_METHODS <<< "${SWEEP_ALL_METHODS:-EMMAX LFMM RDA BLINK}"
COMBINE_WINDOW_KB="${COMBINE_WINDOW_KB:-0}"
TOP_MAX="${TOP_MAX:-5000}"

# Container-side twins of PARAMS_DIR / SWEEP_DIR. The host paths above are what this script
# stats and writes through the bind mount; these are what the R scripts INSIDE the container
# see. They must be set together -- a host override with a stale container path would read
# journal 05's tables while writing journal 06's, which no error would catch.
CONT_PARAMS="${CONT_PARAMS:-/pipeline/benchmarks/mvp_eval/params}"
CONT_SWEEP="${CONT_SWEEP:-/pipeline/benchmarks/mvp_eval/sweep}"
SWEEP_CELLS_TSV="${SWEEP_CELLS_TSV:-benchmarks/mvp_sweep_cells.tsv}"

# PROJ_SUFFIX selects a parallel arm (its own project_name / _results / params tree).
# TRUTH_PROJ names the directory holding that arm's truth tables. They are separate knobs on
# purpose: a MAF arm needs MAF-MATCHED truth (benchmarks/mvp_filter_truth_maf.R), because
# scoring a MAF-filtered VCF against the unfiltered truth divides recall by loci the VCF no
# longer contains -- the defect that disqualified the Laruson benchmark. Defaults keep the
# base arm scoring against its own tables exactly as before.
PROJ_SUFFIX="${PROJ_SUFFIX:-}"

mkdir -p "$SWEEP_DIR"
DOCKER=(nix shell nixpkgs#docker-client -c docker)

if [[ "$SEEDS_ARG" == "all" ]]; then
    mapfile -t SEEDS < <(awk 'NR>1 && NF {print $1}' "$MANIFEST")
else
    IFS=',' read -r -a SEEDS <<< "$SEEDS_ARG"
fi

R() {   # R <script> <args...>   -- run an R script inside the container
    local script="$1"; shift
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
        -e SWEEP_UNIVARIATE="${UNIVARIATE[*]}" -e SWEEP_ALL_METHODS="${ALL_METHODS[*]}" \
        -e SWEEP_DIR="$CONT_SWEEP" -e PARAMS_DIR_CONTAINER="$CONT_PARAMS" \
        -e SWEEP_CELLS_TSV="$SWEEP_CELLS_TSV" \
        --cpus=4 --memory=32g -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
        Rscript "/pipeline/benchmarks/$script" "$@"
}

n_causal() { awk 'NR>1 && $3=="causal"' "$1" | wc -l; }

# ---- TRUTH-JOIN GATE ----------------------------------------------------------------
# The pilot's own gate (MVP_README.md:146): every causal locus in the truth table must be
# present in the scored p-value table. If chromosome naming, coordinate derivation or the
# MAF filter ever drift, causal loci silently vanish from the denominator and every recall
# number becomes optimistic in a way nothing else would reveal. n_causal_testable is emitted
# by eval_detection.R for exactly this purpose.
truth_join_gate() {   # truth_join_gate <rank.tsv> <n_temp> <n_sal> <n_any>
    awk -F'\t' -v nt="$2" -v ns="$3" -v na="$4" '
        $4 == "n_causal_testable" {
            n = split($2, p, "|"); ax = p[n]
            want = (ax == "temp") ? nt : (ax == "sal") ? ns : na
            if (want > 0 && $5 + 0 != want) {
                printf "  MISMATCH %-40s testable=%s expected=%s\n", $2, $5, want; bad++
            }
        }
        END { exit (bad > 0) }' "$1"
}

score_seed() {
    local seed="$1"
    local proj="MVP${seed}${PROJ_SUFFIX}"
    local truth_proj="${TRUTH_PROJ:-$proj}"
    local data="$PIPELINE_ROOT/data/mvp/$truth_proj"
    local out="$SWEEP_DIR/$proj"
    local log="$out/score.log"
    mkdir -p "$out"; : > "$log"

    local n_temp n_sal n_any
    n_temp=$(n_causal "$data/truth_temp.tsv")
    n_sal=$(n_causal "$data/truth_sal.tsv")
    n_any=$(n_causal "$data/truth_any.tsv")
    echo "[$proj] causal: temp=$n_temp sal=$n_sal any=$n_any" | tee -a "$log"

    local rank_out="$CONT_SWEEP/$proj/rank.tsv"
    rm -f "$out/rank.tsv"

    for i in $(seq 1 "$CELLS"); do
        local cd_host="$PARAMS_DIR/$proj/c${i}"
        local cd_cont="$CONT_PARAMS/$proj/c${i}"
        [[ -s "$cd_host/EMMAX_pvalues.tsv" ]] || { echo "[$proj] c$i not harvested, skipping" | tee -a "$log"; continue; }

        # Per-axis tables: keep SNPID/chr/pos plus the one trait column. sweep_thresholds.R
        # treats every non-structural column as a trait and unions calls across them, so the
        # single-axis question can only be asked of a single-axis table.
        for m in "${UNIVARIATE[@]}"; do
            awk -F'\t' 'NR==1{for(j=1;j<=NF;j++) c[$j]=j}
                        {print $c["SNPID"]"\t"$c["chr"]"\t"$c["pos"]"\t"$c["bio_1"]}' \
                "$cd_host/${m}_pvalues.tsv" > "$cd_host/${m}_pvalues_bio_1.tsv"
            awk -F'\t' 'NR==1{for(j=1;j<=NF;j++) c[$j]=j}
                        {print $c["SNPID"]"\t"$c["chr"]"\t"$c["pos"]"\t"$c["bio_2"]}' \
                "$cd_host/${m}_pvalues.tsv" > "$cd_host/${m}_pvalues_bio_2.tsv"
        done

        # ---- rank mode: AUC-PR per method x axis. This is what selects a method's best rung.
        for m in "${UNIVARIATE[@]}"; do
            R eval_detection.R --mode=rank --pvalues="$cd_cont/${m}_pvalues.tsv" \
                --truth="/pipeline/data/mvp/$truth_proj/truth_temp.tsv" --traits=bio_1 \
                --label="${proj}|${m}|c${i}|temp" --out="$rank_out" >> "$log" 2>&1
            if (( n_sal > 0 )); then
                R eval_detection.R --mode=rank --pvalues="$cd_cont/${m}_pvalues.tsv" \
                    --truth="/pipeline/data/mvp/$truth_proj/truth_sal.tsv" --traits=bio_2 \
                    --label="${proj}|${m}|c${i}|sal" --out="$rank_out" >> "$log" 2>&1
            fi
            R eval_detection.R --mode=rank --pvalues="$cd_cont/${m}_pvalues.tsv" \
                --truth="/pipeline/data/mvp/$truth_proj/truth_any.tsv" --traits=all \
                --label="${proj}|${m}|c${i}|any" --out="$rank_out" >> "$log" 2>&1
        done
        R eval_detection.R --mode=rank --pvalues="$cd_cont/RDA_pvalues.tsv" \
            --truth="/pipeline/data/mvp/$truth_proj/truth_any.tsv" --traits=all \
            --label="${proj}|RDA|c${i}|any" --out="$rank_out" >> "$log" 2>&1

        # ---- threshold x combine surface
        local spec_any="" spec_temp="" spec_sal=""
        for m in "${ALL_METHODS[@]}"; do spec_any+="${m}=${cd_cont}/${m}_pvalues.tsv,"; done
        for m in "${UNIVARIATE[@]}"; do
            spec_temp+="${m}=${cd_cont}/${m}_pvalues_bio_1.tsv,"
            spec_sal+="${m}=${cd_cont}/${m}_pvalues_bio_2.tsv,"
        done
        R sweep_thresholds.R --combine-window-kb="$COMBINE_WINDOW_KB" --top-max="$TOP_MAX" --methods="${spec_any%,}"  --truth="/pipeline/data/mvp/$truth_proj/truth_any.tsv" \
            --tag="${proj}_c${i}_any"  --outdir="$CONT_SWEEP/$proj" >> "$log" 2>&1
        R sweep_thresholds.R --combine-window-kb="$COMBINE_WINDOW_KB" --top-max="$TOP_MAX" --methods="${spec_temp%,}" --truth="/pipeline/data/mvp/$truth_proj/truth_temp.tsv" \
            --tag="${proj}_c${i}_temp" --outdir="$CONT_SWEEP/$proj" >> "$log" 2>&1
        R sweep_thresholds.R --combine-window-kb="$COMBINE_WINDOW_KB" --top-max="$TOP_MAX" --methods="${spec_sal%,}"  --truth="/pipeline/data/mvp/$truth_proj/truth_sal.tsv" \
            --tag="${proj}_c${i}_sal"  --outdir="$CONT_SWEEP/$proj" >> "$log" 2>&1

        if (( i == 1 )); then
            if truth_join_gate "$out/rank.tsv" "$n_temp" "$n_sal" "$n_any" >> "$log" 2>&1; then
                echo "[$proj] truth-join gate PASSED" | tee -a "$log"
            else
                echo "[$proj] FATAL: truth-join gate FAILED — causal loci missing from the" \
                     "scored SNP set, every recall number would be inflated. See $log." | tee -a "$log"
                return 1
            fi
        fi
    done

    # ---- ORACLE bundle: each method at its own AUC-PR-best rung, chosen against the truth
    # table. This is an UPPER BOUND, not a configuration any user can reach -- a user has no
    # truth table. Reported as such everywhere. Its counterpart is the PreGEA bundle, which
    # mvp_pick_bundles.R also emits and which IS achievable.
    R mvp_pick_bundles.R --seed="$seed" --cells="$CELLS" > "$out/bundles.txt" 2>> "$log"
    while IFS=$'\t' read -r kind axis truth spec; do
        [[ -z "${spec:-}" ]] && continue
        R sweep_thresholds.R --combine-window-kb="$COMBINE_WINDOW_KB" --top-max="$TOP_MAX" --methods="$spec" --truth="$truth" \
            --tag="${proj}_${kind}_${axis}" \
            --outdir="$CONT_SWEEP/$proj" >> "$log" 2>&1
    done < "$out/bundles.txt"

    echo "[$proj] SCORED $(date -Is)" | tee -a "$log"
}

echo "INFO: scoring ${#SEEDS[@]} seeds, $SEED_JOBS at a time"
for seed in "${SEEDS[@]}"; do
    while (( $(jobs -rp | wc -l) >= SEED_JOBS )); do wait -n; done
    score_seed "$seed" &
done
wait
echo "INFO: scoring finished"
