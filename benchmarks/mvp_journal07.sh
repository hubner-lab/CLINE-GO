#!/usr/bin/env bash
# =============================================================================
# mvp_journal07.sh -- journal 06's protocol, replicated across ALL 14 replicates
#                     in BOTH MAF arms.
#
#   benchmarks/mvp_journal07.sh configs            # write every config, both arms, ONCE
#   benchmarks/mvp_journal07.sh run    [base|m05]  # the 5-cell 11-method ladder, fanned out
#   benchmarks/mvp_journal07.sh truth              # MAF-matched truth from the REAL filtered VCFs
#   benchmarks/mvp_journal07.sh score  [base|m05]  # redundancy -> sweep -> portfolio -> transfer -> wza
#   benchmarks/mvp_journal07.sh aggregate          # the cross-seed tables
#
# WHAT IS DIFFERENT FROM mvp_journal06.sh, AND WHY EACH DIFFERENCE EXISTS
# ----------------------------------------------------------------------
# journal 06 ran ONE seed. Fanning that driver out over 13 seeds would be wrong in three
# specific ways, each of which is silent:
#
# 1. ONE MANIFEST PER ARM, WRITTEN ONCE, UP FRONT.
#    journal06.sh derives CELLS_TSV="benchmarks/mvp_sweep_cells_p11${ARM_TAG}.tsv" -- no seed
#    in the name. Concurrent per-seed invocations would each call mvp_write_sweep_configs.R
#    and clobber that single file; mvp_score_sweep.sh would then read whichever call won the
#    race. So config generation is its own stage here, one call per arm covering every seed,
#    and the run/score stages never write a manifest.
#
# 2. TOP_MAX IS PER SEED, NOT A CONSTANT 100.
#    journal 06 sec.1 states outright that the top<=100 cap is what forced the
#    single-replicate choice: on a seed with 395 causal loci, top-100 caps recall at 0.19 and
#    F1 at ~0.32 no matter which method is used, so a cross-seed comparison at a fixed cap
#    would measure the cap. Rule, applied identically to both arms of a seed:
#
#        TOP_MAX = smallest TOP_GRID value >= 2 x (testable causal, BASE arm), floor 100
#
#    Both arms of one seed share the value so the arm contrast is not also a grid contrast.
#    It is taken from the base arm because that is the larger causal count, so recall 1.0
#    stays reachable in both. On the anchor's m05 arm (43 causal) the rule returns 100 --
#    journal 06's published m05 numbers are therefore reproduced exactly, which is the
#    regression check for this whole change. On its base arm the rule returns 200 rather than
#    journal 06's 100; the anchor is RE-SCORED at 200 here (not re-run) so all 14 replicates
#    share one grid. Journal 06's tables are untouched on disk.
#
#    CONSEQUENCE, and it must be stated wherever these numbers are: absolute F1 is no longer
#    comparable BETWEEN seeds. What is comparable is the delta (combine - best solo), the
#    identity of the winning subset, and its rank.
#
# 3. THE ANCHOR IS NOT RE-RUN.
#    params07{,_m05}/MVP1232548{,M05} are symlinks into journal 06's harvest. Re-running it
#    would move the reference the whole journal is anchored on.
#
# THE TRANSFER TEST -- the reason this journal exists
# ---------------------------------------------------
# `f1_gain` in journal 06's solo_vs_combine.tsv is in-sample on BOTH sides: subset,
# assignment and window are all chosen on the same 43 causal loci the gain is then reported
# on. It cannot support "combining generalises". So the headline here is not the per-seed
# search but the TRANSFER of the anchor's winning schemes, held fixed, onto seeds that had no
# part in choosing them (benchmarks/mvp_anchor_schemes.tsv, evaluated by sweep_portfolio.R
# --schemes). The per-seed re-search still runs -- it is the per-dataset recommendation, and
# the gap between transfer and re-search IS the selection bias, measured.
#
# Docker: /snap/bin/docker cannot start on this host (setuid blocked under nosuid+NoNewPrivs),
# so every invocation goes through a nix-provided client -- see CLAUDE.local.md.
# -e USER=pipeline is required or snakemake's getpwuid() lookup aborts under --user.
# =============================================================================
set -uo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
STAGE="${1:-}"
ARM="${2:-base}"

MANIFEST="$PIPELINE_ROOT/benchmarks/mvp_seeds.tsv"
EVAL="$PIPELINE_ROOT/benchmarks/mvp_eval"
ANCHOR_SEED="1232548"

# The 13 replicates run here, plus the anchor, symlinked in at score time.
SEEDS_ALL="$(awk 'NR>1 && NF {print $1}' "$MANIFEST" | paste -sd,)"
SEEDS_RUN="$(awk -v a="$ANCHOR_SEED" 'NR>1 && NF && $1!=a {print $1}' "$MANIFEST" | paste -sd,)"

METHODS="EMMAX,LFMM,RDA,GLM,MLM,CMLM,ECMLM,SUPER,MLMM,FarmCPU,BLINK"
UNIVARIATE="EMMAX LFMM GLM MLM CMLM ECMLM SUPER MLMM FarmCPU BLINK"
ALL_METHODS="EMMAX LFMM RDA GLM MLM CMLM ECMLM SUPER MLMM FarmCPU BLINK"
SUFFIX="_p11"
CELLS="${CELLS:-5}"
CW="${CW:-0,1,2.5,5}"

# Resource envelope. This is a SHARED host (other users' jobs were holding ~300 GB and
# load ~120 when this was sized), and the user's grant is 120 cpu / 600 GB. Peak measured
# for one 11-method seed in journal 06 was 81 GB (mvp_eval/c1_memtrace.txt), so 5 concurrent
# seeds at a 110 GB cap sits inside the grant with headroom for the other tenants.
SEED_JOBS="${SEED_JOBS:-5}"
CPUS_PER_SEED="${CPUS_PER_SEED:-24}"
MEM_PER_SEED="${MEM_PER_SEED:-110g}"
SNAKE_CORES="${SNAKE_CORES:-6}"
# BLAS_THREADS must stay 16 and IDENTICAL to journals 05/06: unset, the container sees the
# host's 192 cores and OpenBLAS thrashes (40x slower, 81 GB); and thread count reorders
# floating-point summation, so mixing values inside one benchmark is a reproducibility hole.
BLAS_THREADS="${BLAS_THREADS:-16}"

case "$ARM" in
  base) ARM_TAG="";     PROJ_SUFFIX="";    ARM_MAF="" ;;
  m05)  ARM_TAG="_m05"; PROJ_SUFFIX="M05"; ARM_MAF="0.05" ;;
  *)    echo "Unknown arm '$ARM' (expected base|m05)" >&2; exit 2 ;;
esac

PARAMS_REL="benchmarks/mvp_eval/params07${ARM_TAG}"
SWEEP_REL="benchmarks/mvp_eval/sweep07${ARM_TAG}"
PARAMS="$PIPELINE_ROOT/$PARAMS_REL"
SWEEP="$PIPELINE_ROOT/$SWEEP_REL"
CONT_PARAMS="/pipeline/$PARAMS_REL"
CONT_SWEEP="/pipeline/$SWEEP_REL"
CELLS_TSV="benchmarks/mvp_sweep_cells_j07${ARM_TAG}.tsv"
TOPMAX_TSV="$PIPELINE_ROOT/benchmarks/mvp_seed_topmax.tsv"
SCHEMES="/pipeline/benchmarks/mvp_anchor_schemes.tsv"
RUNLOGS="$EVAL/runlogs07${ARM_TAG}"

DOCKER=(nix shell nixpkgs#docker-client -c docker)
R() {   # R <script> <args...>
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
        -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --cpus=16 --memory=64g -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
        Rscript "/pipeline/benchmarks/$1" "${@:2}"
}
topmax_of() { awk -v s="$1" 'NR>1 && $1==s {print $4; exit}' "$TOPMAX_TSV"; }
spec() {   # spec <proj> <cell> <method...>  -> "NAME=PATH,NAME=PATH"
    local p="$1" d="$2"; shift 2
    local out=""
    for m in "$@"; do
        [[ -s "$PARAMS/$p/$d/${m}_pvalues.tsv" ]] || continue
        out+="${m}=$CONT_PARAMS/$p/$d/${m}_pvalues.tsv,"
    done
    echo "${out%,}"
}

# ------------------------------------------------------------------------- configs
if [[ "$STAGE" == "configs" ]]; then
    echo "=== per-seed TOP_MAX ==="
    # Derived from the BASE-arm truth table: the count of causal loci that survive into the
    # scored SNP set. Written to disk rather than recomputed per call, because a cap that
    # silently differed between the run that produced a table and the run that scored it
    # would be invisible in every output.
    {
      printf 'seed\tn_causal_base\tn_causal_x2\ttop_max\n'
      for s in ${SEEDS_ALL//,/ }; do
        # truth_any.tsv is chr/pos/category; `causal` is an exact category value, not a prefix
        # (linked_neutral and background_neutral must not match).
        n=$(awk -F'\t' 'NR>1 && $3=="causal" {c++} END {print c+0}' \
              "$PIPELINE_ROOT/data/mvp/MVP${s}/truth_any.tsv")
        want=$(( n * 2 ))
        tm=100
        for g in 100 200 500 1000 2000 5000; do if (( g >= want )); then tm=$g; break; fi; done
        printf '%s\t%s\t%s\t%s\n' "$s" "$n" "$want" "$tm"
      done
    } > "$TOPMAX_TSV"
    column -t -s$'\t' "$TOPMAX_TSV"

    echo "=== configs, base arm (13 seeds) ==="
    R mvp_write_sweep_configs.R --seeds="$SEEDS_RUN" --methods="$METHODS" \
        --wza-window=1000 --suffix="$SUFFIX" \
        --manifest-out="benchmarks/mvp_sweep_cells_j07.tsv" || exit 1

    echo "=== configs, m05 arm (13 seeds) ==="
    R mvp_write_sweep_configs.R --seeds="$SEEDS_RUN" --methods="$METHODS" \
        --wza-window=1000 --suffix="$SUFFIX" \
        --project-suffix=M05 --maf=0.05 \
        --manifest-out="benchmarks/mvp_sweep_cells_j07_m05.tsv" || exit 1
    exit 0
fi

# ----------------------------------------------------------------------------- run
if [[ "$STAGE" == "run" ]]; then
    mkdir -p "$RUNLOGS"
    echo "=== ladder: arm=$ARM, seeds=$SEEDS_RUN, $SEED_JOBS at a time ==="
    SWEEP_METHODS="$ALL_METHODS" CFG_SUFFIX="$SUFFIX" PARAMS_DIR="$PARAMS" \
        PROJ_SUFFIX="$PROJ_SUFFIX" HARVEST_WZA=1 RUNLOG_DIR="$RUNLOGS" \
        CPUS_PER_SEED="$CPUS_PER_SEED" MEM_PER_SEED="$MEM_PER_SEED" \
        SNAKE_CORES="$SNAKE_CORES" BLAS_THREADS="$BLAS_THREADS" \
        "$PIPELINE_ROOT/benchmarks/mvp_run_sweep.sh" "$SEEDS_RUN" "$SEED_JOBS" "$CELLS"
    rc=$?
    echo "=== ladder arm=$ARM exit=$rc ==="
    exit $rc
fi

# --------------------------------------------------------------------------- truth
# MAF-matched truth for the m05 arm, taken from the ACTUAL filtered VCF the pipeline wrote
# (not from loci_freq): plink's --maf keeps MAF >= threshold, and on the anchor exactly 18
# SNPs sit on the 0.05 boundary, one of them causal. Run AFTER the m05 ladder, because the
# filtered VCF is what mode=processing produces.
if [[ "$STAGE" == "truth" ]]; then
    fail=0
    for s in ${SEEDS_ALL//,/ }; do
        v=$(ls "$PIPELINE_ROOT/MVP${s}M05_results/_work/maf0.05_"*"/MVP${s}.vcf" 2>/dev/null | head -1)
        if [[ -z "$v" ]]; then
            echo "SKIP $s: no filtered VCF yet (m05 ladder not run?)"; fail=1; continue
        fi
        R mvp_filter_truth_maf.R --seed="$s" --maf=0.05 --out-suffix=M05 \
            --vcf="/pipeline/${v#$PIPELINE_ROOT/}" || fail=1
    done
    exit $fail
fi

# --------------------------------------------------------------------------- score
if [[ "$STAGE" == "score" ]]; then
    mkdir -p "$PARAMS" "$SWEEP"
    # The anchor joins the score stage without being re-run.
    [[ -e "$PARAMS/MVP${ANCHOR_SEED}${PROJ_SUFFIX}" ]] || \
        ln -s "../params06${ARM_TAG}/MVP${ANCHOR_SEED}${PROJ_SUFFIX}" \
              "$PARAMS/MVP${ANCHOR_SEED}${PROJ_SUFFIX}"

    for s in ${SEEDS_ALL//,/ }; do
        proj="MVP${s}${PROJ_SUFFIX}"
        data="/pipeline/data/mvp/$proj"
        tm="$(topmax_of "$s")"; tm="${tm:-100}"
        [[ -d "$PARAMS/$proj/c3" ]] || { echo "--- $proj not harvested, skipping"; continue; }
        echo "=== $proj  (top_max=$tm) ==="

        R mvp_method_redundancy.R --methods="$(spec "$proj" c3 $ALL_METHODS)" \
            --truth="$data/truth_any.tsv" --tag="${proj}_c3" \
            --outdir="/pipeline/benchmarks/mvp_eval/redundancy07${ARM_TAG}"

        SWEEP_UNIVARIATE="$UNIVARIATE" SWEEP_ALL_METHODS="$ALL_METHODS" \
            PARAMS_DIR="$PARAMS" SWEEP_DIR="$SWEEP" \
            CONT_PARAMS="$CONT_PARAMS" CONT_SWEEP="$CONT_SWEEP" \
            SWEEP_CELLS_TSV="$CELLS_TSV" PROJ_SUFFIX="$PROJ_SUFFIX" TRUTH_PROJ="$proj" \
            COMBINE_WINDOW_KB="$CW" TOP_MAX="$tm" \
            "$PIPELINE_ROOT/benchmarks/mvp_score_sweep.sh" "$s" 1 "$CELLS"

        # --rank-source pinned to c3 for every cell, so all five rungs search one method pool
        # (journal 05 sec.10 -- re-ranking per cell silently gives each cell a different pool).
        rank="$CONT_SWEEP/$proj/${proj}_c3_any_threshold_sweep.tsv"
        for i in $(seq 1 "$CELLS"); do
            sp="$(spec "$proj" "c$i" $ALL_METHODS)"
            [[ -n "$sp" ]] || { echo "  c$i not harvested, skipping"; continue; }
            R sweep_portfolio.R --methods="$sp" --truth="$data/truth_any.tsv" \
                --tag="${proj}_c${i}_any" \
                --sweep="$CONT_SWEEP/$proj/${proj}_c${i}_any_threshold_sweep.tsv" \
                --rank-source="$rank" --combine-window-kb="$CW" --top-max="$tm" \
                --schemes="$SCHEMES" \
                --outdir="/pipeline/benchmarks/mvp_eval/portfolio07${ARM_TAG}"
        done

        wspec=""
        for m in $ALL_METHODS; do
            [[ -s "$PARAMS/$proj/c3/${m}_wza.tsv" ]] || continue
            wspec+="${m}=$CONT_PARAMS/$proj/c3/${m}_wza.tsv,"
        done
        # WZA is scored on bonf/qval only. Journal 06 sec.13.3: including `top` in the window
        # grid forces a fixed number of calls regardless of evidence, which contradicts the
        # point of windowing and inflated that journal's WZA arm.
        [[ -n "$wspec" ]] && R mvp_wza_gate.R --wza="${wspec%,}" \
            --truth="$data/truth_any.tsv" --loci-freq="$data/loci_freq.tsv" \
            --tag="${proj}_c3" --top-max=0 \
            --outdir="/pipeline/benchmarks/mvp_eval/wza07${ARM_TAG}"
    done
    exit 0
fi

# ----------------------------------------------------------------------- aggregate
if [[ "$STAGE" == "aggregate" ]]; then
    R mvp_crossseed_report.R --seeds="$SEEDS_ALL" --cells="$CELLS" \
        --outdir="/pipeline/benchmarks/mvp_eval/report07"
    exit $?
fi

echo "Usage: $0 {configs|run|truth|score|aggregate} [base|m05]" >&2
exit 2
