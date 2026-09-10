#!/usr/bin/env bash
# =============================================================================
# mvp_journal06.sh -- the journal-06 pipeline, end to end, on ONE replicate.
#
#   benchmarks/mvp_journal06.sh [run|score|all] [SEED]        # ARM=base (default)
#   ARM=m05 benchmarks/mvp_journal06.sh [run|score|all] [SEED]  # the Filter.maf 0.05 arm
#
#   run    generate configs + run the 5-cell 11-method ladder (+ the EMMAX kinship arm)
#   score  everything post-hoc: redundancy -> threshold sweep -> portfolio -> WZA -> profile -> report
#   all    both, in order
#
# WHY A DRIVER AND NOT SIX HAND-TYPED DOCKER LINES
# ------------------------------------------------
# Six invocations that must agree on method lists, truth-table pairings, combine windows and
# the top-N cap is exactly the surface where a benchmark silently scores the wrong thing. The
# scoring contract (benchmarks/MVP_README.md:41-62) is encoded ONCE here, and the journal's
# command chunks quote this file rather than re-deriving it.
#
# ORDERING IS LOAD-BEARING, not convenience:
#   redundancy BEFORE portfolio  -- the family representatives decide which methods the subset
#                                   search is allowed to draw from
#   sweep      BEFORE portfolio  -- sweep_portfolio.R reads the per-method curves to pick each
#                                   method's candidate operating points
#   c3 sweep   BEFORE every other cell's portfolio run -- the method pool is ranked ONCE, on
#                                   the default cell, and reused, or the configuration groups
#                                   stop being comparable (journal 05 sec.10)
#
# Docker: /snap/bin/docker cannot start on this host (setuid blocked under nosuid+NoNewPrivs),
# so every invocation goes through a nix-provided client -- see CLAUDE.local.md.
# =============================================================================
set -uo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
STAGE="${1:-all}"
SEED="${2:-1232548}"
CELLS="${CELLS:-5}"

# The 11 methods. GAPIT models all share one n_pcs per cell by construction
# (_gapit_shared_npcs() in gea.smk hard-errors otherwise), which is what the cell design
# already does. RDA is the only multivariate one and is scored on the `any` axis only.
METHODS="EMMAX,LFMM,RDA,GLM,MLM,CMLM,ECMLM,SUPER,MLMM,FarmCPU,BLINK"
UNIVARIATE="EMMAX LFMM GLM MLM CMLM ECMLM SUPER MLMM FarmCPU BLINK"
ALL_METHODS="EMMAX LFMM RDA GLM MLM CMLM ECMLM SUPER MLMM FarmCPU BLINK"

SUFFIX="_p11"
EVAL="$PIPELINE_ROOT/benchmarks/mvp_eval"

# ---- ARMS -------------------------------------------------------------------------------
# ARM=base   Filter.maf 0.01 (the deposit's own no-op cut), project MVP{seed}
# ARM=m05    Filter.maf 0.05, project MVP{seed}M05, MAF-MATCHED truth tables
#
# The two arms are separate PROJECTS, not two configs of one project: {PROJECT}_results/ is
# keyed on project_name alone and the module tables carry no filter tag, so one project cannot
# hold two MAF values without overwriting itself (CLAUDE.md documents this hazard). Separate
# project names also mean the arms can run CONCURRENTLY -- separate Snakemake locks.
#
# The m05 arm scores against data/mvp/MVP{seed}M05/truth_*.tsv, produced by
# mvp_filter_truth_maf.R. Scoring it against the unfiltered truth would divide recall by loci
# the VCF no longer contains, which is exactly the defect that killed the Laruson benchmark.
ARM="${ARM:-base}"
case "$ARM" in
  base) ARM_TAG="";     PROJ_SUFFIX="";    ARM_MAF="" ;;
  m05)  ARM_TAG="_m05"; PROJ_SUFFIX="M05"; ARM_MAF="0.05" ;;
  *)    echo "Unknown ARM=$ARM (expected base|m05)" >&2; exit 2 ;;
esac
PROJ="MVP${SEED}${PROJ_SUFFIX}"

# Host path and its container-side twin are derived from ONE name each. Writing the two
# independently is how a driver ends up reading journal 05's tables while writing journal
# 06's -- silent, and nothing downstream would flag it.
PARAMS_REL="benchmarks/mvp_eval/params06${ARM_TAG}"
SWEEP_REL="benchmarks/mvp_eval/sweep06${ARM_TAG}"
PARAMS="$PIPELINE_ROOT/$PARAMS_REL"
SWEEP="$PIPELINE_ROOT/$SWEEP_REL"
CONT_PARAMS="/pipeline/$PARAMS_REL"
CONT_SWEEP="/pipeline/$SWEEP_REL"
DATA="/pipeline/data/mvp/$PROJ"
CELLS_TSV="benchmarks/mvp_sweep_cells_p11${ARM_TAG}.tsv"

# top-N is capped at 100 by instruction. On this replicate the cap is free: 70 causal loci,
# so top-100 can still reach recall 1.0. It would censor the polygenic seeds -- which is why
# the single-replicate choice and the cap go together.
TOP_MAX="${TOP_MAX:-100}"
# Agreement windows, calibrated on this replicate's own truth table: linked_neutral sits a
# median 2 351 bp from its nearest causal locus (p75 4 496, p90 6 819); 5 kb is the config's
# GEA.snp_clumping_distance, i.e. what the pipeline itself uses when combining.
CW="${CW:-0,1,2.5,5}"

DOCKER=(nix shell nixpkgs#docker-client -c docker)
R() {   # R <script> <args...>
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
        -e OPENBLAS_NUM_THREADS=16 -e OMP_NUM_THREADS=16 \
        --cpus=16 --memory=64g -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
        Rscript "/pipeline/benchmarks/$1" "${@:2}"
}
spec() {   # spec <cell_dir> <suffix> <method...>  -> "NAME=PATH,NAME=PATH"
    local d="$1" sfx="$2"; shift 2
    local out=""
    for m in "$@"; do
        [[ -s "$PARAMS/$PROJ/$d/${m}_pvalues${sfx}.tsv" ]] || continue
        out+="${m}=$CONT_PARAMS/$PROJ/$d/${m}_pvalues${sfx}.tsv,"
    done
    echo "${out%,}"
}

# ---------------------------------------------------------------------- run
if [[ "$STAGE" == "run" || "$STAGE" == "all" ]]; then
    echo "=== [0/3] arm setup ==="
    if [[ -n "$ARM_MAF" ]]; then
        # MAF-matched truth FIRST: nothing downstream is meaningful without it.
        #
        # On the first pass the filtered VCF does not exist yet, so the truth is derived from
        # loci_freq (MAF >= threshold, matching plink's inclusive `--maf`). Once mode=processing
        # has produced the real filtered VCF, it becomes the authority -- so this step is re-run
        # AFTER the ladder's upstream stage below, and the second call is the one that counts.
        FILT_VCF="$PIPELINE_ROOT/${PROJ}_results/_work/maf${ARM_MAF}_miss0.1_smiss0.5/MVP${SEED}.vcf"
        R mvp_filter_truth_maf.R --seed="$SEED" --maf="$ARM_MAF" \
            --out-suffix="$PROJ_SUFFIX" \
            ${FILT_VCF:+$([[ -s "$FILT_VCF" ]] && echo "--vcf=/pipeline/${PROJ}_results/_work/maf${ARM_MAF}_miss0.1_smiss0.5/MVP${SEED}.vcf")} \
            || exit 1
    fi

    echo "=== [1/3] configs ==="
    R mvp_write_sweep_configs.R --seeds="$SEED" --methods="$METHODS" \
        --wza-window=1000 --suffix="$SUFFIX" \
        ${PROJ_SUFFIX:+--project-suffix="$PROJ_SUFFIX"} ${ARM_MAF:+--maf="$ARM_MAF"} \
        --manifest-out="$CELLS_TSV" || exit 1

    echo "=== [2/3] 5-cell ladder, 11 methods ==="
    SWEEP_METHODS="$ALL_METHODS" CFG_SUFFIX="$SUFFIX" PARAMS_DIR="$PARAMS" \
        PROJ_SUFFIX="$PROJ_SUFFIX" HARVEST_WZA=1 RUNLOG_DIR="$EVAL/runlogs06" \
        CPUS_PER_SEED=48 MEM_PER_SEED=200g SNAKE_CORES=8 \
        "$PIPELINE_ROOT/benchmarks/mvp_run_sweep.sh" "$SEED" 1 "$CELLS" || exit 1

    echo "=== [3/3] EMMAX kinship arm (IBS) ==="
    # Kinship is orthogonal to the MAF question; run it once, on the base arm only.
    if [[ "$ARM" == "base" ]]; then
    # Separate configs, separate harvest tree: the kinship variant writes to the SAME
    # {method}_pvalues_K{k}.tsv path as the BN run -- kinship is not in the filename any more
    # than the swept rung is -- so the two arms cannot share a params directory.
    R mvp_write_sweep_configs.R --seeds="$SEED" --methods=EMMAX --kinship=IBS \
        --suffix="${SUFFIX}ibs" --manifest-out=benchmarks/mvp_sweep_cells_p11ibs.tsv || exit 1
    SWEEP_METHODS="EMMAX" CFG_SUFFIX="${SUFFIX}ibs" PARAMS_DIR="$EVAL/params06_ibs" \
        RUNLOG_DIR="$EVAL/runlogs06" CPUS_PER_SEED=32 MEM_PER_SEED=120g SNAKE_CORES=4 \
        "$PIPELINE_ROOT/benchmarks/mvp_run_sweep.sh" "$SEED" 1 "$CELLS" || exit 1
    fi
fi

# -------------------------------------------------------------------- score
if [[ "$STAGE" == "score" || "$STAGE" == "all" ]]; then
    mkdir -p "$SWEEP/$PROJ"

    echo "=== redundancy (cell c3, the default rung) ==="
    R mvp_method_redundancy.R --methods="$(spec c3 '' $ALL_METHODS)" \
        --truth="$DATA/truth_any.tsv" --tag="${PROJ}_c3" \
        --outdir="/pipeline/benchmarks/mvp_eval/redundancy${ARM_TAG}"

    echo "=== threshold sweep + rank, all cells x 3 axes ==="
    # PARAMS_DIR/SWEEP_DIR are host paths (what this script stats); CONT_* are the
    # container-side twins the R scripts see through the bind mount. Both must move together.
    SWEEP_UNIVARIATE="$UNIVARIATE" SWEEP_ALL_METHODS="$ALL_METHODS" \
        PARAMS_DIR="$PARAMS" SWEEP_DIR="$SWEEP" \
        CONT_PARAMS="$CONT_PARAMS" CONT_SWEEP="$CONT_SWEEP" \
        SWEEP_CELLS_TSV="$CELLS_TSV" PROJ_SUFFIX="$PROJ_SUFFIX" TRUTH_PROJ="$PROJ" \
        COMBINE_WINDOW_KB="$CW" TOP_MAX="$TOP_MAX" \
        "$PIPELINE_ROOT/benchmarks/mvp_score_sweep.sh" "$SEED" 1 "$CELLS"

    echo "=== portfolio: per-method operating points x subsets x windows ==="
    # --rank-source is pinned to c3 for EVERY cell so all cells search the same method pool.
    RANK_SRC="$CONT_SWEEP/$PROJ/${PROJ}_c3_any_threshold_sweep.tsv"
    for i in $(seq 1 "$CELLS"); do
        s="$(spec "c$i" '' $ALL_METHODS)"
        [[ -n "$s" ]] || { echo "  c$i not harvested, skipping"; continue; }
        R sweep_portfolio.R --methods="$s" --truth="$DATA/truth_any.tsv" \
            --tag="${PROJ}_c${i}_any" \
            --sweep="$CONT_SWEEP/$PROJ/${PROJ}_c${i}_any_threshold_sweep.tsv" \
            --rank-source="$RANK_SRC" --combine-window-kb="$CW" --top-max="$TOP_MAX" \
            --outdir="/pipeline/benchmarks/mvp_eval/portfolio${ARM_TAG}"
    done

    echo "=== WZA viability gate + window-level scoring (c3) ==="
    wspec=""
    for m in $ALL_METHODS; do
        [[ -s "$PARAMS/$PROJ/c3/${m}_wza.tsv" ]] || continue
        wspec+="${m}=$CONT_PARAMS/$PROJ/c3/${m}_wza.tsv,"
    done
    if [[ -n "$wspec" ]]; then
        R mvp_wza_gate.R --wza="${wspec%,}" --truth="$DATA/truth_any.tsv" \
            --loci-freq="$DATA/loci_freq.tsv" --tag="${PROJ}_c3" --top-max="$TOP_MAX" \
            --outdir="/pipeline/benchmarks/mvp_eval/wza${ARM_TAG}"
    else
        echo "  no WZA tables harvested -- set HARVEST_WZA=1 on the run stage"
    fi

    echo "=== detection profiling (c3, every method at top_100) ==="
    pts=""; for m in $ALL_METHODS; do pts+="${m}:top_100,"; done
    R mvp_profile_detections.R --methods="$(spec c3 '' $ALL_METHODS)" \
        --truth="$DATA/truth_any.tsv" --loci-freq="$DATA/loci_freq.tsv" \
        --muts="/pipeline/data/mvp/raw/mutations/${SEED}_Rout_muts_full.txt.gz" \
        --tag="${PROJ}_c3" --points="${pts%,}" --combine-window-kb=5 \
        --outdir="/pipeline/benchmarks/mvp_eval/profile${ARM_TAG}"

    echo "=== aggregate ==="
    R mvp_portfolio_report.R --proj="$PROJ" --arm-tag="$ARM_TAG"
fi

echo "DONE $(date -Is)"
