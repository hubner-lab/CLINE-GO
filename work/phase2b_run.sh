#!/usr/bin/env bash
# SS-Clines re-analysis, Phase 2b -- one-fit RDA full re-run, one half of the pooled seed list.
# Plan: ~/.claude/plans/we-are-now-starting-radiant-pillow.md
#
#   work/phase2b_run.sh A 22 40g      # half A, 22 seeds at a time, 40g cap per seed
#   work/phase2b_run.sh B <jobs> <mem>  # half B, sized later from the measured peaks
#
# Seed lists: benchmarks/mvp_eval/onefit_run/seeds_{A,B}.csv -- the 600 SS-Clines primaries
# (mvp_arm.R, MVP_ARM=primary_ssclines, 5 blocks) minus the 3 Phase 2a pilot seeds (already
# harvested into params_onefit/), sorted by n_snps descending and interleaved, so each half
# runs longest-first with the same size mix. Pooled across blocks by user decision 2026-09-29.
#
# Already done before this script: reconversion of all 597 (data + truth byte-identical, only
# conversion_provenance.tsv converted_utc differs), c1 configs --methods=EMMAX,LFMM,RDA
# --cells=1 (differ from the old ones only by the dropped BLINK block).
#
# Per-seed settings are IDENTICAL to the frozen arm and the pilot (--cpus=4, -c2, BLAS 2):
# thread counts reach sNMF and BLAS, so changing them could break LFMM/EMMAX byte-identity
# with params/, which Phase 3's before/after design depends on. Speed comes from concurrency
# only. MEM_PER_SEED is a cgroup cap, not a reservation, and does not change results.
set -uo pipefail
HALF="${1:?half A|B}"; JOBS="${2:?seed jobs}"; MEM="${3:-40g}"
cd /mnt/data/eugene/ADAPTOGENE

export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
export IMAGE=cline-go:latest
export PARAMS_DIR=$PIPELINE_ROOT/benchmarks/mvp_eval/params_onefit
export RUNLOG_DIR=$PIPELINE_ROOT/benchmarks/mvp_eval/runlogs_onefit
export SWEEP_METHODS="EMMAX LFMM RDA"
export SKIP_PREGEA=1
export HARVEST_DIAG=1
export BLAS_THREADS=2
export CPUS_PER_SEED=4
export SNAKE_CORES=2
export MEM_PER_SEED="$MEM"

SEEDS=$(cat "benchmarks/mvp_eval/onefit_run/seeds_${HALF}.csv")
echo "phase2b half=$HALF jobs=$JOBS mem=$MEM start $(date -Is)"
benchmarks/mvp_run_sweep.sh "$SEEDS" "$JOBS" 1
echo "phase2b half=$HALF exit=$? end $(date -Is)"
