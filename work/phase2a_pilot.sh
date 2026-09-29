#!/usr/bin/env bash
# SS-Clines re-analysis, Phase 2a step 7 -- one-fit RDA pilot on 3 seeds.
# Plan: ~/.claude/plans/we-are-now-starting-radiant-pillow.md
#
# Seeds, one per genic level, k_best >= 3 (so condition_pcs > 0 and the removed second fit
# actually ran in the frozen pmax arm):
#   1231418  mod-polygenic     ssclines_nvar_mvar       K=3   56 causal  (Phase 0 pilot seed)
#   1231310  oligogenic        ssclines_ncline_ns       K=4    7 causal
#   1231131  highly-polygenic  ssclines_nequal_mbreaks  K=4  747 causal
#
# Already done by hand before this script (commands in the Phase 2a dossier bullets):
#   pmax-era RDA side tables of 1231418 stashed to onefit_pilot/pmax_diag_1231418/;
#   data sha256 before/after reconvert (onefit_pilot/data_sha256_{before,after}.txt);
#   convert_mvp_all.sh 3 "<seeds>"; mvp_write_sweep_configs.R --cells=1
#   --methods=EMMAX,LFMM,RDA --manifest-out=onefit_pilot/cells.tsv.
#
# Resources match the frozen arm's recorded invocation for all three seeds
# (--cpus=4 --memory=40g, snakemake -c2, image cline-go:latest 93bd6025); BLAS_THREADS=2.
# RUNLOG_DIR is separate: run_seed() truncates its log at seed start, and the default dir
# holds the frozen arm's runlogs. PARAMS_DIR=params_onefit: params/ is frozen (a-w).
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE

export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
export IMAGE=cline-go:latest
export PARAMS_DIR=$PIPELINE_ROOT/benchmarks/mvp_eval/params_onefit
export RUNLOG_DIR=$PIPELINE_ROOT/benchmarks/mvp_eval/onefit_pilot/runlogs
export SWEEP_METHODS="EMMAX LFMM RDA"
export SKIP_PREGEA=1
export HARVEST_DIAG=1
export BLAS_THREADS=2
export CPUS_PER_SEED=4
export MEM_PER_SEED=40g
export SNAKE_CORES=2

benchmarks/mvp_run_sweep.sh "1231418,1231310,1231131" 3 1
