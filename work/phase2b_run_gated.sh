#!/usr/bin/env bash
# SS-Clines re-analysis, Phase 2b -- gated full re-run (replaces the A/B halves, 2026-09-29).
# Plan: ~/.claude/plans/we-are-now-starting-radiant-pillow.md
#
#   work/phase2b_run_gated.sh benchmarks/mvp_eval/onefit_run/seeds_C.csv
#
# Why a gate: mode=gea holds ~25-29 GB for ~18 of its ~48 min (the rda() fit) while every
# upstream mode stays under ~3 GB, and seeds launched together reach rda() together (half A's
# 22 seeds peaked at 491 GB at once). mvp_run_sweep.sh's gea gate (GEA_SLOTS, GEA_STAGGER_S)
# caps how many seeds are in gea at once and spaces their entry, so upstream runs wide on the
# CPUs and gea peaks stop coinciding. Live controls, one integer per file, re-read every check:
#   benchmarks/mvp_eval/params_onefit/.gea_gate/gea_slots      max seeds in gea
#   benchmarks/mvp_eval/params_onefit/.gea_gate/gea_stagger_s  min seconds between gea entries
#   benchmarks/mvp_eval/params_onefit/.gea_gate/seed_jobs      max seeds in flight
# The 30 seeds already in flight when A/B were stopped carry hand-written tickets there and
# are dropped from the count 180 s after their last container exits.
#
# Per-seed settings are identical to the frozen arm and the pilot (--cpus=4, -c2, BLAS 2), so
# LFMM/EMMAX stay byte-comparable with params/. Only concurrency and ordering change.
set -uo pipefail
LIST="${1:?seed list csv}"
cd /mnt/data/eugene/ADAPTOGENE

export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
export IMAGE=cline-go:latest
export PARAMS_DIR=$PIPELINE_ROOT/benchmarks/mvp_eval/params_onefit
export RUNLOG_DIR=$PIPELINE_ROOT/benchmarks/mvp_eval/runlogs_onefit
export GATE_DIR=$PARAMS_DIR/.gea_gate
export SWEEP_METHODS="EMMAX LFMM RDA"
export SKIP_PREGEA=1
export HARVEST_DIAG=1
export BLAS_THREADS=2
export CPUS_PER_SEED=4
export SNAKE_CORES=2
export MEM_PER_SEED=40g
export GEA_SLOTS=30
export GEA_STAGGER_S=60

echo "phase2b gated list=$LIST start $(date -Is)"
benchmarks/mvp_run_sweep.sh "$(cat "$LIST")" 20 1
echo "phase2b gated exit=$? end $(date -Is)"
