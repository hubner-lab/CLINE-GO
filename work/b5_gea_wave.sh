#!/usr/bin/env bash
# Stage A -- b5 (ssclines_nequal_mbreaks) step 5 GEA sweep, one 60-seed half per invocation.
#
# Two waves on DISJOINT halves so the second is a measured decision, not a guess:
# the halves touch different MVP{seed}_results/ trees, hence different Snakemake locks,
# and mvp_run_sweep.sh writes only per-seed files (RUNLOG_DIR/{proj}.log,
# PARAMS_DIR/{proj}/cN) -- no shared state between the two driver processes.
#
#   work/b5_gea_wave.sh h1
#   work/b5_gea_wave.sh h2      # only after the wave-1 RSS reading
#
# PIPELINE_ROOT is MANDATORY: mvp_run_sweep.sh:36 defaults it to /mnt/data/eugene/CLINE-GO,
# which does not exist on this host (killed logs_gea_ssclines_b3.run2_wrongroot.log).
# BLAS_THREADS is MANDATORY: containers see the host's 192 cores, and unset OpenBLAS threads
# to all of them -- ~40x slower, ~80 GB per job. Pinned to 2 for the whole corpus.
set -uo pipefail

HALF="${1:?usage: b5_gea_wave.sh h1|h2}"
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

SEEDS_FILE="work/seeds_b5_${HALF}.txt"
[[ -s "$SEEDS_FILE" ]] || { echo "FATAL: $SEEDS_FILE missing or empty"; exit 1; }
SEEDS=$(paste -sd, "$SEEDS_FILE")
N=$(wc -l < "$SEEDS_FILE")

# 14 seeds x 4 cpu = 56 cores per wave; 112 cores when both waves run.
export CPUS_PER_SEED=4
export SNAKE_CORES=2
export BLAS_THREADS=2
export MEM_PER_SEED=40g

echo "=== b5 GEA wave $HALF : $N seeds, 14 concurrent, $(date -Is)"
benchmarks/mvp_run_sweep.sh "$SEEDS" 14 1 2>&1 | tee "logs_gea_ssclines_b5.${HALF}.log"
echo "=== b5 GEA wave $HALF finished rc=$? $(date -Is)"
