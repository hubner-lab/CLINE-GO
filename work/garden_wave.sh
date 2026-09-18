#!/usr/bin/env bash
# Step 9 garden sweep for one SS-Clines block, one 60-seed half per invocation.
#
#   work/garden_wave.sh <cohort> h1|h2
#   work/garden_wave.sh ssclines_b4 h1
#
# TWO WATCHDOGS WOULD FIGHT. mvp_garden_sweep.sh:56-66 polls host-wide `free -g` and stops
# `docker ps --filter name=gsweep_ | head -1` -- which matches the OTHER half's containers
# too. Both drivers would trip on the same reading and each kill an arbitrary container.
# So h1 carries the real ceiling and h2 gets an unreachable one: exactly one watchdog
# polices the host.
#
# MEM_PER_SEED stays 40g. A docker --memory cap reserves nothing; it only decides what
# happens on exceed (exit 137 + a stale Snakemake lock). b1 ran this stage at 40g/-c2.
set -uo pipefail

COHORT="${1:?cohort}"
HALF="${2:?h1|h2}"
# Seeds run concurrently per wave. b4 ran 16 (32 total) and took 21 h while sharing the box
# with b5's GEA sweep; the RAM guard logged ZERO threshold crossings the whole time, so that
# was under-subscribed. A block whose garden sweep runs alone can go wider.
JOBS="${3:-16}"
# snakemake -cN per seed. Intra-seed parallelism scales BADLY: b3 ran -c8 for 1.51 h/seed
# against b1's -c2 at 2.05 h/seed -- 4x the cores for 1.36x the speed. Seed-level width is
# the real lever, so keep this modest and spend the machine on JOBS instead.
SCORES="${4:-2}"
SHORT="${COHORT#ssclines_}"
export PIPELINE_ROOT=/mnt/data/eugene/ADAPTOGENE
cd "$PIPELINE_ROOT" || exit 1

SEEDS_FILE="work/seeds_${SHORT}_${HALF}.txt"
[[ -s "$SEEDS_FILE" ]] || { echo "FATAL: $SEEDS_FILE missing"; exit 1; }
SEEDS=$(paste -sd, "$SEEDS_FILE")

# SINGLE WATCHDOG. The driver's own guard (mvp_garden_sweep.sh:56-66) stops
# `docker ps --filter name=gsweep_ | head -1` -- an ARBITRARY container, and with two waves
# it can kill the other wave's work. work/ram_guard.sh is strictly better: it picks the
# LARGEST-RSS victim, logs the full candidate list first, and exempts the scoring/tables
# singletons that hold the offset12 lock. Two watchdogs on one reading also double-kill.
# So both waves get an unreachable driver ceiling and ram_guard.sh is the sole authority.
# If ram_guard.sh is not alive, refuse to launch rather than run unguarded.
pgrep -f "work/ram_guard.sh" >/dev/null || { echo "FATAL: work/ram_guard.sh is not running -- refusing to launch an unguarded garden sweep"; exit 1; }
CEIL=100000

export OUT="$PWD/benchmarks/mvp_eval/offset12"
export LOGDIR="$PWD/benchmarks/mvp_eval/offset12/sweeplogs_${COHORT}_${HALF}"
export CONFIG_SUFFIX=_sweep
export BLAS_THREADS=2
export HOST_RAM_CEILING_GB=$CEIL
mkdir -p "$LOGDIR"

echo "=== $COHORT garden wave $HALF : $(wc -l < "$SEEDS_FILE") seeds, $JOBS concurrent, snake -c$SCORES, ceiling=${CEIL}G  $(date -Is)"
benchmarks/mvp_garden_sweep.sh "$SEEDS" "$JOBS" "$SCORES" 40g 2>&1 | tee "logs_garden_${COHORT}.${HALF}.log"
echo "=== $COHORT garden wave $HALF finished  $(date -Is)"
