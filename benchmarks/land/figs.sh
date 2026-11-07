#!/usr/bin/env bash
# figs.sh <land> -- Phase 8 display items per landscape (user, 2026-11-06: Figure 2 + S2 only):
#   dtable  detection_table.tsv from the landscape's remeasure arms (scratch/detection_table_land.R)
#   figs    Figure 2 (A1 + R1), A1/D5/D6, R1/R2, S2 (B4) + numbers.tsv (scratch/fig2_s2_land.R)
# Both copies reproduce SS-Clines' published detection_table, 7 SVGs and 807 numbers byte for byte
# (tie check 2026-11-06, work/mvp_land/ab/land_*). Fixed SS-Clines numbers are soft asserts, logged
# to results/soft_asserts_<land>.tsv.
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE
land=${1:?land}; L=work/mvp_land; SC=/pipeline/$L/scratch; EV=/pipeline/benchmarks/mvp_eval
SEEDF=$L/seeds_$land.txt; OFF=offset15_$land; FIG=figures_land_${land}_$OFF; GAL=figures_land_${land}_gallery
RM=remeasure600_$land; DT=figures_land_${land}_dtable; MS=figures_land_${land}_ms
ARMENV=(-e PIPELINE_ROOT=/pipeline -e MVP_SEEDS_FILE=/pipeline/$SEEDF -e MVP_N_EXPECT=600 -e MVP_N_PER_ARCH=200 -e MVP_SUBTITLE=0)
LOGS=$L/logs/figs_$land; mkdir -p $LOGS; SOFT=$L/results/soft_asserts_$land.tsv; : > $SOFT
die() { echo "FIGS FAIL [$land $1] $(date -Is)"; exit 1; }
run() { local name=$1; shift; echo "--- $land $name $(date -Is)"
        $L/dk.sh figs "mland_figs_${land}_$name" "$@" > "$LOGS/$name.log" 2>&1 || die "$name (log $LOGS/$name.log)"; }
run dtable --memory=16g --cpus=2 "${ARMENV[@]}" -e REM_ROOT=$EV/$RM -e FIG_OUT=$EV/$DT -- Rscript $SC/detection_table_land.R
run figs --memory=16g --cpus=4 "${ARMENV[@]}" -e OPENBLAS_NUM_THREADS=2 -e SOFT_LOG=/pipeline/$SOFT \
    -e GALLERY_DIR=$EV/$GAL -e OFFSET_DIR=$OFF -e RANK_DIR=$EV/$FIG/rank -e RDA_FIG_DIR=$EV/$DT \
    -e REMEASURE_DIR=$EV/$RM/rdaunc_rda2x_ldb -e FIG_OUT=$EV/$MS -- Rscript $SC/fig2_s2_land.R
echo "    soft asserts: $(cut -f1 $SOFT | sort | uniq -c | tr '\n' ';')"
echo "=== figs $land done $(date -Is)"
