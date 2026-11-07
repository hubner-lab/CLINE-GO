#!/usr/bin/env bash
# mvp_land_empty_sets.sh [land ...] -- Job C of the landscapes pack, empty sets at significance thresholds:
# the remeasure harness at its default GRID (bonf, qval, top, custom), each replicate's own LD window,
# on the p-value tables the _ldb arms read, into landscape_ms/empty_sets/<land>/. The _ldb arms ran with
# TOP_SHARES and hold only top rungs, which is why this exists. mvp_land_empty_sets.R reads the result.
#   ssclines          committed benchmarks/mvp_remeasure_600.R, params_rdaunc, p4/seed_bp_out600.tsv
#   ssmtn, estclines  work/mvp_land/scratch/remeasure_600_land.R, params_land_<land>_rdaunc, seed_ld_bp_<land>.tsv
# Settings of each landscape's rdaunc_rda2x_ldb provenance.tsv, minus TOP_SHARES.
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE
EV=/pipeline/benchmarks/mvp_eval; OUT=$EV/landscape_ms/empty_sets; LOGS=benchmarks/mvp_eval/landscape_ms/logs
TAGS=ssclines_nvar_mvar,ssclines_ncline_ns,ssclines_nequal_mconst,ssclines_ncline_ctredge,ssclines_nequal_mbreaks
COMMON=(--cpus=8 --memory=64g -e PIPELINE_ROOT=/pipeline -e MVP_N_EXPECT=600 -e MVP_N_PER_ARCH=200 -e MVP_SUBTITLE=0
        -e ARM=rdaunc_thresholds_ldb -e A3_METHODS=LFMM,EMMAX -e SHARE_MULT=RDA=2 -e MVP_SEED_WINDOW=1 -e CELL=c1 -e NCORES=8)
LANDS=("$@"); [[ ${#LANDS[@]} -eq 0 ]] && LANDS=(ssclines ssmtn estclines)
for land in "${LANDS[@]}"; do
    echo "--- $land $(date -Is)"
    mkdir -p "benchmarks/mvp_eval/landscape_ms/empty_sets/$land"
    if [[ $land == ssclines ]]; then
        benchmarks/mvp_land_dk.sh empty_$land "${COMMON[@]}" -e MVP_ARM=primary_ssclines -e MVP_ADDED=$TAGS \
            -e PARAMS_DIR=$EV/params_rdaunc -e MVP_LINKED_BP=/pipeline/work/ld_window_lab/p4/seed_bp_out600.tsv \
            -e OUT_DIR=$OUT/$land -- Rscript benchmarks/mvp_remeasure_600.R > $LOGS/empty_$land.log 2>&1
    else
        benchmarks/mvp_land_dk.sh empty_$land "${COMMON[@]}" -e MVP_SEEDS_FILE=/pipeline/work/mvp_land/seeds_$land.txt \
            -e PARAMS_DIR=$EV/params_land_${land}_rdaunc -e A3_REF=$EV/detection600_$land/aucpr_per_seed.tsv \
            -e MVP_LINKED_BP=/pipeline/work/mvp_land/results/seed_ld_bp_$land.tsv \
            -e OUT_DIR=$OUT/$land -- Rscript work/mvp_land/scratch/remeasure_600_land.R > $LOGS/empty_$land.log 2>&1
    fi
    rc=$?
    echo "    exit $rc $(date -Is)"
    [[ $rc -eq 0 ]] || { echo "EMPTY SETS FAIL $land (log $LOGS/empty_$land.log)"; exit 1; }
done
echo "=== empty sets done $(date -Is)"
