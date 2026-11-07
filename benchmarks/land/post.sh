#!/usr/bin/env bash
# post.sh <land> [step...] -- Phase 8 of plan nifty-prancing-grove: per-landscape scoring, the steps of
# work/ld_window_lab/p45_post.sh plus the p-value harvest, the landscape AUC-PR reference and both
# remeasure arms, which SS-Clines ran separately. Committed scripts where they are landscape-neutral,
# the copies in work/mvp_land/scratch/ (make_scratch_land.py) where they are not.
# Steps (default all, in order): complete lookup score panels detb harvest aucref remeasure oracle
#   gallery rank reliability.  Stops at the first failure.
# The slice is MVP_SEEDS_FILE=work/mvp_land/seeds_<land>.txt (benchmarks/mvp_arm.R), never MVP_ARM.
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE
land=${1:?land}; shift
[[ $land == ssmtn || $land == estclines ]] || { echo "unknown land $land"; exit 2; }
L=work/mvp_land; SC=/pipeline/$L/scratch; EV=/pipeline/benchmarks/mvp_eval; EVH=benchmarks/mvp_eval
SEEDF=$L/seeds_$land.txt; SEEDS=$(paste -sd, $SEEDF); OFF=offset15_$land
LDBP=$L/results/seed_ld_bp_$land.tsv
FIG=figures_land_${land}_$OFF; GAL=figures_land_${land}_gallery
RM=remeasure600_$land; PU=params_land_${land}_rdaunc; PO=params_land_${land}_onefit
ARMENV=(-e PIPELINE_ROOT=/pipeline -e MVP_SEEDS_FILE=/pipeline/$SEEDF -e MVP_N_EXPECT=600 -e MVP_N_PER_ARCH=200 -e MVP_SUBTITLE=0)
LOGS=$L/logs/post_$land; mkdir -p "$LOGS" $L/results
STEPS=("$@"); [[ ${#STEPS[@]} -eq 0 ]] && STEPS=(complete lookup score panels detb harvest aucref remeasure oracle gallery rank reliability)
die() { echo "POST FAIL [$land $1] $(date -Is)"; exit 1; }
run() { local name=$1; shift; echo "--- $land $name $(date -Is)"
        $L/dk.sh post "mland_post_${land}_$name" "$@" > "$LOGS/$name.log" 2>&1 || die "$name (log $LOGS/$name.log)"; }
pv() { local f; f=$(ls $1/GEA/tables/methods/$2/$2_pvalues_K*.tsv 2>/dev/null | grep -v '_sig_\|_wza'); [[ $(wc -w <<< "$f") -eq 1 ]] && echo "$f"; }

for step in "${STEPS[@]}"; do case $step in
complete)
    bad=""
    while read -r s; do for st in up gea off; do
        [[ "$(cat $L/$st/MVP${s}hd/STATUS 2>/dev/null)" == DONE* ]] || bad+="$st:$s "; done
        [[ -d $EVH/$OFF/gardens/$s ]] || bad+="gardens:$s "
    done < $SEEDF
    echo "complete: $(wc -l < $SEEDF) seeds ${bad:+missing: $bad}"; [[ -z $bad && $(wc -l < $SEEDF) -eq 600 ]] || die complete ;;
lookup)   # built from ld_window_by_replicate.tsv; must equal the window every gea set used
    n=0; bad=0
    while IFS=$'\t' read -r s d; do [[ $s == seed ]] && continue; n=$((n+1))
        for set in gea.solo_lfmm gea.solo_emmax gea.solo_rda gea.union gea.best gea.intersect3; do
            y=$(awk '/^clumping_distance:/{print $2; exit}' MVP${s}hd_results/_intermediate/snp_sets/$set/snp_set.yaml)
            [[ "$y" == "$d" ]] || bad=$((bad+1)); done
    done < $LDBP
    echo "lookup: $n seeds x 6 sets, $bad mismatch(es) against snp_set.yaml"; [[ $n -eq 600 && $bad -eq 0 ]] || die lookup ;;
score)
    run score --memory=200g --cpus=96 -e PIPELINE_ROOT=/pipeline -- Rscript benchmarks/eval_offset_lind.R \
        --seeds=$SEEDS --outdir=$EV/$OFF --ncores=96 ;;
panels)
    run panels --memory=32g --cpus=8 -e PIPELINE_ROOT=/pipeline -e MVP_PANELS=hd -e MVP_LINKED_BP=/pipeline/$LDBP -- \
        Rscript $SC/panel_tables_land.R --seeds=$SEEDS --outdir=$EV/$OFF ;;   # copy: no control filter (eval_offset_lind.R:191)
detb)
    run detb --memory=64g --cpus=32 -e MVP_LINKED_BP=/pipeline/$LDBP -- \
        Rscript $SC/detection_b_land.R /pipeline/$SEEDF hd /pipeline/$L/results/detection_b_$land 32 ;;
harvest)  # the c1 layout the remeasure reads: rdaunc = hd's three tables; onefit = hd LFMM/EMMAX + hdrc RDA
    n=0
    while read -r s; do
        H=MVP${s}hd_results; R=MVP${s}hdrc_results
        mkdir -p $EVH/$PU/MVP$s/c1 $EVH/$PO/MVP$s/c1
        for m in LFMM EMMAX RDA; do f=$(pv $H $m) || die "harvest $s $m"
            ln -f "$f" $EVH/$PU/MVP$s/c1/${m}_pvalues.tsv; done
        for m in LFMM EMMAX; do ln -f "$(pv $H $m)" $EVH/$PO/MVP$s/c1/${m}_pvalues.tsv; done
        f=$(pv $R RDA) || die "harvest $s hdrc RDA"; ln -f "$f" $EVH/$PO/MVP$s/c1/RDA_pvalues.tsv
        cmp -s $EVH/$PU/MVP$s/c1/RDA_pvalues.tsv $EVH/$PO/MVP$s/c1/RDA_pvalues.tsv && die "harvest $s: hd and hdrc RDA identical"
        n=$((n+1))
    done < $SEEDF
    echo "harvest: $n seeds -> $PU, $PO (hard links)"; [[ $n -eq 600 ]] || die harvest ;;
aucref)
    run aucref --memory=48g --cpus=16 "${ARMENV[@]}" -e PARAMS_DIR=$EV/$PU -e OUT_DIR=$EV/detection600_$land -e NCORES=16 -- \
        Rscript $SC/detection_600_land.R ;;
remeasure)
    for arm in rdaunc onefit; do P=$PU; [[ $arm == onefit ]] && P=$PO
        run remeasure_$arm --memory=96g --cpus=16 "${ARMENV[@]}" -e PARAMS_DIR=$EV/$P -e ARM=${arm}_rda2x_ldb \
            -e OUT_DIR=$EV/$RM/${arm}_rda2x_ldb -e A3_REF=$EV/detection600_$land/aucpr_per_seed.tsv -e A3_METHODS=LFMM,EMMAX \
            -e TOP_SHARES=0.001,0.0025,0.005,0.01,0.02 -e SHARE_MULT=RDA=2 -e MVP_LINKED_BP=/pipeline/$LDBP -e MVP_SEED_WINDOW=1 \
            -e CELL=c1 -e NCORES=16 -- Rscript $SC/remeasure_600_land.R
    done ;;
oracle)
    run oracle --memory=48g --cpus=8 "${ARMENV[@]}" -e OPENBLAS_NUM_THREADS=4 -e OFFSET_DIR=$OFF \
        -e FIG_OUT=$EV/$FIG -e STATS_DIR=$EV/$FIG -- Rscript benchmarks/mvp_oracle_stats.R ;;
gallery)
    run gallery --memory=48g --cpus=8 "${ARMENV[@]}" -e OPENBLAS_NUM_THREADS=4 -e OFFSET_DIR=$OFF \
        -e J15_DIR=$EV/$FIG -e AUC_FILE=$EV/$RM/rdaunc_rda2x_ldb/rank_metrics.tsv \
        -e FIG_OUT=$EV/$GAL -- Rscript $SC/j16_gallery_land.R ;;
rank)
    run rank --memory=16g --cpus=4 "${ARMENV[@]}" -e OFFSET_DIR=$OFF -e FIG_OUT=$EV/$FIG/rank -- Rscript $SC/panel_rank_land.R ;;
reliability)
    run reliability --memory=64g --cpus=8 "${ARMENV[@]}" -e OFFSET_DIR=$OFF -e FIG_OUT=$EV/$FIG/rank -- \
        Rscript benchmarks/mvp_rank_reliability.R ;;
*) die "unknown step $step" ;;
esac; done
echo "=== post $land done $(date -Is)"
