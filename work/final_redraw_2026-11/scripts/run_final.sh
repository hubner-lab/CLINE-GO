#!/usr/bin/env bash
# run_final.sh <repro|new> <step ...> -- the final redraw's drawing steps (lab, 2026-11-11), one container each through
# frdk.sh (named fr_<variant>_<step>). repro = the copy that drew the published file, unchanged (only its output folder
# redirected where it was hardcoded); new = this job's copy. Every run goes through dump_layers.R, so the plotted
# values of both are on disk (scratch/<variant>/<step>/dump/) for compare_dumps.py.
# Outputs: scratch/<variant>/<step>/ ; logs: logs/<variant>_<step>.log.
#   j1 s_fp_pca (paper) | j2 s_qc (latest) | j34 s_structure_details + s_random_sets_climate (paper)
#   j5 Fig3 panel c (paper) | j68 Fig5 + s_predictor_change (latest) | j7 s_rank_architecture parts (paper)
#   j9 s_emmax_by_predictor (latest)
set -uo pipefail
cd /mnt/data/eugene/ADAPTOGENE
V=${1:?variant: repro | new}; shift
[[ $V == repro || $V == new ]] || { echo "variant must be repro or new"; exit 2; }
W=work/final_redraw_2026-11
P=/pipeline/$W
S=$P/scripts
DK=$W/scripts/frdk.sh
DUMP=(Rscript $S/dump_layers.R)
die() { echo "FINAL FAIL [$V $1] $(date -Is)"; exit 1; }
run() { local img=$1 name=$2; shift 2
        local out=$W/scratch/$V/$name; mkdir -p "$out"
        echo "--- $V $name $(date -Is)"
        "$DK" "$img" "${V}_$name" "$@" > "$W/logs/${V}_$name.log" 2>&1 || die "$name (log $W/logs/${V}_$name.log)"
        echo "    ok $(date -Is); SOFT ASSERT lines: $(grep -c 'SOFT ASSERT' "$W/logs/${V}_$name.log"); Warning lines: $(grep -c '^Warning' "$W/logs/${V}_$name.log")"; }
D150=ANALYSIS/Trifolium86MainK5d150/fig3/offset_minou
MINOU=(-e AN=$D150 -e PROJ=Trifolium86MainK5d150 -e PALETTE=minou
       -e SITE_NAMES=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.tsv -e KCHOICE_PANELS=a
       -e DELTA_CORR=1 -e NO_LINEAR=1 -e MAIN_MS=1 -e SUPP_MS=1 -e MAIN6_MS=1 -e SUPP2_MS=1)

for step in "$@"; do O=$P/scratch/$V/$step; case $step in
j1)  # s_fp_pca: benchmarks/mvp_land_tracking_fig.R in FP_PCA_OUT mode (Data S1 Job B); repro = git HEAD's copy
     if [[ $V == repro ]]; then SCR=$S/orig/mvp_land_tracking_fig.R; else SCR=/pipeline/benchmarks/mvp_land_tracking_fig.R; fi
     run paper j1 --memory=16g --cpus=2 -e FP_PCA_OUT=$O -e DUMP_SCRIPT=$SCR -e DUMP_OUT=$O/dump -e DUMP_OBJECTS=p6 \
         -- "${DUMP[@]}" ;;
j9)  # s_emmax_by_predictor: build_fig3.R only=as_pred14_v2 (closeout step s13); repro = the closeout's patched copy
     if [[ $V == repro ]]; then SCR=/pipeline/work/closeout_lab/results/scripts/patched/build_fig3.R; else SCR=$S/build_fig3_s11.R; fi
     run latest j9 --memory=32g --cpus=4 -e CLOSEOUT_OUT=$O -e DUMP_SCRIPT=$SCR -e DUMP_OUT=$O/dump -e DUMP_OBJECTS=p \
         -e DUMP_ARGS=only=as_pred14_v2 -- "${DUMP[@]}" ;;
j68) # Fig5 (fig5_offset_ms) + s_predictor_change + s_random_sets: build_paper_pack_d150.R with the MINOU environment
     # (closeout step s17); repro = the closeout's patched copy
     if [[ $V == repro ]]; then SCR=/pipeline/work/closeout_lab/results/scripts/patched/build_paper_pack_d150.R; else SCR=$S/build_paper_pack_final.R; fi
     run latest j68 --memory=60g --cpus=8 "${MINOU[@]}" -e SOFT_LOG=$O/soft_asserts.tsv -e CLOSEOUT_OUT=$O \
         -e DUMP_SCRIPT=$SCR -e DUMP_OUT=$O/dump -e DUMP_OBJECTS=fig5,p6,p7 -- "${DUMP[@]}" ;;
j34) # s_structure_details + s_random_sets_climate: benchmarks/mvp_land_data_s1_trifolium.R (Data S1 F3/F4), MINOU env;
     # repro = git HEAD's copy. Members (CHECK_OUT) are compared byte for byte with the published ones.
     if [[ $V == repro ]]; then SCR=$S/orig/mvp_land_data_s1_trifolium.R; else SCR=/pipeline/benchmarks/mvp_land_data_s1_trifolium.R; fi
     run paper j34 --memory=60g --cpus=8 "${MINOU[@]}" -e FIG_OUT=$O -e CHECK_OUT=$O/members -- Rscript $SCR ;;
j7)  # s_rank_architecture_parts_{a,b,c}: benchmarks/mvp_land_figures.R FIG_MODE=archtitle (new); the untitled pages
     # (CHECK_OUT) are FIG_MODE=supp's call and are compared byte for byte with the published parts. new only.
     run paper j7 --memory=32g --cpus=4 -e FIG_MODE=archtitle -e FIG_OUT=$O -e CHECK_OUT=$O/members \
         -- Rscript /pipeline/benchmarks/mvp_land_figures.R ;;
j5)  # Fig3 panel c (as_miami_v2): the 150 kb build_fig3.R only=as_miami_v2 (Phase 7 p7_gea_figs.sh), OUT redirected;
     # composed into Fig3_structure_association by compose_fig3.sh
     if [[ $V == repro ]]; then SCR=$S/build_fig3_d150_repro.R; else SCR=$S/build_fig3_d150.R; fi
     run latest j5 --memory=40g --cpus=8 -e FR_OUT=$O -e SOFT_LOG=$O/soft_asserts.tsv -e DUMP_SCRIPT=$SCR -e DUMP_OUT=$O/dump \
         -e DUMP_OBJECTS=p -e DUMP_ARGS=only=as_miami_v2 -- "${DUMP[@]}" ;;
j2)  # s_qc: the audit 3b copy of fig3/plot_qc_ms.R with run_qc_ms.sh's 25 pipeline arguments + out stem + counts;
     # new = plot_qc_ms_final.R, + argument 28 (inputs/snp_missing_57.lmiss.gz from j2_snp_missing.sh)
     R=/pipeline/Trifolium86MainK5ERpc5_results; WK=$R/_work/maf0.1_miss0.3_smiss0.8_rel0.999
     A=($R/Processing/tables/vcf_raw_summary.tsv $R/_intermediate/samples/samples_missing_stats.tsv
        $R/_intermediate/qc/snp_missingness_raw.lmiss $R/_intermediate/qc/maf_raw.frq $R/_intermediate/qc/maf_filtered.frq
        $R/_intermediate/qc/het_raw.het $R/_intermediate/qc/het_filtered.het $R/_intermediate/qc/snp_density_raw.snpden
        $R/_intermediate/qc/snp_density_filtered.snpden $R/Processing/tables/metadata.tsv $WK/Trifolium86_chr1-6.vcf
        $WK/ld0.2_win100_step10/Trifolium86_chr1-6.vcf 0.8 0.3 0.1 NULL TRUE $R/_intermediate/qc/depth_per_sample.idepth
        $R/_intermediate/qc/depth_per_site.ldepth.mean $R/Processing/plots/ $R/Processing/tables/
        $R/Processing/tables/relatedness_pairs.tsv 0.999 $R/Processing/tables/relatedness_removed.tsv remove
        $O/s_qc $O/s_qc_counts.tsv)
     if [[ $V == repro ]]; then SCR=/pipeline/work/ld_window_lab/scratch_audit3b/qc/plot_qc_ms.R
     else SCR=$S/plot_qc_ms_final.R; A+=($P/inputs/snp_missing_57.lmiss.gz); fi
     ARGS=$(IFS=$'\t'; echo "${A[*]}")
     run latest j2 --memory=48g --cpus=4 -e DUMP_SCRIPT=$SCR -e DUMP_OUT=$O/dump -e "DUMP_ARGS=$ARGS" \
         -e DUMP_OBJECTS=$([[ $V == repro ]] && echo p_smiss,p_rel,p_mds,p_density,p_maf || echo p_smiss,p_rel,p_mds,p_density,p_snpmiss,p_maf) \
         -- "${DUMP[@]}" ;;
j7merge) # check only: FIG_MODE=merge after page_title() was hoisted; its pages must equal the Data S1 published ones
     run paper j7merge --memory=32g --cpus=4 -e FIG_MODE=merge -e FIG_OUT=$O -e CHECK_OUT=$O/members \
         -- Rscript /pipeline/benchmarks/mvp_land_figures.R ;;
*) die "unknown step $step" ;;
esac; done
echo "=== $V done $(date -Is)"
