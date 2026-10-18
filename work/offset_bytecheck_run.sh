#!/usr/bin/env bash
# Run the offset scripts into one scratch mirror (2026-10-18 byte gate). Inside the image:
#   bash /pipeline/work/offset_bytecheck_run.sh <default|minou> <pre|post>
# default = no flags (as fig3/logs/bytecheck_v1.log): build_site_maps.R + build_paper_pack.R.
# minou   = the env and the three scripts of fig3/run_offset_minou.sh; post adds MAIN6_MS=1.
set -euo pipefail
cd /pipeline
env=$1; tag=$2
export AN=work/offset_bytecheck_2026-10-18/${env}_${tag}
[ -d "/pipeline/$AN" ] || { echo "no mirror $AN"; exit 1; }
V1=ANALYSIS/Trifolium86MainK5/maladaptationER
if [ "$env" = minou ]; then
    export PROJ=Trifolium86MainK5ERpc5 PALETTE=minou
    export SITE_NAMES=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.tsv KCHOICE_PANELS=a DELTA_CORR=1 NO_LINEAR=1 \
        MAIN_MS=1 SUPP_MS=1
    [ "$tag" = post ] && export MAIN6_MS=1
fi
Rscript $V1/build_site_maps.R
Rscript $V1/build_paper_pack.R
[ "$env" = minou ] && Rscript ANALYSIS/Trifolium86MainK5ERpc5/maladaptation/fig_k_choice.R
echo done
