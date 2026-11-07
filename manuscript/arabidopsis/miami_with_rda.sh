#!/usr/bin/env bash
# miami_with_rda.sh: the GEA x GWAS Miami plot with RDA, drawn outside the pipeline for the manuscript.
#
# workflow/rules/gea_x_gwas.smk:25 passes the climate predictors (bio_1 to bio_19) as the GEA traits, so
# plot_miami.R finds nothing under RDA's one column, climate_multivariate, and leaves RDA out (filed as a
# [bug] in docs/pipeline_improvement_requests.md). This runs the image's own plot_miami.R with the rule's
# arguments as `snakemake -n -p --forcerun miami_plot` prints them on the run's image, except two: the
# GEA traits end in climate_multivariate, the list combine_predictors("GEA") gives (common.smk:2318), and
# the plot goes to ANALYSIS/ instead of the results tree. No pipeline code changes; the results tree and
# the app keep the plot without RDA.
#
# Run from the repo root, with no mode running:
#   manuscript/arabidopsis/miami_with_rda.sh
# Returns once the container is started. It writes, under ANALYSIS/Arabidopsis1135/:
#   figures/miami/miami_combined_K6{.png,.svg,_background.png,_coords.json}
#   runlogs/miami_with_rda.log, and runlogs/miami_with_rda.exit when the plot is done
set -euo pipefail

IMAGE=cline-go:6cdedcc
NAME=arab-miami-rda
RESULTS=/data/Arabidopsis1135_results
OUT=ANALYSIS/Arabidopsis1135/figures/miami
RUNLOGS=ANALYSIS/Arabidopsis1135/runlogs

GEA="EMMAX:bonf_0.05:$RESULTS/GEA/tables/methods/EMMAX/EMMAX_pvalues_K6.tsv"
GEA+=",LFMM:bonf_0.05:$RESULTS/GEA/tables/methods/LFMM/LFMM_pvalues_K6.tsv"
GEA+=",RDA:bonf_0.05:$RESULTS/GEA/tables/methods/RDA/RDA_pvalues_K6.tsv"
GWAS="EMMAX:bonf_0.05:$RESULTS/GWAS/tables/methods/EMMAX/EMMAX_pvalues_K6.tsv"
GWAS+=",BLINK:bonf_0.05:$RESULTS/GWAS/tables/methods/BLINK/BLINK_pvalues_K6.tsv"
GEA_TRAITS="bio_1,bio_2,bio_3,bio_4,bio_5,bio_6,bio_7,bio_8,bio_9,bio_10,bio_11,bio_12,bio_13,bio_14"
GEA_TRAITS+=",bio_15,bio_16,bio_17,bio_18,bio_19,climate_multivariate"
GWAS_TRAITS="FT10,FT16"

docker_cli() { nix shell nixpkgs#docker-client -c docker "$@"; }
if docker_cli ps --format '{{.Names}}' | grep -qx "$NAME"; then
    echo "miami_with_rda.sh: container $NAME is already running" >&2
    exit 1
fi

mkdir -p "$OUT" "$RUNLOGS"
rm -f "$RUNLOGS/miami_with_rda.exit"
# The plot rule runs on CPU threads; scattermore's blur reads OMP_NUM_THREADS (scattermore_threads.R),
# and the render is the same at any thread count.
docker_cli run -d --rm --name "$NAME" --user "$(id -u):$(id -g)" -e USER=pipeline -e HOME=/tmp \
    -e OMP_NUM_THREADS=16 -e OPENBLAS_NUM_THREADS=16 --memory 120g \
    -v "$PWD":/data -w /data "$IMAGE" \
    bash -c 'Rscript /pipeline/scripts/plot_miami.R "$1" "$2" "$3" "$4" 6 "/data/$5/" > "/data/$6/miami_with_rda.log" 2>&1
             echo $? > "/data/$6/miami_with_rda.exit"' \
    _ "$GEA" "$GWAS" "$GEA_TRAITS" "$GWAS_TRAITS" "$OUT" "$RUNLOGS"
