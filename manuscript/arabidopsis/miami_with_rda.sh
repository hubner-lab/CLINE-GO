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
# Run from the repo root, with no mode of <project> running:
#   manuscript/arabidopsis/miami_with_rda.sh <project> [rda-second]
# rda-second puts climate_multivariate (RDA) second in the GEA trait order instead of last, so plot_miami.R's
# trait palette (turbo by position above 8 traits) gives RDA an indigo instead of a red next to the GWAS
# traits', and RDA's cloud is drawn second, under most per-predictor clouds. Not first: plot_miami.R:68-70
# takes its coordinate system from the first trait of the first method (EMMAX), which has no
# climate_multivariate column, so RDA first drew an empty plot with exit 0 (filed in
# docs/pipeline_improvement_requests.md, 2026-11-15). Outputs: figures/miami_rda_second/ and
# runlogs/miami_with_rda_second.{log,exit}.
# K is the project's sNMF.k_best (the K in the p-value tables' names). The methods, predictors and
# traits below are the Arabidopsis configs' (GEA EMMAX/LFMM/RDA on 19 BIO, GWAS EMMAX/BLINK on FT10
# and FT16). Returns once the container is started. It writes, under ANALYSIS/<project>/:
#   figures/miami/miami_combined_K<k>{.png,.svg,_background.png,_coords.json}
#   runlogs/miami_with_rda.log, and runlogs/miami_with_rda.exit when the plot is done (_second for rda-second)
set -euo pipefail

PROJECT=${1:?usage: miami_with_rda.sh <project> [rda-second]}
VARIANT=${2:-}
[ -z "$VARIANT" ] || [ "$VARIANT" = rda-second ] || { echo "miami_with_rda.sh: unknown variant $VARIANT" >&2; exit 2; }
CONFIG=config_$PROJECT.yaml
[ -f "$CONFIG" ] || { echo "miami_with_rda.sh: no $CONFIG in $PWD" >&2; exit 2; }
K=$(awk '/^sNMF:/ {s = 1; next} s && /^[^ ]/ {s = 0} s && $1 == "k_best:" {print $2}' "$CONFIG")
[[ "$K" =~ ^[0-9]+$ ]] || { echo "miami_with_rda.sh: no integer sNMF.k_best in $CONFIG" >&2; exit 2; }

IMAGE=cline-go:6cdedcc
NAME=arab-miami-rda-$PROJECT${VARIANT:+-$VARIANT}
RESULTS=/data/${PROJECT}_results
OUT=ANALYSIS/$PROJECT/figures/miami${VARIANT:+_rda_second}
RUNLOGS=ANALYSIS/$PROJECT/runlogs

GEA="EMMAX:bonf_0.05:$RESULTS/GEA/tables/methods/EMMAX/EMMAX_pvalues_K$K.tsv"
GEA+=",LFMM:bonf_0.05:$RESULTS/GEA/tables/methods/LFMM/LFMM_pvalues_K$K.tsv"
GEA+=",RDA:bonf_0.05:$RESULTS/GEA/tables/methods/RDA/RDA_pvalues_K$K.tsv"
GWAS="EMMAX:bonf_0.05:$RESULTS/GWAS/tables/methods/EMMAX/EMMAX_pvalues_K$K.tsv"
GWAS+=",BLINK:bonf_0.05:$RESULTS/GWAS/tables/methods/BLINK/BLINK_pvalues_K$K.tsv"
for spec in ${GEA//,/ } ${GWAS//,/ }; do
    table=${spec##*:}
    [ -f "${table#/data/}" ] || { echo "miami_with_rda.sh: no ${table#/data/}" >&2; exit 2; }
done
GEA_TRAITS="bio_1,bio_2,bio_3,bio_4,bio_5,bio_6,bio_7,bio_8,bio_9,bio_10,bio_11,bio_12,bio_13,bio_14"
GEA_TRAITS+=",bio_15,bio_16,bio_17,bio_18,bio_19,climate_multivariate"
[ -z "$VARIANT" ] || GEA_TRAITS="bio_1,climate_multivariate,${GEA_TRAITS#bio_1,}"
[ -z "$VARIANT" ] || GEA_TRAITS=${GEA_TRAITS%,climate_multivariate}
LOG=miami_with_rda${VARIANT:+_second}
GWAS_TRAITS="FT10,FT16"

docker_cli() { nix shell nixpkgs#docker-client -c docker "$@"; }
if docker_cli ps --format '{{.Names}}' | grep -qx "$NAME"; then
    echo "miami_with_rda.sh: container $NAME is already running" >&2
    exit 1
fi

mkdir -p "$OUT" "$RUNLOGS"
rm -f "$RUNLOGS/$LOG.exit"
# The plot rule runs on CPU threads; scattermore's blur reads OMP_NUM_THREADS (scattermore_threads.R),
# and the render is the same at any thread count.
docker_cli run -d --rm --name "$NAME" --user "$(id -u):$(id -g)" -e USER=pipeline -e HOME=/tmp \
    -e OMP_NUM_THREADS=16 -e OPENBLAS_NUM_THREADS=16 --memory 120g \
    -v "$PWD":/data -w /data "$IMAGE" \
    bash -c 'Rscript /pipeline/scripts/plot_miami.R "$1" "$2" "$3" "$4" "$7" "/data/$5/" > "/data/$6/$8.log" 2>&1
             echo $? > "/data/$6/$8.exit"' \
    _ "$GEA" "$GWAS" "$GEA_TRAITS" "$GWAS_TRAITS" "$OUT" "$RUNLOGS" "$K" "$LOG"
