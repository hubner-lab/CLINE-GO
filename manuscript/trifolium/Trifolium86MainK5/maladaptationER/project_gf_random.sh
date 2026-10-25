#!/usr/bin/env bash
# Project each arm's size-matched RANDOM Gradient Forest model onto the future climate,
# with the pipeline's own scripts/gradient_forest_offset.R (the rule only projects the
# adaptive model; the random one feeds the cumimp plot alone). Output per arm/set:
#   <arm>/gf_random/<set>/{offset_raster.tif, genetic_offset_map.tsv, genetic_offset_site.tsv, diagnostics.tsv}
# A random model written as a sentinel (no SNP with R2 > 0) makes the script stop; that
# is recorded in <arm>/gf_random/<set>/STATUS and the loop goes on.
# Run inside the image from /pipeline:
#   bash ANALYSIS/Trifolium86MainK5/maladaptationER/project_gf_random.sh [arm ...]
set -u
cd /pipeline
B=${AN:-ANALYSIS/Trifolium86MainK5/maladaptationER}   # env AN/PROJ: v2 = PC5 rebuild (2026-10-16)
R=${PROJ:-Trifolium86MainK5ER}_results
SC=year2061-2080_ssp585
FUT=$R/climate/tables/future/climate_future_${SC}_all.tsv
PRES=$R/climate/tables/present/climate_present_all.tsv
RAST=$R/climate/rasters/present/climate_present_rasterstack.tif
SAMP=$R/_intermediate/samples/metadata_climate_valid.tsv
declare -A PRED=([arm16]="bio_1,bio_2,bio_3,bio_4,bio_5,bio_6,bio_7,bio_8,bio_9,bio_10,bio_11,bio_12,bio_13,bio_15,bio_16,bio_19"
                 [arm3]="bio_1,bio_3,bio_5" [arm12]="bio_12,bio_3,bio_5")
# arms as argv (default: arm16 arm3, the first run)
ARMS=${*:-arm16 arm3}
for arm in $ARMS; do
  for set in gea.ER_overlap gea.ER_rda1e4 gea.ER_rda1e3; do
    model=$B/$arm/intermediate/gradient_forest/${set}_nospatial/random_model.qs
    out=$B/$arm/gf_random/$set
    mkdir -p "$out"
    if [ ! -s "$model" ]; then echo "missing model $model" > "$out/STATUS"; continue; fi
    if Rscript scripts/gradient_forest_offset.R "$model" "${PRED[$arm]}" "$FUT" "$PRES" "$RAST" "$SAMP" \
        "$out/offset_raster.tif" "$out/genetic_offset_map.tsv" "$out/genetic_offset_site.tsv" \
        "$SC" TRUE "$out/diagnostics.tsv" > "$out/project.log" 2>&1; then
      echo ok > "$out/STATUS"
    else
      echo "failed: $(tail -2 "$out/project.log" | tr '\n' ' ')" > "$out/STATUS"
    fi
    echo "$arm $set: $(cat "$out/STATUS")"
  done
done
