#!/usr/bin/env bash
# Structure panels for Fig 3 v2 (plan ~/.claude/plans/fancy-scribbling-garden.md, Phase 1), drawn by the derived copies
# of the pipeline's own plot scripts with the same inputs as rules pca_structure_plot / structure_barplot
# (workflow/rules/prestructure.smk:62-90) and the piemap rule, except that `site` is the short name
# (fig3/inputs/*.short.tsv, written by fig3/site_names.py).
# v3 (idempotent-conjuring-noodle.md Phase 2): PCA without site labels (user); piemap labels 7 pt halo text.
# v4 (2026-10-18, plan let-s-implement-this-plan-mighty-bubble.md): Fig 3a PCA at 80 mm wide (legend below, theme_ms() default) (D13);
# piemap raster legend title = display name, not the layer name (D14).
# v5 (2026-10-18, user review): Fig 3a PCA without its legend, the piemap's Clusters legend keys both (D20).
# Run INSIDE the image: bash /pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/run_structure.sh
set -euo pipefail
R=/pipeline/Trifolium86MainK5ERpc5_results
F=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3
D=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/manuscript_fig3/draft
PCA=$R/_work/maf0.1_miss0.3_smiss0.8_rel0.999/ld0.2_win100_step10/Trifolium86_chr1-6.pca/Trifolium86_chr1-6
CL=$F/inputs/clusters_K5.short.tsv

Rscript $F/plot_pca_structure_ms.R $CL $PCA.projections $PCA.eigenvalues 5 $D/st_pcapie12 80x112 1 2 F none
Rscript $F/plot_pca_structure_ms.R $CL $PCA.projections $PCA.eigenvalues 5 $D/st_pcapie34 88x80 3 4 F
Rscript $F/plot_structure_ms.R $CL 5 $D/st_bars_pipeline 180x45 pipeline
Rscript $F/plot_structure_ms.R $CL 5 $D/st_bars_split 180x45 split
Rscript $F/plot_structure_ms.R $CL 5 $D/st_bars_split_site 180x45 split_site

O=$D/piemap/
declare -A LEG=([bio_3]=$'Isothermality\nBIO3 (%)')   # display title (D14); the layer name stays arg 2
for b in bio_3; do   # bio_2 dropped (user review v2, 2026-10-17)
    Rscript $F/plot_piemap_ms.R \
        $R/climate/rasters/present/climate_present_rasterstack.tif "$b" "${LEG[$b]}" \
        $F/inputs/metadata_climate.short.tsv $CL NULL NULL \
        0.6 T 7 "$O" "${O}_qs/" "piemap_K5_$b" NULL 1.0 Clusters
done
