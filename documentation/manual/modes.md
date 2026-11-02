# The modes

The pipeline runs in ten modes, one `clinego run <mode>` each. Every mode writes into
`<project>_results/` and logs under `<project>_logs/`, one folder per module. The commands below use the config of
[Run on your own data](../own-data.md), `config_myproject.yaml`.

| Mode | Needs `sNMF.k_best` | Needs climate (`Regime.mode: standard`) |
|---|---|---|
| [processing](#processing) | no | no |
| [prestructure](#prestructure) | no | no |
| [structure](#structure) | yes | no; climate extraction only in `standard` |
| [climate](#climate) | yes | yes |
| [traits](#traits) | no | no; one plot only in `standard` |
| [pregea](#pregea) (optional) | yes | yes |
| [gea](#gea) | yes | yes |
| [gwas](#gwas) | yes | no; trait maps only in `standard` |
| [gea_x_gwas](#gea_x_gwas) | yes | no; the Miami plot only in `standard` |
| [maladaptation](#maladaptation) | yes | yes |

Run them in the order of the table, one at a time for a project: modes share the project's
results folder and its lock. A mode that needs `sNMF.k_best` stops at once without it; a
climate mode stops at once in a `gwas_only` project. Settings are described in the
[configuration reference](configuration.md).

## processing

Filters the VCF (sample and SNP missingness, minor allele frequency, optional depth and
heterozygosity filters, relatedness), prunes it for linkage disequilibrium, and writes QC plots
of what each step removed.

- **Needs:** the VCF and the metadata named in `Input.*`; the GFF too, when given, which it
  normalizes to the VCF's chromosome names.
- **Main outputs:**
  - `Processing/tables/filtering_summary.tsv`: samples and SNPs after each filter
  - `Processing/tables/metadata.tsv`: the samples kept
  - `Processing/plots/`: missingness, heterozygosity, allele frequency, SNP density,
    relatedness and the filtering attrition

```{code-block} text
|run_own_processing|
```

## prestructure

Runs PCA and sNMF on the LD-pruned SNPs, for every *K* from `sNMF.k_start` to `sNMF.k_end`.
This is where you choose `sNMF.k_best`.

- **Needs:** `processing`.
- **Main outputs:**
  - `PreStructure/tables/cross_entropy_K<start>-<end>.tsv` and its plot: cross-entropy for
    every *K* and run, to choose *K*
  - `PreStructure/plots/pca.png`, and `tracy_widom.png` with
    `PreStructure/tables/tracy_widom.tsv`: the PCA and the significance of each axis
  - `PreStructure/plots/K<k>/` and `PreStructure/tables/K<k>/clusters_K<k>.tsv`: ancestry bar
    plots, the PCA coloured by ancestry, and the ancestry coefficients, for each *K*

```{code-block} text
|run_own_prestructure|
```

## structure

Imputes missing genotypes from the ancestry at `sNMF.k_best`, measures how linkage
disequilibrium decays with distance, and in a `standard` project extracts the climate at each
sample and maps ancestry onto the sampling sites.

- **Needs:** `prestructure` and `sNMF.k_best`. In `standard`, the WorldClim rasters, which it
  downloads into `Input.dir` the first time (or `Climate.source: custom` tables).
- **Main outputs:**
  - `Structure/tables/ld_decay_half_distances.tsv` and `Structure/plots/ld_decay/`: the LD
    decay distances, which set the distance used to merge significant SNPs in `gea` and `gwas`
  - `climate/tables/present/climate_present_site.tsv`: the climate at each site (`standard`)
  - `Structure/plots/piemap/`: ancestry pie charts over each bioclimatic layer (`standard`)
  - with `Population.calc_stats: true`, Tajima's D, nucleotide diversity, AMOVA and (in
    `standard`) a Mantel test, under `Structure/`

```{code-block} text
|run_own_structure|
```

## climate

Characterizes the climate predictors: their distributions, correlations and any that do not
vary across sites, and splits the genetic variation into climate, structure and geography
(variance partitioning with dbMEM spatial variables).

- **Needs:** `structure`, `Regime.mode: standard`, `Climate.predictors`. It stops if a
  predictor does not vary across the samples.
- **Main outputs:**
  - `climate/plots/correlation_heatmap.png` and `climate/tables/present/climate_correlation.tsv`
  - `climate/tables/present/climate_invariant_predictors.tsv`
  - `climate/tables/varpart/variance_partition.tsv` and `climate/plots/varpart/`
  - `climate/tables/spatial/dbmem_vectors.tsv`: the spatial variables a spatial Gradient
    Forest uses
- Run it before `maladaptation` when the Gradient Forest is spatial
  (`spatial_correction: with` or `both`): nothing checks this order for you.

```{code-block} text
|run_own_climate|
```

## traits

Characterizes the measured traits: distributions, correlations, a pairs grid, and any trait
with a single value. It is the phenotype counterpart of `climate` and runs in both regimes.

- **Needs:** trait columns in the metadata (column 5 on); not `sNMF.k_best`.
- **Main outputs:**
  - `Traits/tables/trait_summary.tsv` and `Traits/tables/trait_invariant.tsv`
  - `Traits/plots/`: trait densities, the trait correlogram and the pairs grid (a placeholder
    above `Traits.pairs_max_factors` traits)
  - in `standard`, the traits-by-climate correlogram, which needs `structure` to have run

```{code-block} text
|run_own_traits|
```

## pregea

Optional. Explores the settings the association methods depend on, on the LD-pruned SNPs: the
number of latent factors for LFMM (`sNMF.k_best` plus or minus `PreGEA.k_offset`), and the
number of principal components for EMMAX and for RDA's structure correction (0 to
`PreGEA.n_pcs_max`). It recommends one value per method.

- **Needs:** `structure`, `sNMF.k_best`, `Regime.mode: standard`.
- **Main outputs:**
  - `PreGEA/tables/pregea_recommendations.tsv`: the recommended value per method, which the
    app's GEA settings offer to apply
  - `PreGEA/plots/` and `PreGEA/tables/`: the ladders behind each recommendation
- The mode has not yet been verified end to end; read its recommendations as a guide, not a
  rule.

```{code-block} text
|run_own_pregea|
```

## gea

Genotype-environment association: tests every SNP against every climate predictor with each
method in `GEA.configs`, then merges the significant SNPs into regions and annotates them with
genes.

- **Needs:** `structure`, `sNMF.k_best`, `Regime.mode: standard`, `Climate.predictors` and
  `GEA.configs`. With `snp_clumping_distance: auto_genome_wide` (the default), the LD decay
  table from `structure`.
- **Main outputs:**
  - `GEA/tables/methods/<method>/`: p-values, significant SNPs, the method's resolved settings,
    genomic inflation and thresholds
  - `GEA/tables/selected_snps.tsv`, `regions_per_trait.tsv`, `regions_combined.tsv`
  - `GEA/tables/genes_combined.tsv` and `genes_per_region.tsv` (with a GFF)
  - `GEA/plots/manhattan/`: Manhattan and QQ plots per method, and combined
  - `_intermediate/snp_sets/gea.<name>/`: each SNP set declared under `GEA.snp_sets`, for
    `maladaptation`

```{code-block} text
|run_own_gea|
```

## gwas

Genome-wide association with the measured traits, with each method in `GWAS.configs`, then the
same regions and genes as `gea`.

- **Needs:** `structure`, `sNMF.k_best`, trait columns in the metadata and `GWAS.configs`.
  Runs in both regimes. `GWAS.missing_strategy` decides what happens to samples missing a
  trait (`DROP` leaves them out of that trait's test).
- **Main outputs:**
  - `GWAS/tables/phenotype_missing_summary.tsv`
  - `GWAS/tables/methods/<method>/`, `GWAS/tables/selected_snps.tsv`, the regions and genes
    tables, and `GWAS/plots/manhattan/`, as for `gea`
  - in `standard`, `GWAS/plots/piemap/`: trait values mapped onto the sites

```{code-block} text
|run_own_gwas|
```

## gea_x_gwas

Compares the association results: which pairs of traits and predictors share significant SNPs,
and a Miami plot of environmental against phenotypic signal.

- **Needs:** both `GEA.configs` and `GWAS.configs` in the config, and `gea` and `gwas` run
  first. The mode's rules are loaded only when both lists are present, so keep the template's
  `GEA.configs` even in a `gwas_only` project (but no `GEA.snp_sets` there: in `gwas_only` a
  declared GEA set names no predictors, and every mode stops when the run starts).
- **Main outputs:**
  - `GEAxGWAS/tables/pairwise_overlap_table.tsv` and `pairwise_collapsed_snps.tsv`: shared
    SNPs per pair of traits, from `gwas` alone in a `gwas_only` project
  - `GEAxGWAS/plots/miami_combined_K<k>.png`: GEA above, GWAS below (`standard` only)
- The overlapping regions and their genes are computed in the app, on the GEA x GWAS tab.

```{code-block} text
|run_own_gea_x_gwas|
```

## maladaptation

Genomic offset: fits each method under `Maladaptation.methods` (Gradient Forest, geometric
offset, RDA offset) on each candidate SNP set, and projects it onto future climate.

- **Needs:** `gea` with at least one SNP set of 3 or more SNPs (declared under `GEA.snp_sets`,
  or saved from the app's GEA tab; see the
  [FAQ](../faq.md#maladaptation-stops-with-no-curated-snp-sets)), `sNMF.k_best`,
  `Regime.mode: standard`, and the future climate in `Future.*` (it downloads the WorldClim
  projections of `Future.models` the first time).
- Run `climate` first for a spatial Gradient Forest (`spatial_correction: with`, the
  template's choice, or `both`). A spatial Gradient Forest also needs more sites than
  `Climate.predictors` plus one; the mode stops at once otherwise.
- **Main outputs:**
  - `Maladaptation/tables/<method>/<set>_<spatial>/<scenario>/genetic_offset_site.tsv` and
    `genetic_offset_map.tsv`: the offset per site and per map cell, where `<spatial>` is
    `spatial` or `nospatial`
  - `Maladaptation/plots/<method>/<set>_<spatial>/`: offset maps and, for Gradient Forest,
    predictor importance
  - `Maladaptation/tables/snp_sets_summary.tsv`: every stored SNP set and whether this run
    used it

```{code-block} text
|run_own_maladaptation|
```
