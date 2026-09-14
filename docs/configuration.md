# Configuration Reference

All pipeline parameters are defined in `config.yaml`. This document covers every parameter group, GO enrichment details, and the output directory structure.

---

## Parameter Reference

### Project setup

```yaml
project_name: MyProject
Input:
  dir: data/
  vcf: your_genotypes.vcf
  metadata: your_metadata.tsv
  gff: your_annotation.gff3
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `project_name` | Name for output directories (`{name}_results/`, `{name}_logs/`) | — |
| `Input.dir` | Directory containing input files | `"data/"` |
| `Input.vcf` | VCF filename (`.vcf` or `.vcf.gz`) | — |
| `Input.metadata` | TSV with columns: site, sample, latitude, longitude [, traits...] | — |
| `Input.gff` | GFF3 annotation file | — |

CPU cores are auto-detected (`max(1, cpu_count - 2)`). Override with optional `cpu: N`.

---

### Filter — VCF quality filtering

```yaml
Filter:
  maf: 0.05
  snp_miss: 0.1
  sample_miss: 0.5
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Filter.maf` | Minor allele frequency threshold | `0.05` |
| `Filter.snp_miss` | Max SNP missingness (fraction) | `0.1` |
| `Filter.sample_miss` | Max sample missingness (fraction) | `0.5` |

---

### LD — Linkage disequilibrium pruning

```yaml
LD:
  window: 100
  step: 20
  r2: 0.2
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `LD.window` | Sliding window size (kb) | `100` |
| `LD.step` | Step size (kb) | `20` |
| `LD.r2` | R² threshold for pruning | `0.2` |

---

### sNMF — Population structure

```yaml
sNMF:
  k_start: 3
  k_end: 7
  k_best: 5
  repeats: 10
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sNMF.k_start` | Minimum K to test | `3` |
| `sNMF.k_end` | Maximum K to test | `7` |
| `sNMF.k_best` | Selected K (set after reviewing cross-entropy from `mode=prestructure`) | — |
| `sNMF.repeats` | Number of sNMF repetitions per K | `10` |

Ploidy is hardcoded to 2 (diploid). The pipeline does not support haploid organisms.

---

### Map — Geographic extent and climate resolution

```yaml
Map:
  climate_extent: auto
  gap: 0.5
  resolution: 0.5
  zoom_extent: "NULL"
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Map.climate_extent` | Geographic extent: `"auto"` (from sample coords) or `[min_lon, max_lon, min_lat, max_lat]` | `"auto"` |
| `Map.gap` | Buffer around samples (degrees) | `0.5` |
| `Map.resolution` | WorldClim resolution (arc-minutes: 0.5, 2.5, 5, 10) | `0.5` |
| `Map.zoom_extent` | Zoom coordinates: `"xmin,xmax,ymin,ymax"` or `"NULL"` | `"NULL"` |

**Regional zoom**: When `Map.zoom_extent` is set (e.g., `"34.8,35.8,32.0,33.3"`), the pipeline generates both full-extent and zoomed plots for all piemaps and genetic offset maps. Zoomed plots are organized in coordinate-named subdirectories under `zoom/`.

---

### Climate — Bioclimatic variables

```yaml
Climate:
  enabled: true
  predictors: bio_1,bio_12,bio_15
  Varpart:
    response: pcs
    response_var_cutoff: 0.8
    response_max_pcs: 20
    response_min_pcs: 2
    structure_table: qmatrix
    permutations: 999
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Climate.enabled` | Enable climate download and climate-dependent analyses. Set `false` to skip GEA, maladaptation, and gea_x_gwas modes. Structure runs without climate plots/piemaps (only imputation + pop stats). GWAS runs without piemaps. | `true` |
| `Climate.predictors` | Comma-separated WorldClim bioclimatic variables for GEA/GWAS analyses (bio_1 through bio_19). All 19 are shown in the Structure tab piemaps; this list controls which are used as predictors in association. Ignored when `enabled: false`. | — |

**`Climate.Varpart.*`** configures the `mode=climate` variance-partitioning block: `Y ~ X_clim + X_struct + X_geo`, where Y is the response matrix below, X_clim is `Climate.predictors` (GEA tab), X_struct is the structure covariate below, and X_geo is the dbMEM spatial eigenvectors (always site-level, no config — see the dbMEM note below). Reused by the spatial Gradient Forest in `mode=maladaptation`.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Climate.Varpart.response` | Response matrix Y for `varpart()`: `pcs` (LEA genomic PCs — documented scalability adaptation) or `snps` (raw LD-pruned SNP matrix — literal Lasky et al. 2012 methodology, impractical at WGS scale) | `"pcs"` |
| `Climate.Varpart.response_var_cutoff` | Keep the leading genomic PCs that together explain at least this fraction of genetic variance (when `response: pcs`) | `0.8` |
| `Climate.Varpart.response_max_pcs` | Ceiling on the number of response PCs retained regardless of variance cutoff | `20` |
| `Climate.Varpart.response_min_pcs` | Floor on response PCs so varpart always has a multivariate Y | `2` |
| `Climate.Varpart.structure_table` | Structure covariate X_struct: `qmatrix` (sNMF Q at `k_best`, last column dropped — Q rows sum to 1 so all K columns are collinear) or `none` (drops the structure fraction; 2-table climate-vs-geography varpart only). Q is used rather than PCs here because Y is already PCs — an ancestry-PC covariate would make X = Y, a tautology. | `"qmatrix"` |
| `Climate.Varpart.permutations` | `varpart` `anova.cca` permutations per testable fraction. Since 2026-09-13 every block (Y, climate, Q, MEMs) is collapsed to **one row per sampling site** before any fit, so the permutation unit and Ezekiel's n are the sites (`variance_partition.tsv` records `unit`, `n_units`, `df_env = n_sites − 1`); with fewer sites than climate predictors + 2 the status is `insufficient_sites`, and a 3-way partition that would leave no residual df degrades to the climate + structure 2-way (model label says so). The confounding flag clamps negative adjusted fractions at 0 and fires only when the joint climate + geography model is significant (`joint_p` in `climate_confounding.tsv`) | `999` |

The climate–geography confounding diagnostic (flagged when the shared fraction exceeds either unique fraction) is always computed — not a config switch.

`Climate.Varpart.ordir2step_pin` (`0.01`) and `.r2_permutations` (`999`) are **not** config keys — they are fixed constants (`CLIMATE_VP_PIN` / `CLIMATE_VP_R2PERMS` in `common.smk`), the Blanchet double-stopping-rule standard for forward selection.

**dbMEM** (the geography table X_geo) has **no config keys at all** — it is always computed at site level (one point per site, the coordinate centroid of that site's samples, down to n=1 sample), which avoids the pseudoreplication that sample-level dbMEM would introduce (climate/spatial predictors are constant within a site). Below 3 distinct sites, `mode=climate` writes a skip record instead of crashing (surfaced as a warning badge in the Shiny Climate tab, with the suggestion to expand sampling or run Gradient Forest with `Maladaptation.methods.gradient_forest.spatial_correction: without`), rather than exposing a min-sites threshold as a config knob. `dbMEM.spatial_level` / `dbMEM.min_sites` (previously configurable) are removed; a stale `Climate.dbMEM:` block left in an older config file is simply not read.

---

### Population — Population statistics (optional)

```yaml
Population:
  calc_stats: false
  window_size: 10000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Population.calc_stats` | Calculate Tajima's D, Pi diversity, IBD, Mantel test, AMOVA | `false` |
| `Population.window_size` | Window size (bp) for Tajima's D and Pi diversity | `10000` |

---

### LDdecay — LD decay analysis

```yaml
LDdecay:
  group_by: site
  min_samples: 5
  max_distance: 500
  scope: both
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `LDdecay.group_by` | Group samples by: `"site"` or `"cluster"` (sNMF K) | `"site"` |
| `LDdecay.min_samples` | Minimum samples per group to include | `5` |
| `LDdecay.max_distance` | Maximum pairwise SNP distance (kb) to compute | `500` |
| `LDdecay.scope` | Compute: `"genome_wide"`, `"per_chr"`, `"per_pop"`, or `"both"` (genome + pop) | `"both"` |

LD decay is used to determine `GEA.region_distance`: set `region_distance: auto` to use the r²<0.2 distance from LD decay results.

---

### Piemap — Pie map display

```yaml
Piemap:
  alpha: 0.8
  show_labels: true
  label_size: 3
  pie_scale: 0.5
  use_points: false
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Piemap.alpha` | Transparency of pie charts | `0.8` |
| `Piemap.show_labels` | Show site labels | `true` |
| `Piemap.label_size` | Label font size | `3` |
| `Piemap.pie_scale` | Base pie chart size multiplier | `0.5` |
| `Piemap.use_points` | Use points instead of pies (for very many sites) | `false` |

---

### PreGEA — Hyperparameter exploration (optional, `mode=pregea`)

```yaml
PreGEA:
  predictors: bio_1,bio_2,bio_3
  k_offset: 2
  n_pcs_max: 10
  Advanced:
    collinearity_r: 0.7
    vif_max: 10.0
    axis_alpha: 0.05
    permutations: 199
  TransferGuard:
    enabled: no
    lfmm_k: auto
    emmax_n_pcs: auto
```

`mode=pregea` sweeps hyperparameter ladders on the cheap LD-pruned dataset — LFMM-K, EMMAX-#PC, and RDA-Condition()-PC — and writes one recommended value per (method, parameter) to `PreGEA/tables/pregea_recommendations.tsv`, pre-filled into the Shiny GEA tab's method editor. LFMM/EMMAX/RDA always run together (no per-block switches).

| Parameter | Description | Default |
|-----------|-------------|---------|
| `PreGEA.predictors` | Comma-separated predictors for the ladders; falls back to `Climate.predictors` when unset | `Climate.predictors` |
| `PreGEA.k_offset` | LFMM K sweep range = `sNMF.k_best ± k_offset` | `2` |
| `PreGEA.n_pcs_max` | Shared ceiling for the EMMAX/GAPIT #PCs **and** RDA Condition()-PC ladders (both sweep `0..n_pcs_max`) | `10` |
| `PreGEA.Advanced.collinearity_r` | `\|r\|` pre-screen threshold applied before the RDA fit. The correlation matrix is computed **site-wise** (one row per distinct environment), not per sample — climate values repeat for every sample at a site, so a per-sample correlation would weight each site by its sample size and, on an unbalanced design, select a different predictor set. Falls back to per-sample with a `WARNING` if fewer than 3 distinct environments. | `0.7` |
| `PreGEA.Advanced.vif_max` | Post-fit `vif.cca` cutoff; a rung with `max_vif ≥ vif_max` gets `status="vif_exceeded"` and is excluded from recommendations | `10.0` |
| `PreGEA.Advanced.axis_alpha` | RDA constrained-axis significance level | `0.05` |
| `PreGEA.Advanced.permutations` | RDA `anova.cca` permutations (`999` is the production value; the default is a faster testing value) | `199` |
| `PreGEA.TransferGuard.enabled` | Opt-in full-set λ re-check at the recommended hyperparameters — the only PreGEA rule that pulls the full-set imputation chain into this otherwise LD-pruned-only mode | `false` |
| `PreGEA.TransferGuard.lfmm_k` | LFMM K for the guard; `auto` resolves from `pregea_recommendations.tsv` | `"auto"` |
| `PreGEA.TransferGuard.emmax_n_pcs` | EMMAX #PCs for the guard; `auto` resolves from `pregea_recommendations.tsv` | `"auto"` |

Many former `PreGEA.*` fields are now fixed constants rather than config keys (seed = 42, kinship = BN, PreGEA `genomic_control` = FALSE, FDR = 0.1, `min_predictors` = 2, `lambda_tol` = 0.15, `deflation_floor` = 0.90, `ordiR2step` `Pin`/`R2permutations` = 0.01/999). See the README's **Fixed constants (PreGEA / Climate)** section for the full list and rationale.

---

### GEA — Climate association analysis

```yaml
GEA:
  configs:
    - method: EMMAX
      adjust: bonf
      threshold: '0.05'
    - method: LFMM
      adjust: qvalue
      threshold: '0.05'
    - method: BLINK
      adjust: bonf
      threshold: '0.05'
    - method: RDA
      adjust: bonf
      threshold: '0.01'
      params:
        condition_pcs: 3
        axes: auto
  combine_gap: 200000
  region_distance: 2000000
  sig_snp_distance: 10000
  top_regions: 10
  promoter_length: 3000
  scattermore_threshold: 30000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GEA.configs` | List of `{method, adjust, threshold, params}` entries (see below) | — |
| `GEA.combine_gap` | Gap (bp) for merging per-method SNP sets before region clustering | `200000` |
| `GEA.region_distance` | Region merging distance: `auto_per_chromosome` (default, per-chr LD decay), `auto_genome_wide`, or fixed bp integer | `"auto_per_chromosome"` |
| `GEA.region_r2_threshold` | r² level at which LD is considered background; Hill-Weir curve inverted at this value | `0.2` |
| `GEA.region_ld_decay_group` | Sample group from LD decay table for distance computation | `"All"` |
| `GEA.promoter_length` | Promoter length (bp) upstream of gene start for SNP counting | `3000` |
| `GEA.combine_method` | How per-method sig-SNP tables merge into the **static** `selected_snps.tsv` (and what the pipeline derives from it): `All` (union over every configured method), `Overlap` (only SNPs called by ≥2 methods within `snp_clumping_distance`), `MethodOverlap` (as `Overlap`, same predictor required), or a single configured method name. Does not affect the Shiny filter bar's strategy — see **Combine strategy** below | `"All"` |
| `GEA.scattermore_threshold` | SNP count above which scattermore is used for Manhattan background rendering | `30000` |

**configs format**: Each entry specifies:
- `method`: `EMMAX`, `LFMM`, `RDA`, or any GAPIT3 model (`GLM`, `MLM`, `CMLM`, `ECMLM`, `SUPER`, `MLMM`, `FarmCPU`, `BLINK`)
- `adjust`: p-value correction — `bonf` (Bonferroni), `qvalue` (FDR), or a numeric threshold (top-N selection)
- `threshold`: significance threshold after correction
- `params` (optional): per-method hyperparameters. Source of truth is the method registry (`workflow/methods/gea.py`), also surfaced in the Shiny GEA sidebar's method editor. Omit to use registry defaults (existing configs need no changes). Per-method keys:
  - `EMMAX`: `n_pcs` (int, default = sNMF `k_best`), `kinship` (`BN`|`IBS`, default `BN`)
  - `LFMM`: `K` (int, default = sNMF `k_best`), `genomic_control` (bool, default `TRUE` — `lfmm2.test(genomic.control=...)` for the production run; PreGEA's own K-ladder always fits with this OFF so λ_GC stays informative)
  - `GAPIT models`: `n_pcs` (int, default = sNMF `k_best`)
  - `RDA`: `axes` (`auto` or int ≥2), `axis_alpha` (float, default `0.05`), `condition_pcs` (int, default = sNMF `k_best`), `predictor_set` (`auto` or comma-separated subset of `Climate.predictors`), `permutations` (int, default `999`), `fit_mode` (`auto`|`full`|`pruned` — `pruned` is an engineering fallback, not the literature default), `full_fit_max_snps`/`full_fit_max_gb` (size gate for `fit_mode: auto`), `seed` (int, default `42`)

**RDA** is multivariate: it emits one p-value column (`climate_multivariate`, a pseudo-trait standing in for the whole predictor set) rather than one column per predictor, so the Bonferroni denominator stays the full marker count. Per-SNP diagnostics (Mahalanobis distance, q-value, axis loadings, per-predictor correlation, assigned predictor) live in a side table, `GEA/tables/methods/RDA/RDA_candidates_K{k}.tsv`. Model diagnostics (aliased terms, retained axis count, adj. R², GIF λ, marker-envelope status) are in `RDA_diagnostics_K{k}.tsv`; full/by-axis/by-margin ANOVA in `RDA_anova_K{k}.tsv`. Requires `Climate.enabled: yes` and ≥2 predictors.

**Combine strategy** — there are two, and they are separate on purpose:

- **Interactive (set in Shiny, not config)** — the GEA/GWAS filter bar's strategy selector: `All` (union of all methods), `Overlap` (≥2 methods agree within the SNP clumping distance), or `MethodOverlap` (same, plus the two methods must agree on the same predictor). It re-derives from the per-method sig-SNP tables, drives the live region table/rectangles, and is stamped as provenance onto any SNP set you save for maladaptation — so **every genetic-offset run already carries the strategy you picked there**. Persisted per project in `_intermediate/region_params.json`, not in this YAML.
- **Static (config)** — `GEA.combine_method` / `GWAS.combine_method` (above). Governs the pipeline-written `selected_snps.tsv`, the `regions_*`/`genes_*`/enrichment tables Snakemake derives from it, the `mode=gea_x_gwas` pairwise overlap table, and the combined/Miami background superset. Default `All` is a plain union with no family-wise control across methods — with ~10 methods configured it is the least selective option, and it is deliberately the default because these tables are the recall-maximal superset the interactive path curates down from. Set it only if you want the *published static* tables themselves narrowed.

Note that `Overlap` is a ≥2-of-N rule over whatever methods are configured — with 11 methods, ≥2-of-11 is much weaker than the ≥2-of-3 panel used in the MVP benchmark. To get that scheme, configure three methods, not eleven.

---

### GFF — Feature parsing

```yaml
GFF:
  feature: mRNA
  gene_name: Name
  biotype: biotype
  go_field: ontology
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GFF.feature` | GFF feature type to extract genes (`gene`, `mRNA`, etc.) | `"gene"` |
| `GFF.gene_name` | Attribute for gene display name | `"Name"` |
| `GFF.biotype` | Attribute for biotype filtering | `"biotype"` |
| `GFF.go_field` | GFF attribute containing GO terms (or `"NULL"` to skip enrichment) | `"NULL"` |

---

### Enrichment — GO enrichment visualization

```yaml
Enrichment:
  top_terms: 20
  plot_width: 12
  plot_height: 10
  cnet_label: gene_id
  top_plot_regions: 5
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Enrichment.top_terms` | Max GO terms to show in plots | `20` |
| `Enrichment.plot_width` | Plot width (inches) | `12` |
| `Enrichment.plot_height` | Plot height (inches) | `10` |
| `Enrichment.cnet_label` | Gene label in cnetplot: `"gene_id"`, `"Name"`, or `"description"` | `"gene_id"` |
| `Enrichment.top_plot_regions` | Top regions per trait to generate plots for (`0` = all) | `5` |

GO enrichment is computed on-demand in the Shiny dashboard when the user clicks a region.

---

### Future — Climate projections

```yaml
Future:
  ssp: ssp370
  year: 2061-2080
  models:
    - ACCESS-CM2
    - CNRM-CM6-1
    - MIROC6
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Future.ssp` | Shared Socioeconomic Pathway scenario (`ssp126`, `ssp245`, `ssp370`, `ssp585`) | `"ssp370"` |
| `Future.year` | Future time period (`2021-2040`, `2041-2060`, `2061-2080`, `2081-2100`) | `"2061-2080"` |
| `Future.models` | List of CMIP6 models to ensemble-average | — |

---

### Maladaptation — Gradient Forest

```yaml
Maladaptation:
  methods:
    gradient_forest:
      ntree: 500
      cor_threshold: 0.5
      spatial_correction: both
      random_model: false
      extrap: true
  snp_sets: all
```

Gradient Forest settings are nested under `Maladaptation.methods.gradient_forest`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ntree` | Number of trees in Gradient Forest | `500` |
| `cor_threshold` | Correlation threshold for variable selection | `0.5` |
| `spatial_correction` | Include PCNM spatial variables: `"with"`, `"without"`, or `"both"` | `"both"` |
| `random_model` | Also run a neutral (random SNP) GF model for comparison. When no random SNP has a positive R² the model is written as an `empty_forest` sentinel (the null behaving as a null) and the importance/cumimp plots show the adaptive model alone | `false` |
| `extrap` | How `predict.gradientForest()` projects the offset onto raster cells whose climate lies **outside the training range**: `true` = linear continuation of each turnover function at its boundary slope (the gradientForest default, and what every published GF-offset code base runs by omission); `false` = the predictor is clamped to the range limit, so the turnover function plateaus (conservative, underestimates offset in novel climate). The manual calls extrapolation "an experimental feature"; Rellstab et al. 2021 ask that novel-climate predictions be flagged, which `gradient_forest_diagnostics.tsv` does (per-predictor training envelope, % of present / future cells outside it). Added 2026-09-13 (audit SC1) | `true` |

**`Maladaptation.snp_sets`** (sibling of `methods`, not nested under `gradient_forest`): `"all"` or a list of SNP set names to run maladaptation on. SNP sets are named and saved from the GEA tab's "Save SNP set for maladaptation" action in Shiny — there is no pipeline-side `run_label`/`combine_method`/`combine_gap` selection anymore; SNP selection happens interactively before the set is saved.

---

### GWAS — Phenotype association

```yaml
GWAS:
  configs:
    - method: EMMAX
      adjust: bonf
      threshold: '0.05'
  missing_strategy: DROP
  traits: height,flowering_time
  region_distance: 2000000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GWAS.configs` | Same format as `GEA.configs` | — |
| `GWAS.missing_strategy` | How to handle missing values: `DROP`, `MEAN`, or `MEDIAN` | `"DROP"` |
| `GWAS.traits` | Comma-separated trait column names to include (omit to use all cols 5+) | all |
| `GWAS.region_distance` | Region merging distance (defaults to `GEA.region_distance`) | inherited |
| `GWAS.region_r2_threshold` | r² background threshold (defaults to `GEA.region_r2_threshold`) | inherited |
| `GWAS.region_ld_decay_group` | LD decay group (defaults to `GEA.region_ld_decay_group`) | inherited |
| `GWAS.combine_gap` | Gap for merging SNP sets (defaults to `GEA.combine_gap`) | inherited |
| `GWAS.combine_method` | Same as `GEA.combine_method`, for the static GWAS tables. **Not** inherited from GEA — the two sources have different method sets, so a single-method name valid for one is invalid for the other | `"All"` |
| `GWAS.promoter_length` | Promoter length (defaults to `GEA.promoter_length`) | inherited |

**Missing value strategies**:
- **DROP** — Per-trait sample subsetting. Samples missing a trait value are excluded from that trait's analysis. Separate VCF, TPED, and kinship are computed per trait.
- **MEAN** — Missing values replaced with trait mean. All samples used together.
- **MEDIAN** — Missing values replaced with trait median. All samples used together.

---

### GEAxGWAS — GEA + GWAS overlap

```yaml
GEAxGWAS:
  pairwise:
    window_size: 500000
    min_snps: 1
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GEAxGWAS.pairwise.window_size` | Window (bp) around SNPs for window-based overlap scoring | `500000` |
| `GEAxGWAS.pairwise.min_snps` | Minimum overlapping SNPs to report a pair | `1` |

`region_distance` for overlap region computation is set interactively in the Shiny dashboard (defaults to `max(GEA.region_distance, GWAS.region_distance)`).

---

## GO Enrichment Analysis

When a GFF file with GO annotations is provided (`GFF.go_field` parameter), the pipeline can perform functional enrichment analysis on genes associated with significant regions.

### How it Works

Enrichment is computed **on-demand** in the Shiny dashboard when a region is selected:
1. **Gene finding** — Genes overlapping the region boundaries are identified from the GFF.
2. **Enrichment testing** — Per-region GO enrichment via clusterProfiler (Fisher's exact test).
3. **Visualization** — Three plot types generated: dotplot, emapplot, cnetplot.
4. **Persistence** — Results saved to `{MODULE}/plots/enrichment/` and `{MODULE}/tables/enrichment/` for instant reuse.

### GFF Requirements

Your GFF file must include GO terms in the attributes column. Specify the attribute via `GFF.go_field`:

```gff
chr1  AUGUSTUS  gene  1000  5000  .  +  .  ID=gene001;Name=GENE1;ontology=GO:0009414,GO:0042538
```

```yaml
GFF:
  go_field: ontology    # attribute containing GO terms
```

Set `GFF.go_field: "NULL"` to skip enrichment.

### Visualization Types

**Dotplot** — Enriched GO terms ranked by significance, with gene ratio and adjusted p-value. Works with ≥1 enriched term.

**Emapplot** — GO term similarity network. Nodes = GO terms sized by gene count. Requires ≥2 enriched terms.

**Cnetplot** — Gene-concept bipartite network. Gene labels configurable via `Enrichment.cnet_label` (`gene_id`, `Name`, or `description`). Requires ≥2 enriched terms.

---

## Output Structure

```
{project_name}_results/
├── Processing/
│   └── tables/                            # metadata.tsv, sample_missing_stats.tsv
│
├── PreStructure/
│   ├── plots/                             # pca.png, tracy_widom.png, cross_entropy_K{s}-{e}.png
│   │   └── K{k}/                          # structure_K{k}.png, pca_structure_K{k}.png
│   └── tables/
│       └── K{k}/                          # clusters_K{k}.tsv (Q-matrices)
│
├── climate/
│   ├── plots/                             # density_plot_present, correlation_heatmap, density_plot_future_*
│   ├── plots/{spatial,varpart}/           # dbmem_screeplot, varpart_venn, px_barplot (mode=climate)
│   ├── tables/present/                    # climate_present_{all,site,site_scaled}.tsv, climate_invariant_predictors.tsv
│   ├── tables/spatial/                    # dbmem_vectors, dbmem_diagnostics (mode=climate)
│   ├── tables/varpart/                    # variance_partition (unit/n_units/df_env), climate_confounding (joint_p/n_sites), px_per_variable, dbmem_selected (mode=climate; fitted on one row per SITE since 2026-09-13)
│   ├── tables/future/                     # climate_future_year{Y}_ssp{S}_{site,all}.tsv
│   └── rasters/{present,future}/          # WorldClim .tif rasters
│
├── Structure/
│   ├── plots/piemap/                      # piemap_{bio}.png/svg + zoom/
│   ├── plots/piemap/{tajima_d,pi_diversity}/   # trait-scaled piemaps (optional)
│   ├── plots/pop_stats/                   # mantel_test, amova (optional)
│   ├── plots/ld_decay/                    # ld_decay_genome_wide, per_chr, per_pop
│   └── tables/pop_stats/                  # tajima_d_by_pop, pi_diversity_by_pop, amova
│
├── GEA/
│   ├── GAPIT_native_output/{model}/       # raw GAPIT output files
│   ├── plots/manhattan/
│   │   ├── {method}/                      # manhattan_{trait}_K{k}_{adjust}.png/svg
│   │   │                                  # qq_{trait}_K{k}_{adjust}.png/svg
│   │   └── combined/                      # manhattan_combined_K{k}.png/svg
│   ├── plots/enrichment/{trait}/          # region_{id}_dotplot/emapplot/cnetplot (Shiny on-demand)
│   └── tables/
│       ├── methods/{method}/              # {method}_pvalues_K{k}.tsv, _sig_snps_{adjust}.tsv
│       ├── selected_snps.tsv, selected_snps_per_trait.tsv, regions_per_trait.tsv, regions_combined.tsv
│       ├── genes_per_region.tsv, genes_per_region_collapsed.tsv, genes_combined.tsv
│       └── enrichment/{trait}/            # GO enrichment TSVs (Shiny on-demand)
│
├── GWAS/                                  # same structure as GEA/
│   └── plots/piemap/                      # phenomap_{trait}.png/svg + zoom/
│
├── GEAxGWAS/
│   ├── plots/miami_combined_K{k}.{png,svg,_background.png,_coords.json}
│   └── tables/pairwise_{overlap_table,collapsed_snps}.tsv
│
├── Maladaptation/
│   ├── plots/gradient_forest/{run_label}_{spatial_tag}/
│   │   └── cumulative_importance, overall_importance, genetic_offset_piemap[_{tajima_d,pi_diversity}]
│   └── tables/gradient_forest/{run_label}_{spatial_tag}/
│       └── genetic_offset_{map,site}.tsv
│
├── pipeline_summary.tsv                   # all modes append here
│
├── _work/maf{}_miss{}_smiss{}/ld{}_win{}_step{}/    # parameterized intermediates
│   ├── filtered.vcf, pruned.vcf
│   ├── {name}.geno, .lfmm, .vcfsnp
│   └── {name}.snmfProject/
│
└── _intermediate/                         # internal: samples/, annotation/, {method}/, qs_cache/
    └── region_params.json                 # per-region Shiny parameter persistence

{project_name}_logs/{module}/              # per-module log directories
```

**Path parameterization**:
- `_work/maf0.05_miss0.1_smiss0.5/` — VCF filtered with these thresholds
- `_work/.../ld0.2_win100_step20/` — LD-pruned with R²=0.2, 100kb window, 20kb step
- `Maladaptation/.../my_run_spatial/` — GF run_label + spatial_correction tag
