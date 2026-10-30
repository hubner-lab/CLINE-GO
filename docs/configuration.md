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
  imputed_maf: null
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sNMF.k_start` | Minimum K to test | `3` |
| `sNMF.k_end` | Maximum K to test | `7` |
| `sNMF.k_best` | Selected K (set after reviewing cross-entropy from `mode=prestructure`) | — |
| `sNMF.repeats` | Number of sNMF repetitions per K | `10` |
| `sNMF.imputed_maf` | MAF floor on the sNMF-**imputed** matrices (LD-pruned and full), i.e. the SNP set every GEA/GWAS method tests. Blank/`null` = same as `Filter.maf`; `0` = filter off; must be in [0, 0.5) | `null` |

Ploidy is hardcoded to 2 (diploid; `PLOIDY` in `workflow/rules/common.smk`). There is no config key: a
`sNMF.ploidy` other than 2 is refused at parse time. Haploid and polyploid organisms are not supported.

**Diploid likelihood on a panel without heterozygotes.** sNMF counts each individual as two independent
allele draws. On a panel with (almost) no heterozygous calls — complete selfers, inbred lines, or a VCF
whose heterozygous calls were masked to missing at calling — every homozygote is counted twice, so Q and
the cross-entropy choice of K come out more confident than the data allow. Nothing in the run warns about
it: read `Processing/tables/sample_heterozygosity.tsv` (or the het-vs-missingness QC plot) before
trusting K. (Trifolium86Main has zero heterozygous calls panel-wide.)

**Why a second MAF filter.** `Filter.maf` runs as plink `--maf` on *observed* calls. Imputation then
fills every missing call, mostly with the major allele, so a SNP whose minor allele sits in one small
site can pass `Filter.maf` and be tested well below it (Trifolium86Main: 16.7% of tested SNPs had
imputed MAF < 0.1). `filter_imputed_maf_ld` (structure mode) and `filter_imputed_maf_full` (gea/gwas)
recompute MAF on the imputed matrix and keep SNPs with MAF ≥ the threshold. Before/after counts go to
`pipeline_summary.tsv` (`imputed_maf_threshold`, `snps_before_imputed_maf`, `snps_after_imputed_maf`;
step `structure` = LD-pruned set, `gea`/`gwas` = full set) and to the GEA/GWAS tab KPI strip.

**Limitation.** MAF is computed over the full post-filter cohort. Steps that test a sample *subset* —
the climate-valid subset (samples with NA climate dropped) and GWAS `missing_strategy: DROP` per-trait
subsets — can see a tested MAF slightly below the threshold.

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
| `Map.climate_extent` | Geographic extent: `"auto"` (from sample coords) or one string `"min_lon,max_lon,min_lat,max_lat"` (not a YAML list) | `"auto"` |
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

AMOVA computes its distance on at most 10,000 SNPs of the LD-pruned, imputed VCF, drawn with seed 42.
`amova.tsv` ends with the rows that make its percentages reproducible: `n_snps_vcf`, `snp_cap`,
`n_snps_sampled`, `n_loci_used` and `seed`.

---

### LDdecay — LD decay analysis

```yaml
LDdecay:
  group_by: site
  min_samples: 10
  max_distance: 500
  scope: both
  max_pairs: 20000000
  background_markers_per_chr: 50
  background_quantile: 0.95
  r2_threshold: 0.2
  loess_span: 0.3
  min_signal_ratio: 2
  seed: 1
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `LDdecay.group_by` | Group samples by: `"site"` or `"cluster"` (sNMF K) | `"site"` |
| `LDdecay.min_samples` | Minimum samples per group to include | `10` |
| `LDdecay.max_distance` | Maximum pairwise SNP distance (kb) to compute | `500` |
| `LDdecay.scope` | `"genome_wide"`, `"per_chromosome"`, or `"both"` | `"both"` |
| `LDdecay.max_pairs` | SNP pairs PLINK computes per group; above it, a seeded sample of anchor SNPs is paired with every SNP of their window (all distances thinned evenly) | `2e7` |
| `LDdecay.background_markers_per_chr` | Markers sampled per chromosome for the unlinked (between-chromosome) pairs | `50` |
| `LDdecay.background_quantile` | Quantile of sqrt(r²) of unlinked pairs whose square is the critical r² | `0.95` |
| `LDdecay.r2_threshold` | r² level of the `r2_threshold` decay rule | `0.2` |
| `LDdecay.loess_span` | Span of the LOESS curve through the binned mean r² | `0.3` |
| `LDdecay.min_signal_ratio` | A curve is reliable when its short-range r² is at least this many times the critical r² of unlinked markers (must be > 1); an `auto_*` distance stops on an unreliable curve | `2` |
| `LDdecay.seed` | Seed for the anchor SNPs and the unlinked markers | `1` |
| `LDdecay.exclude_chromosomes` | Optional: contigs left out of the decay curve and its unlinked background, such as organellar or sex contigs, which do not recombine like autosomes. One comma-separated string (`"X,Y,MT,Pt"`, what the app writes) or a list; case, a leading `chr` and M/MT are ignored, and a name with no contig is reported in the log. A YAML boolean is refused: in a list, write `Y` in quotes | none |

r² is PLINK 1.9's `--r2` (squared allele-count correlation, missing calls excluded pairwise) on the filtered, unimputed VCF, for pairs on one chromosome (the decay curve) and on different chromosomes (its background). The mean r² is binned by distance (10 bp below 100 bp, log-spaced beyond) and smoothed by LOESS, and three distances are read off the curve per group and scope (`Structure/tables/ld_decay_half_distances.tsv`):

- `decay_critical_bp` — where the curve meets the critical r² of unlinked markers (Breseghello & Sorrells 2006; Laidò et al. 2014);
- `decay_half_bp` — where r² falls to half its short-range value (Zhao et al. 2011; Vos et al. 2017);
- `decay_r2thr_bp` — where r² falls to `LDdecay.r2_threshold`.

Every row also gets `ld_status` (`reliable` / `unreliable`) and `ld_status_reason`. A curve is unreliable when its short-range r² is less than `LDdecay.min_signal_ratio` (default 2) times the critical r² of unlinked markers: below 2 the unlinked background lies above half the short-range LD, so the curve meets noise before it has halved and every distance read off it is set by the background rather than by decay. A small sample (unlinked r² grows roughly as 1/n) and a few inbred, strongly differentiated lineages both compress the ratio. A curve is also unreliable when it has too few distance bins, or no background (a single chromosome). Measured 2026-10-27 (genome-wide, group All): Arabidopsis 1001G (n 1,049) 14.5; wild barley WBDC (580) 4.7; barley GBS panels (303, 582) 3.9 and 4.4; 57-plant subsamples of Arabidopsis and WBDC 2.7–5.4; 600 SS-Clines simulations (1,000) 4.8–36.5; Trifolium (57 plants, five selfing lineages) 1.75, the one unreliable curve.

An `auto_genome_wide` / `auto_per_chromosome` clumping or region distance is read from one of these columns, chosen by `GEA.clumping_ld_rule` (below). On an unreliable genome-wide curve it **stops** with the reason, whatever the rule: declare `snp_clumping_distance` (or `region_clumping_distance`, or a SNP set's `clumping_distance`) in bp instead. Under `auto_per_chromosome` an unreliable chromosome takes the genome-wide distance (with a warning) when the genome-wide curve is reliable, and stops otherwise. The WZA `window_size: auto_*` keeps its fallback window on an unreliable curve. A table written before the check (no `ld_status`) is read as before, with a notice to re-run `mode=structure`. PLINK runs with `--chr-set 95 --allow-extra-chr`: numeric chromosomes up to 95 stay autosomes, and `X`, `MT` and other names (`Pt`, scaffolds) enter the curve and the background with an autosome's r². `Y` gives no pairs, because a VCF carries no sex and PLINK treats non-males' Y calls as missing. To leave such contigs out, list them in `LDdecay.exclude_chromosomes`.

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
      adjust: qval
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
  snp_clumping_distance: auto_genome_wide
  clumping_ld_rule: critical
  top_regions: 10
  promoter_length: 3000
  scattermore_threshold: 30000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GEA.configs` | List of `{method, adjust, threshold, params}` entries (see below) | — |
| `GEA.snp_clumping_distance` | ONE distance for the agreement between methods, the cross-trait overlap annotation and the region merge: fixed bp, or `auto_genome_wide` / `auto_per_chromosome` read from the LD-decay table (legacy names `region_distance` / `combine_gap` migrate to it). Unset or null, it is `auto_genome_wide` and the parse says so on stderr; until 2026-10-26 the default was a fixed `100000`. Under `LDdecay.scope: per_chromosome` the default is refused: set bp or `auto_per_chromosome` | `auto_genome_wide` |
| `GEA.region_clumping_distance` | Advanced: a separate distance for the region merge and flank only (bp or `auto_*`); unset = `snp_clumping_distance` | unset |
| `GEA.clumping_ld_rule` | Which LD-decay distance an `auto_*` distance takes: `critical`, `half_decay` or `r2_threshold` (see LDdecay). `auto_per_chromosome` needs `LDdecay.scope` per_chromosome or both, `auto_genome_wide` (and a declared SNP set's auto gap) genome_wide or both. A leftover `clumping_r2_threshold` no longer picks the rule (a warning says so) | `"critical"` |
| `GEA.ld_decay_group` | Sample group from the LD-decay table for the distance (legacy `region_ld_decay_group`) | `"All"` |
| `GEA.promoter_length` | Promoter length (bp) upstream of gene start for SNP counting | `3000` |
| `GEA.combine_method` | How per-method sig-SNP tables merge into the **static** `selected_snps.tsv` (and what the pipeline derives from it): `All` (union over every configured method), `Overlap` (only SNPs called by ≥2 methods within `snp_clumping_distance`), `MethodOverlap` (as `Overlap`, same predictor required), or a single configured method name. Does not affect the Shiny filter bar's strategy — see **Combine strategy** below | `"All"` |
| `GEA.scattermore_threshold` | SNP count above which scattermore is used for Manhattan background rendering | `30000` |

**configs format**: Each entry specifies:
- `method`: `EMMAX`, `LFMM`, `RDA`, or any GAPIT3 model (`GLM`, `MLM`, `CMLM`, `ECMLM`, `SUPER`, `MLMM`, `FarmCPU`, `BLINK`)
- `adjust`: the significance rule — `bonf` (Bonferroni α), `qval` (FDR q target), `top` (the N smallest p-values per predictor), `toppct` (the top share of SNPs per predictor, in percent: `top` with N = ceiling(share × the table's rows), NA rows included, so methods tested on the same SNPs get the same N), or `custom` (a raw p cutoff)
- `threshold`: the rule's value — α or q in (0, 1], a SNP count ≥ 1 for `top`, a percentage in (0, 100] for `toppct`, a p > 0 for `custom`
- `params` (optional): per-method hyperparameters. Source of truth is the method registry (`workflow/methods/gea.py`), also surfaced in the Shiny GEA sidebar's method editor. Omit to use registry defaults (existing configs need no changes). Per-method keys:
  - `EMMAX`: `n_pcs` (int, default = sNMF `k_best`), `kinship` (`BN`|`IBS`, default `BN`)
  - `LFMM`: `K` (int, default = sNMF `k_best`), `genomic_control` (bool, default `TRUE` — `lfmm2.test(genomic.control=...)` for the production run; PreGEA's own K-ladder always fits with this OFF so λ_GC stays informative)
  - `GAPIT models`: `n_pcs` (int, default = sNMF `k_best`)
  - `RDA`: `axes` (`auto` or int ≥2), `axis_alpha` (float, default `0.05`), `condition_pcs` (int, default = sNMF `k_best`), `predictor_set` (`auto` or comma-separated subset of `Climate.predictors`), `permutations` (int, default `999`), `fit_mode` (`auto`|`full`|`pruned` — `pruned` is an engineering fallback, not the literature default), `full_fit_max_snps`/`full_fit_max_gb` (size gate for `fit_mode: auto`), `seed` (int, default `42`)

**RDA** is multivariate: it emits one p-value column (`climate_multivariate`, a pseudo-trait standing in for the whole predictor set) rather than one column per predictor, so the Bonferroni denominator stays the full marker count. Per-SNP diagnostics (Mahalanobis distance, q-value, axis loadings, per-predictor correlation, assigned predictor) live in a side table, `GEA/tables/methods/RDA/RDA_candidates_K{k}.tsv`. Model diagnostics (aliased terms, retained axis count, adj. R², GIF λ, marker-envelope status) are in `RDA_diagnostics_K{k}.tsv`; full/by-axis/by-margin ANOVA in `RDA_anova_K{k}.tsv`. Requires `Climate.enabled: yes` and ≥2 predictors.

**Combine strategy** — there are two, and they are separate on purpose:

- **Interactive (set in Shiny, or declared per SNP set)** — the GEA/GWAS filter bar's strategy selector: `All` (union of all methods), `Overlap` (≥2 methods agree within the SNP clumping distance), or `MethodOverlap` (same, plus the two methods must agree on the same predictor). It re-derives from the per-method sig-SNP tables, drives the live region table/rectangles, and is stamped as provenance onto any SNP set you save for maladaptation — so **every genetic-offset run already carries the strategy you picked there**. Persisted per project in `_intermediate/region_params.json`. The same selection, strategy included, can be written into config as a declared SNP set (`GEA.snp_sets`, below), which runs the app's own selection code.
- **Static (config)** — `GEA.combine_method` / `GWAS.combine_method` (above). Governs the pipeline-written `selected_snps.tsv`, the `regions_*`/`genes_*`/enrichment tables Snakemake derives from it, the `mode=gea_x_gwas` pairwise overlap table, and the combined/Miami background superset. Default `All` is a plain union with no family-wise control across methods — with ~10 methods configured it is the least selective option, and it is deliberately the default because these tables are the recall-maximal superset the interactive path curates down from. Set it only if you want the *published static* tables themselves narrowed.

Note that `Overlap` is a ≥2-of-N rule over whatever methods are configured — with 11 methods, ≥2-of-11 is much weaker than the ≥2-of-3 panel used in the MVP benchmark. To get that scheme, configure three methods, not eleven.

**Declared SNP sets (`GEA.snp_sets`, mirrored by `GWAS.snp_sets`)** — the config twin of the GEA tab's "Save SNP set for maladaptation". Each entry names one curated SNP set. `mode=gea` (or `mode=gwas`) builds it into `_intermediate/snp_sets/gea.<name>/` (`gwas.<name>/`), where `mode=maladaptation` picks it up exactly like a set saved from the app. The selection is the app's own code (`scripts/R/lib/snp_set_select.R`), so a declared set reproduces a saved one SNP for SNP. Omitted entirely, nothing changes.

```yaml
GEA:
  snp_sets:
    - name: default                # each method's own GEA.configs rule, every column, Union
    - name: strict
      threshold_type: bonf
      threshold_value: 0.01        # every cell; RDA keeps its registry rule (bonf 0.01)
    - name: tuned
      threshold_type: bonf
      threshold_value: 0.05
      methods: [EMMAX, LFMM]       # omitted = every GEA.configs method
      traits: [bio_1, bio_12]      # omitted = every column those methods emit
      overrides:                   # applied method -> trait -> trait::method
        LFMM: {value: 0.01}        # no type: keeps the type of the rule it overrides
        bio_12: {type: custom, value: 1.0e-5}
        "bio_1::EMMAX": {type: qval, value: 0.1}
      strategy: Cross-method       # Union (default) | Cross-method | Cross-method per-trait | Cross-method >=K | one method
      clumping_distance: 100000    # Cross-method window (bp); default GEA.snp_clumping_distance
    - name: picked
      pairs: ["bio_1::LFMM", "bio_12::EMMAX", "climate_multivariate::RDA"]
```

| Key | Meaning | Default |
|-----|---------|---------|
| `name` | `^[A-Za-z0-9_]+$`, unique within the group. The `gea.` / `gwas.` prefix is added for you, and the app refuses to save a set under either prefix | required |
| `threshold_type` / `threshold_value` | One rule for every cell (a cell is one `trait::method` p-value column). Types: `bonf`, `qval`, `top`, `toppct`, `custom` (see `adjust` above). Type alone takes that type's default value (`0.05` / `0.05` / `100` / `1` / `1e-5`); value alone takes the type of the **first** `configs` entry, as the app's master rule does | neither: each method's own `configs` rule |
| `overrides` | Per-cell exceptions keyed by a method, a trait, or a quoted `"trait::method"`. Applied in that order, so the most specific wins. A rule gives `type`, `value` or both: a missing `value` takes the type's default, a missing `type` keeps the type it overrides. A bare key that is both a trait and a method name is refused | none |
| `methods`, `traits` | Select the `traits × methods` rectangle | all `configs` methods × every column they emit |
| `pairs` | An explicit `trait::method` list instead of a rectangle. Exclusive with `methods`/`traits` | — |
| `strategy` | `Union` (`All`, `Sum`), `Cross-method` (`Overlap`), `Cross-method per-trait` (`MethodOverlap`, `PairOverlap`), `Cross-method >=K`, or one method name from the set's own pairs. `Cross-method >=K` (k-of-n agreement, declared sets only, written exactly so) keeps a SNP when at least K of the set's methods call within `clumping_distance` of it, its own call included; methods are counted, not cells, and calls do not chain. `2 <= K <=` the set's method count; `Cross-method >=2` selects what `Cross-method` does | `Union` |
| `clumping_distance` | Window within which two methods' calls count as agreeing. Used only by the Cross-method strategies. A positive integer of bp, or `auto_genome_wide`, read from the genome-wide LD-decay row under the group's `clumping_ld_rule` and `ld_decay_group` (`auto_per_chromosome` is read genome-wide: the gap is one number) | the group's `snp_clumping_distance` (bp or auto) |
| `regime` | Only `snp`. `wza` is refused (not declarable yet) | `snp` |

The rules a set resolves to, and how they differ from the app's sidebar:

- **No threshold declared** — every cell, RDA included, uses its method's own `configs` rule. That makes `gea.default` the same rule set as the static `selected_snps.tsv`.
- **Any threshold declared** — a multivariate method (RDA) keeps its registry rule (`bonf 0.01`, `workflow/methods/gea.py`), as in the app. An `RDA:` override lifts that pin.
- **RDA emits one column, `climate_multivariate`**. With RDA configured, `traits: [bio_1]` and no `methods:` is refused, because RDA has no `bio_1` column. List `methods` without RDA, or add `climate_multivariate` to `traits`.
- **Everything is checked at parse time**: an unknown key, method, trait or override key; a duplicate name or pair; `pairs` together with `traits`/`methods`; a value outside its type's range. A typo fails `snakemake -n`, not the build after every GEA job has run.
- **A set that selects fewer than 3 SNPs is not an error.** It is written as is (zero SNPs: a header-only TSV with `n_snps: 0` in its record), and `mode=maladaptation` skips it under `snp_sets: all`, because the RDA offset stops below 3 SNPs (`MIN_OFFSET_SNPS`).
- **A one-item list may be a bare string** (`traits: bio_1`). The app's YAML writer collapses `[bio_1]` that way when it re-saves a config.

Each set directory holds `selected_snps.tsv` (`SNPID chr pos min_pvalue`), `declaration.json` (the resolved inputs above, rule input of the build) and `snp_set.yaml`, the reproducibility record (declaration + `k_best`, `n_snps`, `created`). A set saved from the app writes the same record, plus `config_snippet`.

**From the app to config**: saving a set in the GEA tab opens a dialog showing the `GEA.snp_sets` entry that reproduces it, with a copy button. The Maladaptation tab's set details show the same block. Paste it under `GEA.snp_sets`, run `mode=gea`, and `_intermediate/snp_sets/gea.<name>/selected_snps.tsv` matches the saved set. A set saved before records existed shows a block rebuilt from its manifest, with a note that its pair selection is unknown.

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
  predictors: auto
```

**`Maladaptation.predictors`** (sibling of `methods`): the climate predictors every offset model is fitted on — Gradient Forest (adaptive, random, frequency refit, sensitivity check, offset, cumulative-importance plot), geometric offset and RDA offset. `auto` (or unset/blank) = `Climate.predictors`, so an existing config is unchanged. Otherwise a comma string or a list that must be a subset of `Climate.predictors` (the climate tables are built for those columns); anything else fails at parse time. Default `auto`.

It exists because collinearity matters differently by use. The univariate scans (EMMAX, LFMM) test one predictor at a time, so a large correlated `Climate.predictors` does them no harm. The offset models fit all predictors jointly on one row per sample, and climate only varies among sites. With more predictors than sites − 1 the RDA offset aliases the surplus by **column order**: on 16 predictors and 11 sites it silently dropped bio_11, bio_12, bio_13, bio_15, bio_16, bio_19, which is every precipitation predictor. Gradient Forest splits importance across correlated predictors. Pass a non-redundant subset here (e.g. the GEA RDA's `predictor_set`, screened at site level for |r| and VIF). `density_plot_future` (climate characterization) keeps `Climate.predictors`. The predictors are not part of the output path: changing the key reruns the maladaptation rules in place.

Gradient Forest settings are nested under `Maladaptation.methods.gradient_forest`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ntree` | Number of trees in Gradient Forest | `500` |
| `cor_threshold` | Correlation threshold for variable selection | `0.5` |
| `spatial_correction` | Include PCNM spatial variables: `"with"`, `"without"`, or `"both"` | `"both"` |
| `random_model` | Also run a neutral (random SNP) GF model for comparison: as many SNPs as the set, drawn from the SNPs outside it (fewer, with a warning, when fewer lie outside). When no random SNP has a positive R² the model is written as an `empty_forest` sentinel (the null behaving as a null) and the importance/cumimp plots show the adaptive model alone | `false` |
| `extrap` | How `predict.gradientForest()` projects the offset onto raster cells whose climate lies **outside the training range**: `true` = linear continuation of each turnover function at its boundary slope (the gradientForest default, and what every published GF-offset code base runs by omission); `false` = the predictor is clamped to the range limit, so the turnover function plateaus (conservative, underestimates offset in novel climate). The manual calls extrapolation "an experimental feature"; Rellstab et al. 2021 ask that novel-climate predictions be flagged, which `gradient_forest_diagnostics.tsv` does (per-predictor training envelope, % of present / future cells outside it). Added 2026-09-13 (audit SC1) | `true` |

**`Maladaptation.snp_sets`** (sibling of `methods`, not nested under `gradient_forest`): `"all"` or a list of SNP set names to run maladaptation on. A set comes from one of two places, and both land in the same store, `_intermediate/snp_sets/<name>/`:

- **Saved in the app** — the GEA tab's "Save SNP set for maladaptation". Any name without a `gea.`/`gwas.` prefix.
- **Declared in config** — `GEA.snp_sets` / `GWAS.snp_sets` (see the GEA section), named `gea.<name>` / `gwas.<name>`. Built at the end of `mode=gea` / `mode=gwas`, never by `mode=maladaptation` itself.

How the names resolve at parse time:

- `all` = every set in the store plus every declared name.
- A declared set that has not been built yet raises, telling you to run `mode=gea` (or `mode=gwas`) first.
- A set with fewer than 3 SNPs (`MIN_OFFSET_SNPS`, the RDA offset's floor; zero included) is skipped under `all` with an `INFO` line that lists the counts, and raises when named explicitly. The floor is the same whatever offset methods run, and the app's set picker offers only sets above it.
- A directory left behind by a declaration you have since removed stays in `all` until you delete it (from the app's set list or by hand). Deleting a declared set's directory blocks `all` until `mode=gea` rebuilds it.

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
  snp_clumping_distance: auto_genome_wide
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GWAS.configs` | Same format as `GEA.configs` | — |
| `GWAS.missing_strategy` | How to handle missing values: `DROP`, `MEAN`, or `MEDIAN` | `"DROP"` |
| `GWAS.traits` | Comma-separated trait column names to include (omit to use all cols 5+) | all |
| `GWAS.snp_clumping_distance` | As `GEA.snp_clumping_distance` (defaults to it) | inherited |
| `GWAS.region_clumping_distance` | As `GEA.region_clumping_distance` (defaults to it) | inherited |
| `GWAS.clumping_ld_rule` | As `GEA.clumping_ld_rule` (defaults to it) | inherited |
| `GWAS.ld_decay_group` | LD decay group (defaults to `GEA.ld_decay_group`) | inherited |
| `GWAS.combine_method` | Same as `GEA.combine_method`, for the static GWAS tables. **Not** inherited from GEA — the two sources have different method sets, so a single-method name valid for one is invalid for the other | `"All"` |
| `GWAS.snp_sets` | Same grammar as `GEA.snp_sets`, over `GWAS.configs` and the trait columns; sets are named `gwas.<name>` and built by `mode=gwas`. Thresholds, methods and clumping distance come from the GWAS group, never from GEA | none |
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
│   └── tables/                            # tracy_widom.tsv, cross_entropy_K{s}-{e}.tsv (the two plots' numbers)
│       └── K{k}/                          # clusters_K{k}.tsv (Q-matrices)
│
├── climate/
│   ├── plots/                             # density_plot_present, correlation_heatmap, density_plot_future_*
│   ├── plots/{spatial,varpart}/           # dbmem_screeplot, varpart_venn, px_barplot (mode=climate)
│   ├── tables/present/                    # climate_present_{all,site,site_scaled}.tsv, climate_invariant_predictors.tsv
│   │                                      # (written by mode=structure as well), climate_correlation.tsv
│   ├── tables/spatial/                    # dbmem_vectors, dbmem_diagnostics (mode=climate)
│   ├── tables/varpart/                    # variance_partition (unit/n_units/df_env), climate_confounding (joint_p/n_sites), px_per_variable, dbmem_selected (mode=climate; fitted on one row per SITE since 2026-09-13), dbmem_selection_path
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
│       ├── methods/{method}/              # {method}_pvalues_K{k}.tsv, _sig_snps_{adjust}.tsv, {method}_params_K{k}.tsv,
│       │                                  # {method}_lambda_K{k}.tsv, {method}_thresholds_K{k}.tsv
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
│   ├── tables/gradient_forest/{run_label}_{spatial_tag}/
│   │   └── genetic_offset_{map,site}.tsv
│   └── tables/snp_sets_summary.tsv        # every stored SNP set and whether this run used it
│
├── pipeline_summary.tsv                   # all modes append here
│
├── _work/maf{}_miss{}_smiss{}/ld{}_win{}_step{}/    # parameterized intermediates
│   ├── filtered.vcf, pruned.vcf
│   ├── {name}.geno, .lfmm, .vcfsnp
│   └── {name}.snmfProject/
│
└── _intermediate/                         # internal: samples/, annotation/, {method}/, qs_cache/
    ├── region_params.json                 # per-region Shiny parameter persistence
    ├── run_records/                       # {mode}.json + history.jsonl: what each mode ran with (see Run records)
    └── snp_sets/                          # curated SNP sets read by mode=maladaptation
        ├── manifest.json                  # app-saved sets only (the pipeline never writes it)
        ├── <name>/                        # saved in the GEA tab: selected_snps.tsv, snp_set.yaml
        └── {gea,gwas}.<name>/             # declared in config: selected_snps.tsv, snp_set.yaml, declaration.json

{project_name}_logs/{module}/              # per-module log directories
```

**Path parameterization**:
- `_work/maf0.05_miss0.1_smiss0.5/` — VCF filtered with these thresholds
- `_work/.../ld0.2_win100_step20/` — LD-pruned with R²=0.2, 100kb window, 20kb step
- `Maladaptation/.../my_run_spatial/` — GF run_label + spatial_correction tag

**Method parameters**: the `K{k}` in `{method}_pvalues_K{k}.tsv` and its siblings is `sNMF.k_best`,
not the method's own structure parameter. `{method}_params_K{k}.tsv` beside the table holds the
method's `GEA.configs` / `GWAS.configs` entry as the pipeline resolved it: one row per parameter
(LFMM `K`, EMMAX and GAPIT `n_pcs`, RDA `condition_pcs` and its resolved `predictor_set`, the
significance rule, GWAS `missing_strategy`), each with `set_by` (`config`, `default`,
`sNMF.k_best`, `Climate.predictors` or `registry`).

**Tables beside the plots** (since 2026-10-30). Numbers that used to live only in a figure or a log:

- `{method}_lambda_K{k}.tsv` — the genomic inflation factor of every trait column of the method's
  p-value table (median 1-df chi-square over its null median), with `n_tests` and `correction`:
  `genomic_control` (LFMM with `genomic_control: true`) and `rdadapt_gif` (RDA) mean the method
  already divided its statistics by an inflation factor, so a lambda near 1 is by construction;
  `none` (EMMAX, GAPIT) means the lambda is the scan's own.
- `{method}_thresholds_K{k}.tsv` — per trait column: the rule and value, the threshold it gave,
  `n_tested` (the non-NA p-values: Bonferroni's family, counted per trait column), `n_na_dropped`
  and `n_significant` (the trait's rows in the sig-SNP table). The name carries no rule, so a changed
  adjust rewrites it.
- `PreStructure/tables/tracy_widom.tsv` — LEA's Tracy-Widom test of every PCA axis (eigenvalue,
  statistic, p-value, effective number of markers, the axis's share of the summed eigenvalues as a
  fraction).
- `PreStructure/tables/cross_entropy_K{s}-{e}.tsv` — sNMF's cross-entropy for every K and run;
  `best_run` marks the run each K's Q-matrix comes from and the point the plot draws.
- `climate/tables/present/climate_correlation.tsv`, `Traits/tables/trait_correlation.tsv`,
  `Traits/tables/trait_climate_correlation.tsv` — each correlogram's matrix as pairs (`r`, `n_obs`,
  and `basis`: `sites` when both factors are constant within every site and were correlated over
  one row per site, else `samples`).
- `climate/tables/varpart/dbmem_selected.tsv` and `dbmem_selection_path.tsv` — `selection` says how
  each MEM's flag was decided: `forward_selection`, `fallback_best_marginal_r2` (forward selection
  chose none and the best single axis was kept), `not_selected` or `skipped_insufficient_sites`.
- `Maladaptation/tables/snp_sets_summary.tsv` — every set in the SNP-set store: origin, source
  module, strategy, clumping distance, `n_snps` (counted from its `selected_snps.tsv`),
  `usable_for_offset` (at least 3 SNPs), `used` and `reason` (`used`, `below_min_offset_snps`,
  `not_listed`).

### Run records

Every mode run that has jobs to do writes `_intermediate/run_records/{mode}.json` and appends one
line to `_intermediate/run_records/history.jsonl`. The record names the mode, status, start and
finish (UTC), the command line, each `--configfile` with the sha256 of its bytes, the full config
the run used, the pipeline's git commit with the paths changed since it and a hash of the
uncommitted diff, the Snakemake version, `cpu`, and for `gea`/`gwas` every method's parameters. A
results tree copied elsewhere keeps them. A dry run (`-n`) and a run with nothing to do write
nothing. No author, committer or e-mail address is recorded.

A container cannot see its own image, so the image ID is recorded only when it is passed in:

```bash
docker run ... -e CLINEGO_IMAGE_ID=$(docker images --no-trunc -q cline-go:latest) cline-go:latest snakemake ...
```

Passing it, or any other variable, does not change EMMAX's numbers: `scripts/emmax_run.sh` runs
both EMMAX binaries in a fixed environment (`env -i` with one thread and `MKL_CBWR=COMPATIBLE`).
The bundled `emmax-intel64` links Intel MKL, whose results otherwise follow the size of the process
environment and the thread count, in the fifth significant digit of the p-values.
