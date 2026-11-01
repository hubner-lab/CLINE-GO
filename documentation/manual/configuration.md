# Configuration reference

A project's settings live in one YAML file, `config_<project>.yaml`, in the project folder.
`clinego init <project>` writes it from the pipeline's template (shown in full
[at the end of this page](#the-template)), and the app's settings panel edits the same file.
The keys are grouped, and the spelling, case included, is as shown here.

**Default** below is the value the pipeline uses when a key is absent. Where the template writes
a different value, both are given. Keys marked **required** must be present; a missing one
stops every mode at once. Most values are checked when the run starts, before any job. Settings
used only in benchmarks or rarely changed (`Simulation`, RDA's `fit_mode`, and a few others) are
left out.

## Project

| Key | Default | What it does |
|---|---|---|
| `project_name` | required | Names the results and log folders, `<project_name>_results/` and `<project_name>_logs/`. The app finds a project by its `config_<project_name>.yaml`. Give each dataset its own name. |
| `cpu` | the machine's cores minus 2 | Threads for the steps that use several. Snakemake never gives one job more than `clinego run --cores`. |

## Input

| Key | Default | What it does |
|---|---|---|
| `Input.dir` | required | The folder holding the input files, relative to the project folder, with a trailing slash (`data/`). `""` means the project folder itself. Downloaded climate rasters are kept here too. |
| `Input.vcf` | required | The VCF, `.vcf` or `.vcf.gz`, inside `Input.dir`. |
| `Input.metadata` | required | The metadata TSV: columns `site`, `sample`, `latitude`, `longitude`, then one column per trait (see [Run on your own data](../own-data.md#what-you-need)). |
| `Input.gff` | `""` (none) | A GFF3 gene annotation, optional. Without it, regions are not annotated with genes. |

## Regime

| Key | Default | What it does |
|---|---|---|
| `Regime.mode` | `standard` | `standard`: coordinates and climate, every mode. `gwas_only`: no climate; the climate modes (`climate`, `pregea`, `gea`, `maladaptation`) refuse to run, and `structure` and `gwas` skip their climate maps, `gea_x_gwas` its Miami plot. |

`Climate.enabled` follows from it: `gwas_only` forces it to `false`. Do not set it by hand.

## Filter

| Key | Default | What it does |
|---|---|---|
| `Filter.maf` | required (template 0.05) | Minimum minor allele frequency, 0 to 0.5. |
| `Filter.snp_miss` | required (template 0.1) | Largest share of samples a SNP may be missing in, 0 to 1. |
| `Filter.sample_miss` | 0.5 | Largest share of SNPs a sample may be missing, 0 to 1. |
| `Filter.min_depth`, `Filter.max_depth` | null (off) | Set a genotype to missing when its read depth (`FORMAT/DP`) is below or above this. |
| `Filter.het_outlier_sd` | null (off) | Remove samples whose heterozygosity lies more than this many standard deviations from the mean. |
| `Filter.relatedness` | off (template 0.99) | Allele-sharing (IBS) level above which two samples count as related; 0 to 1. |
| `Filter.relatedness_action` | `keep` | `keep` only marks and counts related pairs; `remove` drops the member of each pair with more missing data. |

## LD

Pruning for linkage disequilibrium, used for PCA, sNMF and the LD-pruned analyses
(plink `--indep-pairwise`).

| Key | Default | What it does |
|---|---|---|
| `LD.window` | required (template 100) | Window size, in SNPs. |
| `LD.step` | required (template 20) | How many SNPs the window moves each step. |
| `LD.r2` | required (template 0.2) | Within a window, one SNP of each pair with r² above this is removed. |

## sNMF

| Key | Default | What it does |
|---|---|---|
| `sNMF.k_start`, `sNMF.k_end` | required (template 2, 6) | The range of *K*, the number of ancestral populations, that `prestructure` tests. |
| `sNMF.repeats` | required (template 10) | Runs per *K*; the best run of each *K* is kept. |
| `sNMF.k_best` | null | The *K* every later mode uses, chosen after `prestructure`; inside the tested range. |
| `sNMF.imputed_maf` | null (= `Filter.maf`) | Minor allele frequency filter applied again after imputation; 0 turns it off. |

sNMF runs diploid only; there is no ploidy key.

## Map

| Key | Default | What it does |
|---|---|---|
| `Map.climate_extent` | `auto` | The area of the climate rasters: `auto` (the sampling sites plus `Map.gap`), `world`, or one string `"min_lon,max_lon,min_lat,max_lat"`. |
| `Map.gap` | 0.5 | Margin around the sites for `auto`, in degrees. |
| `Map.resolution` | 0.5 | WorldClim resolution in arc-minutes: 0.5, 2.5, 5 or 10. At 0.5 the download is very large; the bundled test config uses 2.5. |
| `Map.zoom_extent` | `"NULL"` | An optional zoomed map: `"xmin,xmax,ymin,ymax"`. |

## Climate

| Key | Default | What it does |
|---|---|---|
| `Climate.predictors` | none | The climate variables every GEA method tests, comma-separated: WorldClim's `bio_1` to `bio_19`. Required by the climate modes. |
| `Climate.source` | `worldclim` | `worldclim` downloads the rasters. `custom` reads your own climate tables, named under `Climate.custom` (`present_table`, `future_table` or `future_tables`, `columns`, `grid_resolution`, `key`); their columns must be named `bio_<N>`. |
| `Climate.Varpart.response` | `pcs` | What `climate`'s variance partitioning explains: genomic principal components (`pcs`) or the SNPs themselves (`snps`). |
| `Climate.Varpart.response_var_cutoff` | 0.8 | Keep the leading PCs that explain this share of genetic variance... |
| `Climate.Varpart.response_min_pcs`, `response_max_pcs` | 2, 20 | ...but at least and at most this many. |
| `Climate.Varpart.structure_table` | `qmatrix` | The structure term: the sNMF ancestry at `k_best` (`qmatrix`), or none (`none`). |
| `Climate.Varpart.permutations` | 999 | Permutations for the significance tests. |

## Traits and population statistics

| Key | Default | What it does |
|---|---|---|
| `Traits.pairs_max_factors` | 8 | Above this many traits, `traits` writes a placeholder instead of the pairs grid. |
| `Population.calc_stats` | `false` | Compute Tajima's D, nucleotide diversity and AMOVA in `structure` (needs at least 3 samples per site). |
| `Population.window_size` | 10000 (template 3000000) | Window for Tajima's D and diversity, in bp. |
| `Population.custom_trait_file` | `"NULL"` | A table of `site` plus a value, to size the pie charts by. |

## Piemap

How sites are drawn on the maps.

| Key | Default | What it does |
|---|---|---|
| `Piemap.alpha` | 0.6 | Opacity of the pies, 0 to 1. |
| `Piemap.show_labels` | `false` | Write the site names. |
| `Piemap.label_size` | 10 (template 8) | Size of those names. |
| `Piemap.pie_scale` | 1.0 | Size of the pies. |
| `Piemap.use_points` | `false` | Points instead of pies, for the haplotype maps the app draws. The pipeline's maps always come both ways, and the app switches between them. |

## LDdecay

How `structure` measures linkage-disequilibrium decay. The decay distance sets the distance at
which significant SNPs are merged (see `GEA.snp_clumping_distance`).

| Key | Default | What it does |
|---|---|---|
| `LDdecay.group_by` | `site` | Measure LD within each sampling site (`site`) or each sNMF cluster (`cluster`). |
| `LDdecay.min_samples` | 10 | Groups with fewer samples are skipped. |
| `LDdecay.max_distance` | 500 | Longest SNP distance measured, in kb. |
| `LDdecay.scope` | `both` | `genome_wide`, `per_chromosome` or `both`. `auto_genome_wide` distances need the genome-wide rows, `auto_per_chromosome` the per-chromosome ones. |
| `LDdecay.r2_threshold` | 0.2 | The r² level of the `r2_threshold` rule. |
| `LDdecay.min_signal_ratio` | 2 | A decay curve counts as reliable when its short-range r² is at least this many times the r² of unlinked markers. An `auto_*` distance stops on an unreliable genome-wide curve; then give a distance in bp. |
| `LDdecay.exclude_chromosomes` | none | Contigs left out of the curve, as one string: `"X,Y,MT,Pt"`. |
| `LDdecay.max_pairs`, `background_markers_per_chr`, `background_quantile`, `loess_span`, `seed` | 2e7, 50, 0.95, 0.3, 1 | Sampling and smoothing details; the defaults suit most data. |

## GEA

Genotype-environment association, run by `mode=gea`.

| Key | Default | What it does |
|---|---|---|
| `GEA.configs` | none | The methods to run, one entry each (below). Required by `gea`. |
| `GEA.snp_clumping_distance` | `auto_genome_wide` | One distance for the agreement between methods and for merging significant SNPs into regions: `auto_genome_wide` (the LD decay distance), `auto_per_chromosome` (one per chromosome), or a number of bp. |
| `GEA.clumping_ld_rule` | `critical` | Which decay distance an `auto_*` value reads: where LD meets the r² of unlinked markers (`critical`), where r² halves (`half_decay`), or where it falls to `LDdecay.r2_threshold` (`r2_threshold`). |
| `GEA.region_clumping_distance` | unset | Advanced: a separate distance for the region merge only. |
| `GEA.combine_method` | `All` | How the methods' significant SNPs are joined in the pipeline's own tables (`selected_snps.tsv` and the regions and genes built from it): `All` (every SNP any method calls; also `Union`), `Cross-method` (called by at least two methods within the clumping distance; also `Overlap`), `Cross-method per-trait` (two methods, same predictor; also `MethodOverlap`), or one method's name. The app's GEA tab chooses its own. |
| `GEA.promoter_length` | 10000 | Length of the promoter counted upstream of each gene, in bp. |
| `GEA.wza.window_size`, `GEA.wza.fallback_window_bp` | `auto_genome_wide`, 10000 | Windows for the WZA window-based scores. |
| `GEA.snp_sets` | none | Candidate SNP sets for `maladaptation` (below). Leave it out in a `gwas_only` project, where it stops every mode. |

Each entry of `GEA.configs` names a `method`, a significance rule `adjust` and its
`threshold`, and may set `params`:

```{code-block} yaml
GEA:
  configs:
    - method: "LFMM"
      adjust: "bonf"
      threshold: '0.05'
      params:
        K: 4
    - method: "RDA"
      adjust: "bonf"
      threshold: '0.01'
```

- **Methods:** `LFMM`, `EMMAX`, `RDA`, and the GAPIT models `GLM`, `MLM`, `CMLM`, `ECMLM`,
  `SUPER`, `MLMM`, `FarmCPU`, `BLINK`. Each method appears once.
- **Rules (`adjust`)** and their `threshold`:
  - `bonf`: Bonferroni, a probability in (0, 1]
  - `qval`: q-value (false discovery rate), a probability in (0, 1]
  - `top`: the N SNPs with the lowest p-values, N ≥ 1
  - `toppct`: the lowest-p share of the tested SNPs, a percentage in (0, 100]
  - `custom`: a raw p-value cutoff above 0
- **`params`**, per method (each defaults to the value given):
  - `LFMM`: `K`, the number of latent factors (`sNMF.k_best`); `genomic_control` (`true`)
  - `EMMAX`: `n_pcs`, principal components as covariates (`sNMF.k_best`); `kinship`, `BN` or
    `IBS` (`BN`)
  - `RDA`: `condition_pcs`, principal components conditioned out (`sNMF.k_best`; 0 for no
    structure correction); `axes` (`auto`); `axis_alpha` (0.05); `predictor_set` (`auto`, all of
    `Climate.predictors`, or a comma-separated subset of at least two); `permutations` (999);
    `seed` (42)
  - GAPIT models: `n_pcs` (`sNMF.k_best`), the same for every GAPIT model in the list
- `mode=pregea` recommends values for `K`, `n_pcs` and `condition_pcs`.

### Declared SNP sets

`GEA.snp_sets` lists candidate SNP sets for `maladaptation`. `mode=gea` builds each one into
`_intermediate/snp_sets/gea.<name>/`, the same way the app's "Save SNP set for maladaptation"
does:

```{code-block} yaml
GEA:
  snp_sets:
    - name: default                # each method's own GEA.configs rule, joined
    - name: strict
      threshold_type: bonf         # one rule for every method and predictor
      threshold_value: 0.01
      strategy: Cross-method       # keep SNPs that two methods call
    - name: tuned
      traits: [bio_1, bio_12]      # only these predictors; methods: [...] likewise
      overrides:                   # by method, by predictor, or "predictor::method"
        LFMM: {value: 0.01}
        "bio_12::EMMAX": {type: custom, value: 1.0e-5}
```

- `name`: letters, digits and underscores; the set is stored as `gea.<name>`.
- `threshold_type`, `threshold_value`: one rule (the `adjust` rules above) for every method
  and predictor, instead of each method's own. RDA keeps its own rule.
- `methods`, `traits`: a subset of the methods in `GEA.configs` and of the predictors; all by
  default. Or `pairs`, a list of `"predictor::method"`.
- `overrides`: a rule (`type`, `value` or both) for one method, one predictor, or one
  `"predictor::method"`; the most specific wins.
- `strategy`: `Union` (the default: every SNP any method calls), `Cross-method` (at least two
  methods within `clumping_distance`), `Cross-method per-trait` (two methods on the same
  predictor), `Cross-method >=K` (at least K methods), or one method's name.
- `clumping_distance`: in bp, or `auto_genome_wide`; `GEA.snp_clumping_distance` by default.

## GWAS

Association with the measured traits, run by `mode=gwas`. Its keys work as in GEA.

| Key | Default | What it does |
|---|---|---|
| `GWAS.configs` | none | As `GEA.configs`, with `EMMAX` and the GAPIT models (LFMM and RDA are GEA only). Required by `gwas`. |
| `GWAS.missing_strategy` | `DROP` | Samples missing a trait: `DROP` leaves them out of that trait's test; `MEAN` and `MEDIAN` fill the missing value. |
| `GWAS.traits` | all | A comma-separated subset of the metadata's trait columns. |
| `GWAS.combine_method` | `All` | As in GEA, and not inherited from it. |
| `GWAS.snp_clumping_distance`, `clumping_ld_rule`, `region_clumping_distance`, `promoter_length`, `wza.*` | GEA's | Inherited from `GEA.*` unless set here. |
| `GWAS.snp_sets` | none | As `GEA.snp_sets`, over the traits; stored as `gwas.<name>`. |

## GEAxGWAS

| Key | Default | What it does |
|---|---|---|
| `GEAxGWAS.pairwise.window_size` | 500000 | Two traits' significant SNPs are grouped into windows wherever neighbours lie within this many bp. |
| `GEAxGWAS.pairwise.min_snps` | 2 | A window counts as shared when it holds at least this many SNPs of each trait. |
| `GEAxGWAS.snp_clumping_distance` | from GEA and GWAS | The app's starting distance on the GEA x GWAS tab: GEA's when either side is `auto_*`, else the larger of the two. |

## Future

| Key | Default | What it does |
|---|---|---|
| `Future.ssp` | `'585'` | The CMIP6 scenario: 126, 245, 370 or 585. |
| `Future.year` | `'2061-2080'` | The period: `2021-2040`, `2041-2060`, `2061-2080` or `2081-2100`. |
| `Future.models` | none | Climate models, comma-separated (`'MPI-ESM1-2-HR,UKESM1-0-LL,IPSL-CM6A-LR'`). Required by `maladaptation` with WorldClim. |

## Maladaptation

| Key | Default | What it does |
|---|---|---|
| `Maladaptation.methods` | `gradient_forest` only (template: three) | The offset methods to run, each a key with its settings below: `gradient_forest`, `geometric_offset`, `rda_offset`, `rda_offset_corrected`. A method listed runs. |
| `Maladaptation.snp_sets` | `all` | The SNP sets to fit on: `all` (every stored set with at least 3 SNPs), or a list of names such as `[gea.default]`. |
| `Maladaptation.predictors` | `auto` | The predictors the offset models use: `auto` (all of `Climate.predictors`) or a subset of them. The RDA offset can use at most one fewer than the number of sites. |
| `Maladaptation.emit_plots` | `true` | `false` writes the tables only. |

`gradient_forest`:

| Key | Default | What it does |
|---|---|---|
| `ntree` | 500 (template 1000) | Trees per SNP. |
| `cor_threshold` | 0.5 | Correlation level used for the conditional importance. |
| `spatial_correction` | `with` | Add the dbMEM spatial variables from `climate` as predictors (`with`), not (`without`), or both runs (`both`). `with` needs more sites than `Climate.predictors` plus one. |
| `random_model` | `true` | Also fit a model on random SNPs, as a neutral comparison. |
| `extrap` | `true` | Beyond the training climate: continue each curve linearly (`true`, the gradientForest default) or hold it flat (`false`). |
| `sensitivity_check` | `true` | Refit on site allele frequencies from observed calls only, to check the effect of imputation. |

`geometric_offset`: `scale` (`true`, standardize the predictors) and `k` (blank: `sNMF.k_best`
latent factors).

`rda_offset`: `axes` (`auto`, or an integer of at least 2), `axis_alpha` (0.05),
`permutations` (999), `condition_pcs` (0: no structure correction) and `seed` (42).
`rda_offset_corrected` takes the same keys, inherits `rda_offset`'s, and conditions on 1
principal component by default.

## PreGEA

Settings of the optional `mode=pregea`.

| Key | Default | What it does |
|---|---|---|
| `PreGEA.predictors` | `Climate.predictors` | The predictors the ladders use. |
| `PreGEA.k_offset` | 2 | LFMM's *K* is tried from `sNMF.k_best` minus this to plus this. |
| `PreGEA.n_pcs_max` | 10 | EMMAX's and RDA's principal components are tried from 0 to this. |
| `PreGEA.Advanced.collinearity_r`, `vif_max`, `axis_alpha` | 0.7, 10, 0.05 | RDA's predictor screen and axis test. |
| `PreGEA.Advanced.permutations` | 199 (template 999) | Permutations for RDA's tests. |
| `PreGEA.TransferGuard.enabled` | `false` | Also check the recommendations on the full SNP set (slower). |

## GFF, enrichment and haplotypes

| Key | Default | What it does |
|---|---|---|
| `GFF.feature` | `mRNA` | The GFF feature type counted as a gene. |
| `GFF.go_field` | `"NULL"` | The GFF attribute holding GO terms; `"NULL"` turns GO enrichment off. |
| `GFF.gene_name`, `GFF.biotype` | `description`, `biotype` | Attributes the app's region plot shows. |

The app alone reads the rest, for what it computes on request: `Enrichment.top_terms`,
`plot_width`, `plot_height` and `cnet_label` for the GO enrichment plots, and `haplotype.*` for
the haplotype scan.

## The template

`clinego init` copies this file, `scripts/clinego.app/inst/config_default.yaml` in the
repository, and sets `project_name`. Its values are the template's choices, which differ from
the defaults above in places. Some of its comments predate the launcher:

- `Input.dir` is relative to the project folder, not to `/pipeline/`.
- `LD.window` counts SNPs, not kb.
- `docs/configuration.md`, which a comment names for the SNP-set grammar, is not in the
  repository; [Declared SNP sets](#declared-snp-sets) above has it.
- The `GEAxGWAS` block belongs to `mode=gea_x_gwas` (the comment's `mode=overlapping` is its
  old name), and `threshold_type` also takes `toppct`.

Keys the template leaves out take their defaults: `cpu`, `Climate.source`,
`GEA.combine_method`, `GWAS.combine_method` and `GEA.wza`. A `run_label` under
`gradient_forest` is no longer read (a warning says so): runs are named after their SNP set.

```{literalinclude} ../../scripts/clinego.app/inst/config_default.yaml
:language: yaml
```
