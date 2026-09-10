# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Workflow (MANDATORY)

**Every new feature or change MUST follow this process. No exceptions.**

### 1. Discuss Before Implementing

Before writing any code, **always** start with a discussion:
- Clarify requirements with the user via questions and options
- Identify which existing outputs, scripts, and rules can be reused
- Agree on the approach before touching any files
- Use plan mode for anything non-trivial

### 2. Snakemake First, Scripts Second

Implementation order is strict:

1. **Design the Snakemake rules first** — Define inputs, outputs, and the DAG structure. This is where most effort goes: avoiding redundancy, reusing existing rules/scripts with different inputs, and keeping the workflow clean.
2. **Reuse existing scripts wherever possible** — Many scripts are generic (e.g., `plot_piemap.R` works for ancestry, traits, Tajima's D, genetic offset — just different inputs). Before writing a new script, check if an existing one can handle the job with different arguments.
3. **Add placeholder comments for new scripts** — If a new script is genuinely needed, add the rule with a comment describing what the script should do, its expected inputs/outputs, and argument order. This serves as a specification.
4. **Write scripts last** — Only after the Snakemake structure is finalized and approved, implement the actual R/Python scripts.

### 3. Avoid Redundancy Aggressively

- If two rules differ only in their inputs, use **one generic rule** (or a shared script with different arguments) rather than duplicating logic.
- Check the existing Snakefile for patterns that already solve the problem (e.g., piemap plotting, Manhattan plots, gene finding, enrichment — these are already parameterized).
- When adding a new mode or analysis, map its steps onto existing scripts before proposing new ones.

## Project Overview

CLINE-GO is a Dockerized, Snakemake-based bioinformatics pipeline for population genomics performing VCF preprocessing, population structure analysis (PCA, sNMF), GWAS/GEA (EMMAX, LFMM), and maladaptation assessment (Gradient Forest).

## Pipeline Philosophy

**Core Principle: Reduce Data Uncertainty**

Biological data is messy. Different tools and file formats introduce inconsistencies (chr1 vs 1 vs 1H, silent format changes, mismatched identifiers).

**CLINE-GO's approach**:
1. **Normalize early, normalize consistently** - Standardize at the earliest step
2. **Enforce consistency across all outputs** - All outputs maintain same standards
3. **Fail loudly, not silently** - Error early rather than propagate bad data
4. **Document all transformations** - Log normalization steps

**Example**: LEA's `vcf2lfmm()` strips "chr" prefix, so we normalize BOTH VCF and GFF early (in processing mode) ensuring all downstream outputs have consistent chromosome names.

**Code Structure Rules**:
1. **Compute once, reuse downstream** — Each stage produces finished outputs. Downstream scripts use them as-is, never re-derive or re-extend.
2. **No pass-through parameters** — Don't pass values a script doesn't use. If a parameter only matters at creation time, only the creation script should receive it.
3. **No dead code paths** — Don't add conditional logic for cases that never occur (e.g., `if distance > 0` when always called with `0`). Remove it instead.
4. **Single source of truth** — Region boundaries live in the regions table. Gene distances live in the genes table. Don't recompute from raw parameters downstream.
5. **One responsibility per script** — Read upstream outputs, do one thing, write downstream outputs. Don't mix concerns.
6. **Design for WGS density** — The pipeline processes whole-genome sequencing data with hundreds of thousands to millions of SNPs. Never use raw data points (geom_point, scatter) in plots where the number of observations scales with SNP count or pairwise comparisons. Use binned summaries, smoothed curves (geom_smooth/loess), or density representations instead. Individual points become unreadable noise at WGS scale.
7. **Prefer the most intuitive data transformation** — When multiple approaches produce equivalent output, choose the one a user can follow step-by-step through the analysis without needing to know implementation details. For example: filtering data before plotting is clearer than filtering inside a plot layer; explicit `df_background <- plot_data_all` (all SNPs) is clearer than a hidden `filter(method == first)` that silently drops methods. When a single transformation rule applies across multiple scripts or contexts (per-method, combined, Miami), use the **same rule everywhere** — inconsistency between scripts with no documented reason creates unpredictable output and maintenance debt.

   **Example — Manhattan background**: Background PNG = all SNPs, all methods (no filtering). Interactive overlay = sig SNP highlight markers. Toggle = highlight on/off. Every sig SNP always has a ghost dot behind it. Same rule in `plot_manhattan.R`, `plot_manhattan_combined.R`, `plot_miami.R`. Never use a "first-method-only" de-dup shortcut in the background — it makes behavior differ unpredictably by method order.

8. **Instructional text belongs in Shiny, not baked into plots** — Plots (PNG/SVG) should carry only structural elements: title, axis labels, legend, geometry that directly encodes data (bar/point value labels, `geom_vline`/`geom_hline` reference lines). Never add floating `annotate("text", ...)`/`labs(subtitle=...)` commentary that *explains* a count, threshold, or what a dashed line means (e.g. "0 SNPs below threshold", "Dashed lines: missingness threshold + het mean ± SD") — even when the number itself is data-derived. Those overflow small plot cards, can't be styled/positioned responsively, and duplicate what the app already says. Put that explanation in the `note` reactive passed to `mod_image_card_server()` instead (rendered as the card's hover badge/tooltip via `filter_note()`/`relatedness_note()`/etc. in `utils_ui.R`) — same content, one place, actually readable. If the number needed for the note isn't already computed in Shiny, derive it from a table the app already loads (e.g. `filtering_summary.tsv` via `load_filtering_summary()`) rather than writing a new file just to carry one count.

   **Exception**: an empty-state placeholder that IS the entire plot content when there's no data to show (e.g. "Not enough samples for relatedness MDS") is fine to keep on-plot — it's not commentary layered on top of real data, matching the `plot_placeholder()` convention used elsewhere in the app.

9. **Use the shared theme for new plots** — `scripts/R/utils/theme_clinego.R` (`source()` it, same pattern as `scripts/R/utils/manhattan_utils.R`) defines the pipeline's plot look, rolled out to `plot_qc_processing.R` first and intended to extend to the rest of `scripts/*.R` over time. Reuse it rather than inventing per-script themes/palettes:
   - `theme_clinego()` — "Publication Classic": `theme_classic()` base (axis lines, no gridlines), no `base_family` (Docker ships no fonts; setting one silently falls back to default sans while implying a match that isn't there).
   - Semantic colors (from `wesanderson::wes_palette("Rushmore1")`'s last 3 entries, user preference) — always use these constants, never raw hex, so meaning stays consistent everywhere: `CLINEGO_RETAINED` (green `#0B775E`, also aliased as `CLINEGO_NEUTRAL` for plain/no-flag data points — no grey dots by request), `CLINEGO_REMOVED` (red `#F2300F`, flagged/dropped/discarded), `CLINEGO_THRESHOLD` (dark plum `#35274A`, cutoff/reference lines).
   - `CLINEGO_CATEGORICAL` — Okabe-Ito colorblind-safe palette + `scale_color_clinego()`/`scale_fill_clinego()`, for multi-level categorical grouping (traits, methods, etc.) — converges the palette already duplicated across `manhattan_utils.R`/`fct_manhattan.R`/`mod_gea.R`.

## Build and Run Commands

### Build Docker Image
```bash
docker build -t cline-go .
```
Rebuild required after any Dockerfile or R package version change.

### Run Pipeline (SIMDATA — main testing dataset)
```bash
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline cline-go:latest \
  snakemake -c4 -s Snakefile --config mode=<MODE> --configfile config_SIMDATA.yaml --scheduler greedy
```

### Run Pipeline (additional real-data project — validation only, on demand)
No default config file for this anymore (see Test Datasets below) — create a
`config_<PROJECT>.yaml` naming the project after its actual VCF before running:
```bash
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline cline-go:latest \
  snakemake -c4 -s Snakefile --config mode=<MODE> --configfile config_<PROJECT>.yaml --scheduler greedy
```

### Interactive Docker Entry
```bash
docker run -w /pipeline --user $(id -u):$(id -g) -it -v $PWD:/pipeline cline-go bash
```

### Dry Run (check what would execute)
```bash
docker run --user $(id -u):$(id -g) --rm -v $PWD:/pipeline cline-go:latest \
  snakemake -n -s Snakefile --config mode=<MODE> --configfile config_SIMDATA.yaml --scheduler greedy
```

**Pipeline modes**: `processing`, `prestructure`, `structure`, `climate`, `traits`, `pregea`, `gea`, `gwas`, `gea_x_gwas`, `maladaptation`

*`climate` — predictor characterization: correlation heatmap, density plots, invariant-predictor detection, plus the shared spatial artifacts — dbMEM eigenvectors + climate/structure/geography variance partitioning (`climate/{plots,tables}/{spatial,varpart}/`). The dbMEM/varpart products are reused by the spatial Gradient Forest, so **run `mode=climate` before requesting a spatial (or `both`) Gradient Forest in `mode=maladaptation`** — there is no config-time guard for this ordering anymore (varpart is now unconditional inside `mode=climate`, no `enabled` switch). WorldClim download / custom-climate staging stays in `mode=structure` (climate_site feeds Structure's own piemaps); the climate/structure boundary is "predictor characterization" vs "raw data acquisition".*

*`traits` — the phenotypic counterpart of `mode=climate`: per-trait density plots, a traits-only correlogram, a pairs grid, and trait summary / invariant-trait tables (`Traits/{plots,tables}/`). **Runs in both regimes and needs no coordinates and no climate** — that is the point: `mode=climate` raises when `Climate.enabled: false`, so in `gwas_only` (where traits are the only factors a project has) trait characterization used to be unreachable. The one climate-dependent product is the joint traits x climate correlogram (`correlation_heatmap_traits_climate.png`), emitted only when `Climate.enabled: true`, which means `mode=structure` (climate extraction) must have run first; the Shiny Phenotypic tab toggles between the two correlograms and hides the toggle when the joint one is absent. `Traits.pairs_max_factors` (default 8) caps the pairs grid — above it the rule writes a placeholder PNG+SVG instead of an unreadable k x k grid. Phenotype density moved here from `mode=climate`; the climate correlogram is climate-only now.*

*`pregea` (optional, added alongside RDA integration) — LD-pruned-only hyperparameter exploration: LFMM-K / EMMAX-#PC / RDA Condition()-PC ladders, now ONE decision (LFMM/EMMAX/RDA always run together, no per-block switches). EMMAX #PCs and RDA Condition()-PCs sweep ONE shared range (`PreGEA.n_pcs_max`). Predictor characterization + dbMEM/varpart moved out to `mode=climate`. Writes `PreGEA/tables/pregea_recommendations.tsv`, one row per (method, param); the Shiny GEA tab's method editor reads it for pre-fill/Apply badges. See `docs/rda_research.md` Part C for the design rationale.*

*`haplotype_scan` / `haplotype` modes removed in Phase 3 — haplotype analysis now runs interactively in the Shiny app (Region Explorer → Run Haplotype Scan / Run Haplotype Viz).*

**CRITICAL: Never run multiple modes in parallel for the same project.** All modes share the same `{PROJECT}_results/` working directory and Snakemake lock. Running two modes concurrently (e.g., `association` and `association_phenotypes` for SIMDATA) causes lock conflicts and potential file corruption. Always run modes sequentially. Parallel runs are only safe across **different projects** (e.g., two projects with distinct `project_name` values, each with its own `_results/` directory).

**CRITICAL: `{PROJECT}_results/` is keyed by `project_name` alone — the pipeline never checks that a config's `Input.vcf/metadata/gff` still matches what actually produced the existing outputs in that folder.** Two different config files sharing the same `project_name` (e.g., the CLI's config file and the Shiny app's own `config_{project_name}.yaml`, written separately from the Home tab's Input Files section) will both read/write the same `_results/` directory even if they point at completely different raw data. Non-parameterized paths (`Processing/tables/metadata.tsv`, `PreStructure/plots/K*/`, etc.) get silently overwritten by whichever config ran most recently — this already happened once (2026-07-21, a CLI run and a Shiny-triggered run both targeting `project_name: TEST` with different `Input.vcf`, three minutes apart). **Name projects after their actual VCF/dataset, never reuse a project name across datasets, and never run the CLI and the Shiny "Run" button against the same project concurrently.**

**Snakemake debug flags**: `-n` (dry run), `-R <rule>` (rerun rule + downstream), `--forcerun <rule>` (rerun only rule), `-F` (force all), `-p` (print shell commands)

## Architecture

### Key Files
- `Snakefile` - Thin orchestrator (~24 lines): includes + `rule all`
- `workflow/rules/common.smk` - Config parsing, path dicts, helpers, `get_targets()`
- `workflow/rules/{module}.smk` - Per-module rules: `processing`, `structure`, `structure_k`, `association`, `phenotype_assoc`, `overlapping`, `maladaptation`, `summary`
- `workflow/rules/regionplot.smk` - Only `gff2topr` rule remains (used by Shiny on-demand regionplots); `mode=regionplot` is deprecated
- `Snakefile_old` - Original workflow (**do NOT modify**)
- `config_SIMDATA.yaml` - Main testing dataset config (nested YAML). Additional per-project configs (`config_<PROJECT>.yaml`, named after the actual VCF/dataset) are created on demand for real-data validation runs — see Test Datasets.
- `scripts/*.R` - R scripts (executed at `/pipeline/scripts/` inside Docker)
- `scripts/emmax-intel64`, `scripts/emmax-kin-intel64` - Pre-built EMMAX binaries (checked in, not built from source)
- `scripts/app.R` - Shiny interactive results viewer (~2,500 lines) -- **NOT YET UPDATED** for new paths
- `Dockerfile` - Container with pinned package versions (rocker/shiny:4.5, Bioconductor 3.22)

### Configuration Groups (Nested YAML)

Config uses nested YAML groups. Old flat `UPPER_SNAKE_CASE` keys are auto-migrated via `_migrate_config()` with deprecation warnings.

1. **input** - `dir`, `vcf`, `metadata`, `gff` + top-level `project_name`, `cpu`
1b. **Regime** - `mode`: `standard` (geography + climate: structure, GEA, maladaptation) or `gwas_only` (no coordinates: Home/Processing/PreStructure/Structure/GWAS). Declared once at project creation and **hard-validated** in `common.smk` — an unknown value raises at parse time. `Climate.enabled` is **derived** from it (`gwas_only` forces `false`) on every write the Shiny app makes, and forced to `False` in `common.smk` for a hand-written CLI config; it is no longer user-editable in the app. A project predating the key is read as `gwas_only` only when it sets `Climate.enabled: false` explicitly, else `standard`.
2. **filter** - `maf`, `snp_miss`, `sample_miss`
3. **ld** - `window`, `step`, `r2`
4. **snmf** - `k_start`, `k_end`, `k_best`, `ploidy`, `repeats`
5. **map** - `climate_extent` (was MAP_CROP), `gap`, `resolution`, `zoom_extent` (was MAP_REGIONMAP_EXTENT)
6. **climate** - `predictors` (was CLIMATE_VARS)
7. **pop** - `calc_stats`, `window_size`, `custom_trait_file`
8. **piemap** - `alpha`, `show_labels`, `label_size`, `pie_scale`, `use_points`
9. **association** - `configs` (list of method/adjust/threshold), `combine_method`, `combine_gap`, `sig_snp_distance`, `region_distance`, `top_regions`, `promoter_length`, `go_field`, `scattermore_threshold`
10. **gff** - `feature`, `gene_name`, `biotype`
11. **enrichment** - `top_terms`, `plot_width`, `plot_height`, `cnet_label`, `top_plot_regions`
12. **future** - `ssp`, `year`, `models`
14. **gradient_forest** - `ntree`, `cor_threshold`, `spatial_correction` (was GF_PCNM), `run_label` (was `suffix`), `random_model`
15. **phenotype_association** - `configs`, `missing_strategy`; inherits from `association.*` by default, override only if different
16. **overlap** - `region_distance` (Shiny filter bar default, falls back to `max(association, phenotype_association)`); `pairwise.window_size`, `pairwise.min_snps` for pairwise trait overlap table. Overlap regions/genes/enrichment are computed interactively in Shiny — not pipeline-side.
17. **haplotype** - `scan.regions_source`, `scan.regions_file`, `scan.top_regions`, `scan.min_snps`, `scan.min_group_size` (was MGMIN), `scan.min_haplotype_size` (was MINHAP), `scan.epsilon_range`, `scan.metadata_type`, `epsilon_selected`

### Output Organization

Organized by **module** (matching pipeline modes). Each module owns its plots and tables.

```
{PROJECT}_results/
├── Processing/tables/                     # metadata.tsv, sample_missing_stats.tsv
├── PreStructure/
│   ├── plots/                             # pca.png, tracy_widom.png, cross_entropy_K{start}-{end}.png
│   │   └── K{k}/                          # structure_K{k}.png, pca_structure_K{k}.png, pop_diff_K{k}.png
│   └── tables/
│       └── K{k}/                          # clusters_K{k}.tsv (Q-matrices)
├── climate/                               # mode=climate owns predictor characterization + spatial/varpart
│   ├── plots/                             # density_plot_present, correlation_heatmap, density_plot_future_*
│   ├── plots/{spatial,varpart}/           # dbmem_screeplot, varpart_venn (nested donut), px_barplot, dbmem_selection_path
│   ├── tables/present/                    # climate_present_all.tsv, _site.tsv, _site_scaled.tsv, climate_invariant_predictors.tsv
│   ├── tables/spatial/                    # dbmem_vectors.tsv, dbmem_diagnostics.tsv
│   ├── tables/varpart/                    # variance_partition, climate_confounding, px_per_variable, dbmem_selected, dbmem_selection_path
│   ├── tables/future/                     # climate_future_year{Y}_ssp{S}_site.tsv, _all.tsv
│   └── rasters/{present,future}/          # WorldClim .tif rasters (terra)
├── Traits/                                # mode=traits — phenotypic factor characterization, BOTH regimes
│   ├── plots/                             # density_plot_phenotypes, correlation_heatmap_traits,
│   │                                      # correlation_heatmap_traits_climate (standard only), trait_pairs
│   └── tables/                            # trait_summary.tsv, trait_invariant.tsv
├── Structure/
│   ├── plots/piemap/                      # piemap_{bio}.png/svg/qs + _points.png/svg/qs (clear-map companion) + zoom/
│   ├── plots/piemap/{tajima_d,pi_diversity}/  # trait-scaled piemaps (optional)
│   ├── plots/pop_stats/                   # mantel_test, amova (optional)
│   └── tables/pop_stats/                  # tajima_d_by_pop, pi_diversity_by_pop, ibd_*, amova
├── PreGEA/                                # optional, mode=pregea — grid-level plots only, no per-rung files (EXCEPT RDA per-model artifacts below)
│   ├── plots/{structure,lfmm,emmax,rda,transfer}/
│   ├── plots/rda/models/pc{n}/            # per-Condition()-PC: biplot.png/svg, axis_screeplot.png/svg (Shiny RDA tab selector)
│   ├── tables/{structure,lfmm,emmax,rda}/
│   ├── tables/rda/models/pc{n}/           # per-Condition()-PC: axis_anova.tsv
│   ├── tables/pregea_recommendations.tsv  # one row per (method, param); read by GEA tab's method editor
│   └── tables/pregea_transfer_guard.tsv   # opt-in (PreGEA.TransferGuard.enabled)
├── GEA/
│   ├── GAPIT_native_output/{model}/       # raw GAPIT output files
│   ├── plots/manhattan/
│   │   ├── {method}/                      # manhattan_{trait}_K{k}_{adjust}.png/svg, qq_{trait}_K{k}_{adjust}.png/svg
│   │   └── combined/                      # manhattan_combined_K{k}.png/svg, qq_combined_K{k}.png/svg
│   ├── plots/enrichment/{trait}/          # region_{id}_dotplot/emapplot/cnetplot
│   └── tables/
│       ├── methods/{method}/              # {method}_pvalues_K{k}.tsv, _sig_snps_{adjust}.tsv
│       ├── selected_snps.tsv, regions_per_trait.tsv, regions_combined.tsv
│       ├── genes_per_region.tsv, genes_per_region_collapsed.tsv, genes_combined.tsv
│       └── enrichment/{trait}/            # GO enrichment TSVs per region
├── GWAS/                                  # same structure as GEA/ + phenomap piemaps
├── GEAxGWAS/
│   ├── plots/miami_combined_K{k}.{png,svg,_background.png,_coords.json}
│   └── tables/pairwise_{overlap_table,collapsed_snps}.tsv
│   # NOTE: overlap regions/genes/enrichment are fully interactive in Shiny (not pipeline-computed)
├── Maladaptation/
│   ├── plots/{method}/{SUFFIX}/           # cumulative_importance, overall_importance, genetic_offset_piemap[_{tajima_d,pi_diversity,points}]
│   │   └── zoom/{coords}/                # zoomed piemaps
│   └── tables/{method}/{SUFFIX}/          # genetic_offset_map, genetic_offset_site
├── haplotype_scan/{tag}/                  # clustree plots, selected_regions.tsv, scan_status.tsv
├── haplotype/{tag}/                       # crosshap_viz, boxplots, haplotype piemaps, assignment/frequency tables
├── pipeline_summary.tsv                   # all modes append here
├── _work/maf{}_miss{}_smiss{}/ld{}_win{}_step{}/  # parameterized intermediates
└── _intermediate/                         # internal: samples/, annotation/, enrichment/, {method}/{SUFFIX}/, haplotype/, qs_cache/, flags/

{PROJECT}_logs/{module}/                   # per-module log directories
```

### Workflow Dependencies

**Pipeline flow**: Processing → PreStructure → Structure → Climate / Traits → PreGEA (optional) → GEA/GWAS → GEAxGWAS → Maladaptation

(Climate, Traits and PreGEA are all siblings — none is a hard input to another; the spatial Gradient Forest in Maladaptation depends on Climate's dbMEM/varpart outputs. Traits needs only Processing's `metadata.tsv`, except for its optional joint traits x climate correlogram, which needs Structure's climate extraction.)

Each mode is run separately via `--config mode=<MODE>`.

## Test Datasets (CRITICAL)

**SIMDATA is the main, primary testing dataset (`config_SIMDATA.yaml` / `SIMDATA_results/`).** All routine development and testing runs use SIMDATA — fast, self-contained, no ambiguity about what it contains.

Real-data projects (e.g. a specific WGS/GBS dataset) are **additional, on-demand validation** — not the default. There is no standing "TEST" project or config file anymore (2026-07-21: a leftover `config.yaml`/`TEST_results/` pairing under the generic name `TEST` collided with an unrelated Shiny-side project also named `TEST`, silently mixing two different datasets' outputs — see the CRITICAL note under Build and Run Commands). When a real dataset is needed:
- Create `config_<PROJECT>.yaml` naming the project after its actual VCF/dataset (e.g. `config_WBDC.yaml`, `project_name: WBDC`) — never the generic name `TEST`.
- This applies to both the CLI config and the Shiny app's own per-project `config_{project_name}.yaml` (written from the Home tab's Input Files section) — **keep the same project name and Input files in both, or don't reuse the name at all.**

| | **SIMDATA** |
|---|---|
| **Config** | `config_SIMDATA.yaml` |
| **Size** | 51 samples / 9 sites raw (50 samples after filtering), 350 SNPs raw / 337 filtered / 230 LD-pruned, 5 chr, 3 traits (height, flowering_time, disease_score) |
| **Speed** | Seconds to ~2 min |
| **Purpose** | Primary testing dataset — all routine development |
| **Features** | 3 original pops (Negev/TelAviv/Galilee) + 6 preGEA sites added by `scripts/add_pregea_sites.R` (spatially IDW-interpolated genotypes/phenotypes, inside the original coordinate bounding box so `Climate.climate_extent: auto` reuses the cached WorldClim raster), missing data test, climate-associated SNPs, GO terms, relatedness-test duplicate samples (`*_DUP`, `scripts/add_related_samples.R`) |

**`test_data/config_testdata.yaml` is a separate, tracked copy of SIMDATA for GitHub** (2026-08-19) — `data/` and `config_SIMDATA.yaml` are gitignored, so a fresh clone had no runnable dataset. `test_data/` ships copies of the SIMDATA VCF/metadata/GFF renamed to `testdata.*` (`project_name: testdata`), plus `config_testdata.yaml` (same as `config_SIMDATA.yaml` except `Map.resolution: 2.5` — avoids a 10.4 GB WorldClim download on first clone). **Do not confuse the two**: `config_SIMDATA.yaml` is the local working config (37 GB WorldClim cache in `data/`, 30s resolution) — never edit `test_data/*` expecting it to affect local SIMDATA runs, and never point local dev at `test_data/`. See README "Test Dataset" section.
**Rebuilding SIMDATA from scratch** (it is gitignored — `data/` and `config_*.yaml` both are — so a fresh clone has no fixture; rebuilt 2026-08-15 after exactly that):
1. `Rscript scripts/generate_simdata.R data/` — writes **only** `data/SIMDATA.vcf` and `data/SIMDATA.gff3`
2. **Write `data/SIMDATA_metadata.tsv` by hand** — the generator does not produce it. Columns must be exactly `site sample latitude longitude height flowering_time disease_score` (both injector scripts below hardcode those trait names), sample names `NEG01-10`/`TAV01-10`/`GAL01-10`, and coordinates **identical within a site** (Negev 30.854/34.7826, TelAviv 32.0837/34.7817, Galilee 33.0128/35.4985) — `add_pregea_sites.R:74` asserts exactly 3 distinct `(site, lat, lon)` triples and dies otherwise
3. `Rscript scripts/add_related_samples.R` **then** `Rscript scripts/add_pregea_sites.R` — that order matters: the second excludes `_DUP` samples from its IDW anchor set
4. Write `config_SIMDATA.yaml` from `scripts/clinego.app/inst/config_default.yaml`; `sNMF.k_best: 3` is ground truth (3 simulated ancestral populations), not a cross-entropy guess
5. WorldClim: `data/wc2.1_30s/` is a cached 11 GB global extract — as long as it is present, `mode=structure` does no download

**Testing workflow**:
1. Make code changes to `Snakefile` or `scripts/*.R`
2. Test with SIMDATA (`config_SIMDATA.yaml`)
3. Only spin up a real-data project (its own `config_<PROJECT>.yaml`) for validation the user explicitly requests

### Retired: Láruson et al. 2022 — unusable as a GEA benchmark (2026-08-01)

Archived to `_archive/laruson/`. **Do not revive it for GEA benchmarking.** Two disqualifying properties, both confirmed from the Dryad archive's own simulation parameters and from our runs:

1. **No population structure.** SLiM parameters are `m = 0.2`, `n = 100` → Nm = 20 migrants/deme/generation, expected Fst ≈ 0.012. Measured: PC1 = 0.38% of variance, PC1–3 = 0.94%. Worse, **PC1/PC2 are just a rotation of the two environmental axes** (R²_env = 0.835 / 0.838; PC3 onward ≈ 0.006). There is no structure confounding the environment — the only structure present *is* the adaptive signal. Consequence: structure correction can only destroy signal (EMMAX `n_pcs` 0 → 3 collapses AUC-PR 0.580 → 0.052; RDA `condition_pcs` likewise), so the dataset cannot exercise or validate the pipeline's structure-correction machinery at all. It was built to test Gradient Forest genetic-offset prediction, not GEA.
2. **Truth set incompatible with our MAF filter.** `causal_mutations_pos_filtered.txt` lists causal variants already MAF-filtered by the authors at **MAF ≥ 0.01** (minimum causal MAF is exactly 0.010000), but the shipped VCF is their *unfiltered* one. At our `Filter.maf: 0.05` only **28 of 102** causal loci are testable; at 0.01 all 102 are. Any recall figure against the 102 denominator is capped at ~27% by construction.

Also: LD decays to r²=0.2 at **591 bp** (`r = 1e-5`, ≈5 Morgans over the 500 kb contig), so causal loci cannot be picked up via linked markers — and WZA is useless here (only 50 windows genome-wide).

The reusable parts were kept out of the archive: `benchmarks/lib_detection.R`, `eval_detection.R`, `sweep_thresholds.R`, `sweep_rda.sh`, `score_ladders.sh` are dataset-agnostic (they take any wide p-value table + a `chr/pos/category` truth table) and should be reused by the next simulation benchmark.

### Testing Guidelines

**CRITICAL - Preserve output directories**:
- **DO NOT** remove `SIMDATA_results/`, `SIMDATA_logs/` (or any other project's `_results/`/`_logs/`) between runs
- Snakemake skips rules with valid outputs (saves time)
- `download_climate_present` is very slow - preserve it!

**Force re-runs** (when needed):
- `snakemake -R <rule>` - Re-run rule + downstream
- `snakemake --forcerun <rule>` - Re-run only this rule
- `snakemake -F` - Force all rules

Prefer Snakemake flags over manually removing files.

### What cannot be tested on mac-studio (and the rule for it)

mac-studio is arm64 (colima, Virtualization.framework, 8 CPU / 16 GiB). The image is
**amd64-only** — `rocker/shiny:4.5` publishes no arm64 manifest — so every `docker run`
here is QEMU emulation, and some things simply do not work. **Verified 2026-08-21:**

| What | Symptom on mac-studio | Where it does work |
|---|---|---|
| `emmax-kin-intel64` / `emmax-intel64` | Used to **hang forever** after `Identified N individuals` / `nex = 0` — no crash, no exit, `docker run` sat there until killed. Every call site now goes through `scripts/emmax_run.sh`, which refuses to launch on an emulated x86-64 host and **fails in under a second** with the reason. Override for experiments: `-e ADAPTOGENE_ALLOW_EMULATED_EMMAX=1` (pins Intel OpenMP/MKL to one thread; may still hang). | Linux (x86) |
| `mode=pregea` (end to end) | Fails at `pregea_kinship_pruned` (the EMMAX preflight, immediately — no longer a hang), which is upstream of `pregea_emmax_ladder` → `pregea_recommend`. The LFMM ladder, the RDA setup and the screeplot have no EMMAX dependency, so sub-targets that avoid that branch DO run — e.g. `snakemake … <file>` for `PreGEA/tables/rda/rda_predictor_collinearity.tsv`. **Still never verified end to end on any machine.** | Linux |
| `mode=gwas`, `mode=gea` with EMMAX in `association.configs` | Same preflight failure (was: same hang). GAPIT/LFMM/RDA methods are unaffected. | Linux |
| `docker build` | Random `gcc: internal compiler error: Segmentation fault … cc1` under QEMU, a different package each time. Retry — layers that compiled are cached, so each attempt resumes. ~11+ min. | native x86 |
| `-c4` on a heavy mode | OOM-kills the container (exit 137) inside the 16 GiB VM. Use `--workflow-profile workflow/profiles/lowmem` (cores 2 + a `mem_mb` budget that bounds concurrency; no rule declares `resources:` otherwise, so `-cN` is the only dial and it ignores memory entirely). A kill also leaves a Snakemake lock → `snakemake --unlock` (see Known Quirks). | bigger box |

**Always pass `--name` to `docker run` on this machine.** `--rm` removes a container that *exits*; killing the docker **client** (Ctrl-C, a timeout) leaves the **container** running, still holding the `{PROJECT}_results/` lock. Three orphans once fought over one results directory before anyone noticed. With `--name adaptogene-run` the orphan is one `docker stop adaptogene-run` away; without it, `docker ps` and guesswork.


Per-architecture EMMAX binaries do **not** fix this: architecture is a property of the
image, so inside an amd64 container `uname -m` is always `x86_64` and the arm64 binary
could not execute anyway. The only real fix is a multi-arch image (base would have to
move to `rocker/r-ver:4.5`, which does publish arm64, and the Shiny server layer
re-added by hand) — a project, not a task.

**THE RULE — when you cannot test something here, write it down, in both places:**

1. Add a row to the table above if it is a *new* class of mac-studio limitation.
2. Append a line to the ADAPTOGENE dossier's `## Findings` (`~/Orthidian/AGENTS.md` §13
   grammar) naming the rule/script that went unverified and why.

Never let "could not run it" quietly become "ran fine". A change that shipped untested
must say so where the next session will see it — a green review with an unverified
script in the diff is the exact failure this rule exists to prevent.

## Important Rules

1. **Do NOT read/view image files** (PNG, SVG, JPG) - User checks plots themselves
2. **Do NOT modify `Snakefile_old`** - Reference only
3. **Chromosome normalization** - Pipeline strips "chr" prefix early (chr1→1, chr2H→2H). Both VCF and GFF normalized in processing mode. All outputs use normalized names.
4. **Log every manual workaround to `docs/pipeline_improvement_requests.md`** — see below. This is mandatory in every session, not optional.

## MANDATORY: Log pipeline gaps as you hit them

**Whenever you have to work *around* the pipeline, append an entry to
[`docs/pipeline_improvement_requests.md`](docs/pipeline_improvement_requests.md) — in the same
session, at the moment you hit it, not at the end.**

The trigger is any of:

- **You wrote a helper script to obtain something the pipeline already computed but never wrote to
  disk.** (Canonical case: cross-entropy is rendered as a PNG only, and rule 1 forbids reading
  plots — so `sNMF.k_best` could not be chosen from pipeline output at all without extracting the
  values from the `.snmfProject` by hand.)
- **You hand-built an output the pipeline cannot emit**, even though the logic exists somewhere in
  the codebase (e.g. `combine_sigsnps()` implements the ≥2-method consensus rule, but
  `SIGSNPS_METHOD` is hardcoded so batch mode can only union).
- **A config value was accepted at parse time and rejected at run time**, after other rules had
  already burned compute.
- **An output is misleading, stale, or contradicts a sibling output** in the same directory.
- **You needed a number to make a decision and no table contained it** — especially a dispersion
  or uncertainty measure that turns "the minimum is at X" into "the minimum is inside the noise".
- **Cross-run / cross-project comparison** — nothing in the pipeline compares two `_results/`
  dirs, so every multi-arm design re-invents this.

**Entry format:** what happened → why it cost something → what would fix it. Cite `file:line` when
known. Tag `[workaround]` and name the script when a script was written — those are the strongest
candidates, since the fix already exists and only needs adopting.

**This file is a request queue, not a changelog.** Under `CLAUDE.local.md`'s USE-ONLY mode you log
the request and move on; you do not implement it. The user decides later what graduates into the
pipeline. Do not let a gap live only in the session transcript — that is exactly what this file
exists to prevent.

## Common Development Tasks

### Adding a New Plot
1. Create `scripts/plot_<name>.R`
2. Add rule to appropriate `workflow/rules/{module}.smk`
3. Define output path in `O` dictionary (in `common.smk`)
4. Add to mode's target list in `get_targets()`
5. Test with SIMDATA, then TEST

### Adding a GEA/GWAS Association Method
1. Create `scripts/<method>.R`
2. Add entry to `workflow/methods/gea.py` (`GEA_METHODS` dict) with `engine`, `script`, and capability flags
3. Add config parameters if needed (in `common.smk` config parsing)
4. Test with SIMDATA

### Adding a Maladaptation Method
1. Add entry to `workflow/methods/maladaptation.py` (`MALADAPTATION_METHODS` dict) — fields: `engine`, `model_script`, `offset_script`, `cumimp_script`, `importance_script`, `supports_spatial`, `supports_random_model`
2. Add a new rule block in `workflow/rules/maladaptation.smk` using the `mala_*` template functions (e.g., `mala_model(method, ...)`, `mala_offset_piemap(method, ...)`) — output paths automatically land under `Maladaptation/{plots,tables}/{method}/{run_label}_{spatial_tag}/`
3. Add config block under `Maladaptation.methods.<method_name>:` in the relevant config YAML files
4. Test with SIMDATA

### Debugging Failed Rule
1. Check `{PROJECT}_logs/{module}/<rule_name>.log`
2. Verify input files exist and aren't empty
3. Re-run with `-p` flag for verbose output
4. Test R scripts interactively inside Docker

## Technical Details

### File Format Conversions
- **VCF → GENO/LFMM**: `vcf2lfmm.R` (LEA strips "chr" prefix - why we normalize early)
- **LFMM → VCF**: `lfmm2vcf.R`
- **VCF → TPED/TFAM**: `tped_assoc` rule (plink, separate from emmax.R)
- **VCF → GD/GM**: `vcf_to_gapit_numeric.R` (numeric 0/1/2 for GAPIT)
- **GFF3 → topr**: `gff2topr.py`

### Key R Packages
- **LEA** - sNMF, PCA, imputation
- **vcfR** - VCF manipulation
- **terra** (1.8-5) - Raster operations (replaced `raster` package — 64-bit, no overflow at 30s resolution)
- **gradientForest** - Landscape genomics
- **qvalue** - FDR correction
- **topr** (≥2.0.0) - Regional Manhattan plots
- **geodata** - WorldClim download
- **scattermore** - Fast Manhattan rendering (>30k SNPs)
- **enrichplot** - GO enrichment visualizations
- **ggraph** - Network plot support
- **GAPIT3** (tag GAPIT3.5) - GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK

### Association Workflow
1. **EMMAX** - Mixed model with BN kinship (PCA covariates, full dataset)
2. **LFMM** - Latent factor model (trains on LD-pruned, tests on full)
3. **GAPIT3** - 8 additional models (GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK) via pre-computed GD/GM numeric format
4. **Combine** - Merge methods via `association.combine_method` (Sum/Overlap/single)
5. **Region clustering** - Per-trait and combined climate regions
6. **Gene annotation** - Genes within extended regions
7. **GO enrichment** - Per-region analysis with dotplot/emapplot/cnetplot

GAPIT models are auto-detected from config method names and routed through `gapit.R` with shared BN kinship + LEA PCA covariates. EMMAX/LFMM continue through their existing scripts. All methods produce standardized pvalue TSVs consumed by the same downstream rules.

### Phenotype Association Workflow (`association_phenotypes` mode)
1. **prepare_phenotypes** - Extract traits from metadata columns 5+, handle missing values (MEAN/MEDIAN/DROP)
2. **Per-trait VCF subsetting** (DROP mode) - Subset VCF to samples with non-missing trait values
3. **TPED/Kinship** - Separate Snakemake rules (not inside R script) for TPED conversion and BN kinship computation
4. **EMMAX** - `emmax_phenotypes.R` receives pre-computed TPED, kinship, whole-dataset PCA projections
5. **GAPIT3** - Same 8 models as GEA; Path A: single call with all traits; Path B: per-trait calls with sample subset
6. **Combine** (DROP mode) - `combine_pheno_pvalues.R` merges per-trait results (per method)
7. **Downstream** - Reuses existing scripts: find_sig_snps, create_regions, find_genes, enrichment, manhattan
- Two paths: Path A (MEAN/MEDIAN, single sample set) vs Path B (DROP, per-trait `{pheno_trait}` wildcard rules)
- PCA: uses whole-dataset PCA projections (subsetted to trait samples in `emmax_phenotypes.R`)
- Config: `phenotype_association.*` inherits from `association.*` by default, override only if different
- GAPIT DROP mode uses full-dataset GD + kinship; `gapit.R` subsets internally via `SAMPLES_SUBSET` arg

### association.configs Format
`association.configs` in config YAML is a list of `{method, adjust, threshold}` entries. Parsed into a dict mapping method name to `"adjust_threshold"` string. Controls which association methods run, their p-value adjustment, and significance thresholds.

### Imputation Strategy
- Uses sNMF Q-matrix (ancestry-informed)
- Happens twice: LD-pruned (for PCA/structure/climate) and full (for association)

### Path Management
- `W` - Working files (`_work/`, `_intermediate/`)
- `O` - Organized outputs (module-based: `GEA/plots/`, `PreStructure/tables/`, etc.)
- Module path constants: `MOD_PROCESSING`, `MOD_PRESTRUCT`, `MOD_CLIMATE`, `MOD_STRUCT`, `MOD_PREGEA`, `MOD_GEA`, `MOD_GWAS`, `MOD_GEAXGWAS`, `MOD_MALAD`
- Paths expand dynamically based on config parameters

## Snakefile Internals

### workdir
The main Snakefile includes `workflow/rules/common.smk` which sets `workdir: OUTDIR`, so all rule paths are relative to `{PROJECT}_results/`. The `W` and `O` dictionaries use absolute paths, so this mostly affects `shell:` directives that use relative paths.

### Path Dictionaries
- `W` dict: Working/intermediate file paths. Populated at parse time + dynamically by `add_kbest_paths()`, `add_association_paths()`, `add_maladaptation_paths()`.
- `O` dict: Organized output paths (plots/, tables/).
- Template functions (e.g., `clusters_table(k)`, `manhattan_plot(method, trait, adjust)`): For outputs parameterized by K, trait, or method.

To add a new output: add to the appropriate dict/function, add to `get_targets()` for the relevant mode, then create the rule.

### R Script Conventions
- Args parsed via `commandArgs(trailingOnly=TRUE)` positionally (args[1], args[2], etc.)
- Logging via `message()` (stderr, captured by Snakemake log)
- Libraries loaded at top; never installed in scripts
- `tryCatch()` for optional operations (enrichment plots that may fail with small data)
- All scripts executed as `Rscript /pipeline/scripts/<name>.R` inside Docker
- **CRITICAL: Always use `colClasses` when reading files with sample/site IDs via `fread()`**. Numeric sample IDs (e.g., `88`, `108`) cause `fread()` to infer integer type, breaking `%in%` and `left_join` against character VCF headers. Use `fread(..., colClasses = c("site" = "character", "sample" = "character"))` for metadata/sample files, or `fread(..., colClasses = "character")` for headerless sample lists.

### Known Quirks
- LEA's `vcf2lfmm()` silently strips "chr" prefix — this is why we normalize both VCF and GFF early
- LEA's `pca()` strips file extension and creates `{basename}.pca/` directory
- sNMF creates `.snmfProject` directory — must use `mode="new"` or `mode="continue"`
- `qvalue` can fail with small datasets — scripts fall back to `p.adjust`
- `emapplot` requires `pairwise_termsim()` called first (adds similarity matrix to enrichResult)
- `emapplot` and `cnetplot` need ≥2 enriched terms; `dotplot` works with 1+
- `download_climate_present` (WorldClim) is very slow — always preserve its output
- Dead `assoc_genes()`/`assoc_genes_collapsed()` template functions referencing undefined `GENE_DISTANCE` were removed — gene finding uses `find_genes_around_regions.R` directly
- `scattermore` is used via `association.scattermore_threshold` (default 30,000) to downsample non-sig SNPs in Manhattan plots
- `ld_prune` sed pattern uses `0_0_` while `filter_vcf` and `subset_vcf_pheno` use `0_` — inconsistent but confirmed working; likely because LD pruning goes through an extra plink step that doubles the prefix. Investigate only if errors arise.

## Shiny App — golem Package (`scripts/clinego.app/`)

Interactive results viewer built as a **golem R package** using bslib (Bootstrap 5). The legacy `scripts/app.R` is preserved as reference but is not in active use.

### Dev Mode (no Docker rebuild)

```bash
docker run --user $(id -u):$(id -g) --rm -e USER=pipeline -p 3838:3838 -v $PWD:/pipeline cline-go:latest \
  Rscript /pipeline/scripts/clinego.app/dev.R
```

`dev.R` sources all `R/*.R` files from the mounted volume at startup. Docker rebuild only needed when adding new R package dependencies to DESCRIPTION.

**Always restart the container after any R file change.** Shiny autoreload does not work reliably when files are modified from outside the container (inotify events from the host volume mount are not forwarded). Use `docker stop $(docker ps -q --filter ancestor=cline-go:latest) && docker run ...` — the restart takes ~10s and is the only reliable way to pick up changes.

**Input persistence in `renderUI` (UI rule):** When user-editable inputs live inside `renderUI`, they reset to their default value on every re-render. **Never silently discard a user's chosen parameter value.** For any input the user can modify that lives inside a `renderUI`:
- Store the user's value in a `reactiveVal` (captured with `observeEvent(..., ignoreInit = TRUE)`)
- Use the stored value as the input's `value` in the `renderUI` (fallback to config default when NULL)
- Show a visual "modified" indicator (amber `badge bg-warning` badge) when the value differs from **saved params** (`region_params.json`), NOT the config default. The badge means "current value doesn't match what produced the visible output". No badge when no computation has been done yet for the region.
- This rule applies to any parameter where users pick a specific value to run a computation (epsilon, distance thresholds, etc.)

**Exploratory parameter persistence (UX rule):** Exploratory parameters (region_distance, hap scan/viz params) are stored per-region in `{PROJECT}_results/_intermediate/region_params.json`, separate from the pipeline config YAML. When a computation (haplotype scan/viz) completes successfully, save the params used. When a region is selected, load saved params as input defaults. Config defaults shown only for new regions with no prior computations. Use `read_region_params()`, `save_region_params()`, `get_region_param()`, `set_region_param()`, `get_global_param()`, `set_global_param()` from `fct_region_params.R`. Global params (like `region_distance`) are keyed by module name. Complex parameter groups should have a "Reset to defaults" button (`btn-link text-muted`) that restores config YAML defaults.

**Playwright project switching:** The selectize.js project dropdown cannot be set with `Shiny.setInputValue()` or `playwright-cli select` — it ignores programmatic value changes. To switch projects in playwright-cli: click the dropdown container (`e14`), wait for the snapshot to reveal the option refs, then click the option ref directly (e.g. `playwright-cli click e274`).

### Visual Testing

**Playwright is allowed and encouraged for layout/visual checks — no need to ask per session.**
(Reversed 2026-08-14: the old blanket ban cost more than it saved. One 40-line script caught two
layout bugs — a flex-shrunk nav bar and a 12px page-level horizontal scroll — that parse checks,
Sass compilation and HTML inspection all passed clean.)

**How to invoke:** the `playwright-cli` wrapper named in the global CLAUDE.md is *not* on PATH on
this machine. Use plain `playwright` (1.58.0) or, for anything involving clicks and assertions,
the Python API — `python3` + `from playwright.sync_api import sync_playwright`. Browsers are
cached at `~/.cache/ms-playwright/`.

**What it is good for:** anything measurable in the DOM or visible in a screenshot — element
geometry, overflow/scroll behaviour at several viewport widths, computed styles, which pane is
active, whether a control is visible. Prefer *asserting* over eyeballing: screenshot for the look,
`page.evaluate()` for the facts.

**What it is still bad for:** deep stateful flows. The app's reactive dependencies and selectize.js
dropdowns make long scripted journeys brittle (a known one: the project dropdown ignores
`Shiny.setInputValue()` — you must click the container, wait for the option refs, then click the
option). Keep scripts short and single-purpose; do not build a UI regression suite out of them.

**Gotcha — scope your selectors.** The app nests navsets inside module panels, so a bare
`.tab-content > .tab-pane.active` matches an inner tabset, not the main one. Read
`#main_tabs`'s `data-tabsetid` and query
`.tab-content[data-tabsetid="<id>"] > .tab-pane.active` instead.

The user still reviews the actual design — send screenshots and wait for the verdict on anything
aesthetic. Playwright replaces the *mechanical* half of that loop, not the judgement half.

### Architecture

**Framework**: golem + bslib (Bootstrap 5). All plots served as static PNG/SVG images — no `.qs` dependencies.

**Manhattan performance**: Static background PNG (**all SNPs, all methods**, via scattermore, generated by pipeline) + lightweight plotly overlay (sig SNPs only, ~100-500 points). Pipeline outputs `_background.png` + `_coords.json` alongside normal plots. The background is the complete, static Manhattan cloud — the overlay only adds interactive highlight markers on top. Toggling a trait/method in the matrix removes the overlay marker; a faint ghost dot always remains in the background.

**Module files** (`R/`):
- `app_ui.R`, `app_server.R`, `run_app.R` — top-level app entry point
- `app_theme.R` — bslib Bootstrap 5 theme ("Precision Genomics": primary #1B7A6E, navbar #1A2332)
- `app_config.R` — `get_pipeline_path()` (option → golem-config.yml → env var → /pipeline)
- `fct_paths.R` — all output path functions, MOD_* constants
- `fct_discovery.R` — runtime discovery: `find_projects()`, `find_k_values()`, `find_assoc_methods()`, etc.
- `fct_config.R` — YAML config reader, `config_get()`, `%||%`
- `fct_data_loading.R` — TSV loaders, `load_cached()` with 200MB cachem cache
- `fct_manhattan.R` — `build_manhattan_plotly()`, coord alignment, `add_cum_pos()`
- `fct_labels.R` — `build_region_labels()`, `format_region_id()`, `format_hap_tag()`
- `utils_ui.R` — `plot_placeholder()`, `region_info_bar()`, `mode_status_row()`
- `utils_helpers.R` — `file_ok()`, `resolve_adjust()`, `make_project_data()`, `safe_datatable()`
- `mod_image_card.R` — reusable static image card (full-screen + SVG/PNG download)
- `mod_manhattan_overlay.R` — plotly overlay on background PNG; returns `selected_region` reactive
- `mod_region_detail.R` — GO enrichment, genes, GO table, regionplot, haplotype (used by 3 tabs)
- `mod_piemap_viewer.R` — piemap display with bio/metric/zoom path resolution
- `mod_home.R`, `mod_structure.R`, `mod_structure_k.R`, `mod_maladaptation.R` — simple tabs
- `mod_association.R`, `mod_phenotype.R`, `mod_overlapping.R`, `mod_haplotype.R` — complex tabs

**Region-centric design**: Clicking a sig SNP in any Manhattan → selects region → `mod_region_detail` renders below (enrichment plots, genes table, GO table, regionplot, haplotype viz). Same pattern in Association, Phenotype Association, and Overlapping tabs.

**Haplotype tag resolution**: Tags are `{meta_type}_{source}` (e.g., `site_association`). Each association tab finds its matching haplotype tag by splitting on `_` and matching source part.

**Plotly source scoping**: Each `mod_manhattan_overlay` instance uses `ns("overlay")` as the plotly event source to prevent click events from cross-firing between Combined and Per-Method Manhattans.

### Piemap Pie/Points Toggle

Dense sampling makes pie charts overlap and occlude the raster/geography. Every geo piemap (Structure ancestry, Maladaptation genetic-offset, GWAS phenomap) is generated **twice** by the pipeline: the pie chart (main render) and a `*_points` companion (tiny dark dots marking sample locations only, via `plot_piemap.R`'s `use_points` branch). The Shiny app renders both and lets the user flip between them at runtime with a `bslib::input_switch("Points")` in each tab's `.control-bar` — no pipeline re-run needed.

- Points are **trait/metric-independent** (they only plot lon/lat) — one points file per background raster (per bio for Structure, per method/run_label/spatial_tag for Maladaptation, per project for GWAS phenomap), not per metric/trait/variant.
- `piemap_path()` and `pheno_piemap_path()` (`fct_paths.R`) take a `points = FALSE` param that resolves the companion filename and ignores `metric`/`trait` when TRUE. `gf_offset_piemap_path(..., variant = "points")` already works via its verbatim-variant branch — no change needed there.
- Structure's `_no_spatial_variance.flag` placeholder (climate var invariant across sites) is **bypassed** in points mode — points show geography, which stays meaningful even when the climate signal doesn't vary spatially.
- Maladaptation's points toggle ignores the zoom/variant selectors (the flat `zoom/{tag}.png` naming has no points companion) and always shows the project-level points file.
- `Piemap.use_points` config key still exists but only governs the on-demand haplotype-viz piemaps (`fct_regions.R`) — it is NOT read by the Structure/Maladaptation/GWAS map rules, which always emit both renders.

### Piemap Sizing

Pipeline outputs piemaps at diverse aspect ratios — from small regional maps (e.g., Israel) to global projections (e.g., Arabidopsis 1001 genomes). The app handles this via CSS:
- `max-height: 70vh` prevents tall maps from pushing content off-screen
- `object-fit: contain` preserves aspect ratio in full-screen mode
- `width: auto; max-width: 100%` — natural size up to card width, no stretching
- No fixed pixel dimensions — the pipeline determines the map extent via `map.climate_extent`

When adding new map-type images, wrap the card in `htmltools::div(class = "piemap-container", ...)`.

### Conditional Content Display

When pipeline modules produce optional outputs (pop_stats piemaps, zoom maps, haplotype viz), the app hides controls/sections that have no data rather than showing empty placeholders:
- Build selector `choices` dynamically by checking `file_ok()` on each variant path
- If only one variant exists, hide the selector entirely (`return(NULL)`)
- Prefer "no selector" over "selector with one disabled option"
- Use `plot_placeholder(message, suggestion)` — the `suggestion` param adds a smaller line explaining which pipeline mode to run

**CRITICAL — selectInput inside cards**: Never wrap `selectInput` inside `bslib::card()` without `overflow: visible` on the card. Bootstrap 5 cards have `overflow: hidden` by default (for border-radius clipping), which hides selectize.js dropdown menus. The `.control-bar` CSS class in `custom.scss` handles this automatically — always use that class for inline control bars.

### Running the App

Two run paths, both verified 2026-09-10. Prefer **dev mode** for day-to-day work — it needs no
rebuild after an R file change. Use the **package mode** only to check that the installed
package itself is sound.

**Dev mode** (file path — `dev.R` sources `R/*.R` off the mount, `.onLoad()` never fires):
```bash
docker run --user $(id -u):$(id -g) --rm --name clinego_app -e USER=pipeline -p 3838:3838 \
  -v $PWD:/pipeline cline-go:latest Rscript /pipeline/scripts/clinego.app/dev.R
```

**Package mode** (`library(clinego.app)` → `.onLoad()` in `zzz.R` loads the shared libs):
```bash
docker run --user $(id -u):$(id -g) --rm --name clinego_prod -e USER=pipeline -p 3838:3838 \
  -v $PWD:/pipeline cline-go:latest R -e "clinego.app::run_app(host = '0.0.0.0', port = 3838)"
```

**Pass `host`/`port` directly — NOT `options = list(...)`.** `run_app(...)` already wraps its
`...` in `shinyApp(options = list(...))`, so the old documented form nested to
`options = list(options = list(host, port))`, which Shiny ignores: the app came up on a random
loopback port (observed: `127.0.0.1:7729`) and the published `-p` mapping went nowhere.

`-v $PWD:/pipeline` is required for **both** paths: the image contains no pipeline code (see the
marker comment at the app-install block in `Dockerfile`).

### Key Config Files
- `scripts/clinego.app/inst/golem-config.yml` — sets `pipeline_path: /pipeline`
- `scripts/clinego.app/DESCRIPTION` — package metadata + Imports

### Legacy App
- `scripts/app.R` — shinydashboard monolith (3,466 lines). Preserved for reference. Uses OLD flat path structure and `.qs` files — do NOT use for new development.

## TODO

### Unit Testing

Two testthat roots. **Both must be green before a merge.** One command runs them:

```bash
./tests/run_all.sh                                # the merge gate
./tests/run_all.sh --invariants SIMDATA_results   # + validate a results tree
```

It runs every suite to completion (never stops at the first failure), prints one line per
suite and exits non-zero if any failed. The two suites individually:

```bash
# Tier 1 — scripts/R/lib + scripts/R/utils (the shared science). ~20 s.
docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -v $PWD:/pipeline \
  cline-go:latest Rscript /pipeline/tests/run_tests.R

# The Shiny app package. ~40 s.
docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -v $PWD:/pipeline \
  cline-go:latest Rscript -e 'setwd("/pipeline/scripts/clinego.app/tests"); source("testthat.R")'
```

Baseline as of 2026-09-10 (after Tier 5): `tests/` = **567 passing / 14 skipped**, app =
**213 passing / 1 skipped**. `run_tests.R` exits non-zero on any failure, so it is CI-able
as-is. (The earlier figure of 473 recorded here was stale — it predated Tier 4; the counts
reconcile against the dossier's post-Tier-4 494 plus Tier 5's 73.)

**Tier 5 — app/pipeline equivalence** lives in
`tests/testthat/test-equivalence-app-pipeline.R` (33 tests, 73 assertions, 4 `skip()`ped).
It asserts the Shiny app and the Snakemake pipeline compute the same science across the
**five** divergence surfaces, and is the only place that does. Two facts about it matter
before editing:

- **Both sides load in one session without colliding**, which is what makes the file
  possible: `zzz.R:57-60` sources the shared libs into `asNamespace("clinego.app")` while
  `helper-libs.R:28-41` sources the pipeline libs into the test env. `combine_sigsnps`
  exists on **both** sides with different arity, so the app's must always be written
  `clinego.app:::combine_sigsnps`. Do **not** add `library(clinego.app)` to
  `helper-libs.R` — `test_dir()` sources every helper for the directory, so the whole
  pipeline suite would fail whenever the app is not installed; the skip is file-local.
- **The `qval` equivalence is asserted structurally, not live, on purpose.**
  `helper-libs.R:14` attaches `qvalue` and namespace lookup falls through to the search
  path, so the app's qval branch *works under this suite* while returning `NA` in
  production. The file asserts the cause instead — that `qvalue` is absent from the app's
  DESCRIPTION Imports. A live qval comparison there would pass and prove nothing.

**`--invariants` is opt-in and is EXPECTED to be red**, which is why it is not part of the
default gate. `scripts/check_invariants.R` validates a `{PROJECT}_results/` tree against
`scripts/R/lib/invariants.R` — checks that must hold on ANY dataset (region/SNP consistency,
p-values in [0,1], chromosome names never re-acquiring a `chr` prefix, sample accounting that
closes, and cross-module referential integrity). On `SIMDATA_results` it reports 41 violations,
**every one of them a defect already filed in `docs/pipeline_improvement_requests.md`**:

| count | check | filed as |
|---|---|---|
| 14 | `overlap_traits_includes_own_trait` | `sig_snps.R:152` |
| 14 | `overlap_snps_names_unknown_snp` | `sig_snps.R:159` |
| 10 | `duplicate_column_names` | `genes_per_region_collapsed.tsv` |
| 1 | `pairwise_table_references_unknown_trait` | stale GEAxGWAS after a `mode=gea` re-run |
| 1 | `multiple_threshold_variants_on_disk` | RDA sig table orphaned by a threshold change |
| 1 | `climate_predictor_count_disagrees` | `write_summary.R:362` |

A run that comes back **clean is itself a failure signal** — those six are confirmed present.
Any check name outside that table is new and must be triaged: a real defect gets filed, a
checker misreading a schema gets fixed. Never tune a check down to make the output green.

The checkers are pure (data.tables in, a violations data.table out) and unit-tested against
hand-built fixtures in `tests/testthat/test-invariants.R` — each one paired: a clean fixture
that must return zero violations and a broken fixture that must return exactly the expected
violation. Those tests are **not** quarantined and must never be skipped; they assert that new
code fires, not that a defect exists.

`tests/` is a **non-package** root: it uses `test_dir()` plus
`tests/testthat/helper-libs.R`, which attaches the packages the libs assume (they never
call `library()` themselves — several need a bare `%>%`, `qvalue()` or `covRob()`) and
`source()`s them in dependency order. Add a new lib to that vector when you add one.
Both commands need `-v $PWD:/pipeline`; the image ships no pipeline code. Do **not** route
either through `dev.R` or prepend `.R_libs_dev` — that carries testthat 3.3.2 while the
image pins 3.2.3.

`tests/testthat/test-known-bugs.R` holds correct-behaviour assertions for defects that are
known and deliberately unfixed, each behind a `skip()` naming where it is filed. Fixing one
means deleting a `skip()` line. Never weaken a test there to match current output. The app
suite now uses the same convention — `test-fct_threshold_rules.R` carries one `skip()`ped
correct-behaviour assertion for the qvalue-not-in-Imports defect.

Still missing: **Python tests for `scripts/*.py`** (`design_adequacy.py`, `gff2topr.py`,
`snakemake_progress_handler.py` — no Python test infrastructure exists at all; note the image
has python3.12 + numpy + stdlib `unittest` but **no pytest**, so `unittest` needs no Dockerfile
change), golden-file regression on SIMDATA outputs, and Shiny/pipeline equivalence checks.

**On golden files specifically**: they are blocked, not merely undone. Three committed SIMDATA
outputs are known-wrong and would be frozen as "expected" (`overlap_traits`/`overlap_snps`,
exon/promoter counts), and `find_significant_snps_per_trait()`'s diagnostics are
**core-count-dependent** — 2 rows at `cpu=1`, 0 at `cpu>=2` (`sig_snps.R:44`) — so that output
cannot be a golden at all. Fix or explicitly quarantine those first.

**On Shiny/pipeline equivalence**: the surface is larger than it looks. `zzz.R:39-42` shares
only `regions.R` and `pval_threshold.R` with the pipeline; gene finding (`fct_regions.R:67` vs
`genes_in_regions.R`), combining (`fct_combine.R` vs `combine_sigsnps.R`) and region distance
are **reimplemented** in the app, and nothing asserts the two agree.

See the "Regression tests for scientific outputs" objective in the ADAPTOGENE pipeline dossier
for the full tier plan.

### Exon/Promoter SNP Validation — validated 2026-09-10, and it is WRONG
`.count_snps_in_features()` (`scripts/R/lib/genes_in_regions.R:174`) builds its SNP id from
the plain, feature-side `s` column of `foverlaps(snps, feats)`, so `exon_snps` /
`promoter_snps` report the **feature's start coordinate** and `exon_snp_count` /
`promoter_snp_count` count features hit, not SNPs. Two SNPs at 120 and 130 inside a
100-200 exon yield `1:100` and a count of 1. Invisible on SIMDATA (every cell is empty
there). Not fixed — filed in `docs/pipeline_improvement_requests.md` and quarantined in
`tests/testthat/test-known-bugs.R`.

## Active Obsidian Project
- Project: ADAPTOGENE
- File: ~/Orthidian/projects/ADAPTOGENE/ADAPTOGENE.md

## Permission Guidelines

### Allowed (no approval needed)
- `docker build` or `docker run` commands
- Reading any project files
- Writing/editing: `Snakefile`, `scripts/*.R`, `scripts/*.py`, `scripts/clinego.app/R/*.R`, `config*.yaml`, `CLAUDE.md`, `Dockerfile`
- Writing/editing test data: `data/SIMDATA*`
- Removing individual output files (prefer `-R`/`--forcerun`)

### Requires permission
- Removing entire `*_results/` or `*_logs/` directories
- Git operations (user handles manually)
- Modifying `Snakefile_old`
- Installing new R packages or changing Dockerfile versions

## graphify — DISABLED for this project

Do not use graphify (`graphify-out/`, `graphify query/path/explain/update`) in CLINE-GO. Use normal Grep/Glob/Explore instead, even for architecture/cross-module questions. `graphify-out/` may still exist on disk — ignore it, do not read `GRAPH_REPORT.md` or `wiki/index.md`, do not run `graphify update` after edits.
