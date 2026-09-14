# Shiny App Path Migration Guide

The pipeline restructuring (Phases 1-4) changed output directory structure, file naming, region ID format, and config keys. The Shiny app (`scripts/app.R`) still references old paths and will not work until updated.

This document maps every old path to its new equivalent.

---

## 1. Directory Structure Changes

All paths are relative to `{PROJECT}_results/`.

### Old → New Base Directories

| Old Path | New Path | Notes |
|---|---|---|
| `plots/` | `{module}/plots/` | Split by module |
| `tables/` | `{module}/tables/` | Split by module |
| `tables/association/` | `association/tables/` | |
| `tables/association_phenotypes/` | `phenotype_association/tables/` | Mode renamed |
| `tables/overlapping/` | `overlapping/tables/` | |
| `tables/structure/` | `processing/tables/` | metadata.tsv moved to processing |
| `tables/gradientForest/` | `maladaptation/tables/` | |
| `tables/haplotype_{tag}/` | `haplotype/{tag}/tables/` | Tag no longer prefixed |
| `plots/haplotype_{tag}/` | `haplotype/{tag}/plots/` or `haplotype_scan/{tag}/plots/` | Split scan vs viz |
| `intermediate/` | `_intermediate/` | Underscore prefix |
| `intermediate/haplotype_{tag}/` | `_intermediate/haplotype/{tag}/` | |
| `work/` | `_work/` | Underscore prefix |

### Module Mapping

| Analysis Mode | Old Paths | New Module |
|---|---|---|
| Processing | `tables/structure/metadata.tsv` | `processing/tables/` |
| Structure | `plots/{pca,tracy_widom,cross_entropy,structure_K,pca_structure_K,pop_diff_K}.*` | `structure/plots/` |
| Structure K - Climate | `plots/climate_{res}/`, `intermediate/climate/` | `climate/{plots,tables,rasters}/` |
| Structure K - Piemaps | `plots/{PieMap,piemap}*` | `structure_k/plots/piemap/` |
| Structure K - Pop Stats | `plots/{MantelTest,AMOVA}*`, `tables/{TajimaD,Pi_diversity,IBD}*` | `structure_k/{plots,tables}/pop_stats/` |
| Association (GEA) | `plots/{EMMAX,LFMM}/`, `tables/association/` | `association/{plots,tables}/` |
| Phenotype Association | `plots/association_phenotypes/`, `tables/association_phenotypes/` | `phenotype_association/{plots,tables}/` |
| Overlapping | `tables/overlapping/` | `overlapping/{plots,tables}/` |
| Regionplot | `plots/regionplot_*` | `regionplot/plots/` |
| Maladaptation | `plots/{CumulativeImportance,OverallImportance,GeneticOffsetPieMap}*` | `maladaptation/plots/` |
| Haplotype Scan | `plots/haplotype_{tag}/clustree/` | `haplotype_scan/{tag}/plots/clustree/` |
| Haplotype Viz | `plots/haplotype_{tag}/Region_*` | `haplotype/{tag}/plots/` |

---

## 2. File Naming Changes (CamelCase → snake_case)

### Tables

| Old Filename | New Filename | Used In |
|---|---|---|
| `Selected_SNPs.tsv` | `selected_snps.tsv` | Association, Phenotype, Overlapping |
| `Regions_climate_combined.tsv` | `regions_combined.tsv` | Association |
| `Regions_phenotype_combined.tsv` | `regions_combined.tsv` | Phenotype Association |
| `Regions_all_combined.tsv` | `regions_combined.tsv` | Overlapping |
| `Regions_per_trait.tsv` | `regions_per_trait.tsv` | Association, Phenotype |
| `Regions_per_trait_all.tsv` | `regions_per_trait_all.tsv` | Overlapping |
| `Genes_per_region.tsv` | `genes_per_region.tsv` | All association modes |
| `Genes_per_region_collapsed.tsv` | `genes_per_region_collapsed.tsv` | All association modes |
| `Genes_per_region_combined.tsv` | `genes_combined.tsv` | Phenotype, Overlapping |
| `Overlap_summary.tsv` | `overlap_summary.tsv` | Overlapping |
| `Selected_SNPs_all.tsv` | `selected_snps_all.tsv` | Overlapping |
| `Pipeline_summary.tsv` | `pipeline_summary.tsv` | Top-level |
| `GeneticOffset_site_{sfx}.tsv` | `genetic_offset_site_{sfx}.tsv` | Maladaptation |
| `TajimaD_byPop.tsv` | `tajima_d_by_pop.tsv` | Structure K (pop stats) |
| `Pi_diversity_byPop.tsv` | `pi_diversity_by_pop.tsv` | Structure K (pop stats) |
| `IBD_raw.tsv` | *(removed 2026-09-13 — per-pair Mantel was not an IBD test; see mantel_test.tsv)* | Structure K (pop stats) |
| `IBD_notIsolated.tsv` | *(removed 2026-09-13)* | Structure K (pop stats) |
| `AMOVA.tsv` | `amova.tsv` | Structure K (pop stats) |
| `metadata.tsv` | `metadata.tsv` | Unchanged name, moved to `processing/tables/` |

### Enrichment Tables

| Old Pattern | New Pattern |
|---|---|
| `Region_{rid}_enrichment.tsv` | `region_{rid}_enrichment.tsv` |
| `Region_{rid}_{trait}_enrichment.tsv` | `region_{rid}_{trait}_enrichment.tsv` |
| `Enrichment_{rid}_summary.tsv` | `enrichment_{rid}_summary.tsv` |

### Enrichment Table Directory

| Old Path | New Path |
|---|---|
| `tables/association/enrichment/{trait}/` | `association/tables/enrichment/{trait}/` |
| `tables/association_phenotypes/enrichment/{trait}/` | `phenotype_association/tables/enrichment/{trait}/` |
| `tables/overlapping/enrichment/{trait}/` | `overlapping/tables/enrichment/{trait}/` |

### Association P-value Tables

| Old Pattern | New Pattern | Notes |
|---|---|---|
| `tables/association/{METHOD}/{METHOD}_pvalues_K{k}.tsv` | `association/tables/{METHOD}/{method}_pvalues_K{k}.tsv` | Method name lowercased in filename |
| `tables/association/{METHOD}/{METHOD}_pvalues_K{k}_sig_snps_{adjust}.tsv` | `association/tables/{METHOD}/{method}_pvalues_K{k}_sig_snps_{adjust}.tsv` | |
| `tables/association_phenotypes/{method}_phenotypes_pvalues_K{k}.tsv` | `phenotype_association/tables/{METHOD}/{method}_pvalues_K{k}.tsv` | Removed "_phenotypes" infix |

---

## 3. QS Plot File Changes

### Structure Plots (unchanged names, new location)

| Old Path | New Path |
|---|---|
| `plots/pca.qs` | `structure/plots/pca.qs` |
| `plots/tracy_widom.qs` | `structure/plots/tracy_widom.qs` |
| `plots/cross_entropy_K{s}-{e}.qs` | `structure/plots/cross_entropy_K{s}-{e}.qs` |
| `plots/pca_structure_K{k}.qs` | `structure/plots/pca_structure_K{k}.qs` |
| `plots/pop_diff_K{k}.qs` | `structure/plots/pop_diff_K{k}.qs` |
| `plots/structure_K{k}.qs` | `structure/plots/structure_K{k}.qs` |

### Climate Plots

| Old Pattern | New Pattern |
|---|---|
| `plots/CorrelationHeatmap_present.qs` | `climate/plots/correlation_heatmap.qs` |
| `plots/DensityPlot_bio_{N}.qs` | `climate/plots/density_plot_bio_{N}.qs` |
| `plots/DensityPlot_present.qs` | `climate/plots/density_plot_present.qs` |

### Piemap Plots

| Old Pattern | New Pattern |
|---|---|
| `plots/PieMap_bio_{N}.qs` | `structure_k/plots/piemap/piemap_bio_{N}.qs` |
| `plots/PieMap_bio_{N}_TajimaD.qs` | `structure_k/plots/piemap/tajima_d/piemap_bio_{N}_tajima_d.qs` |
| `plots/PieMap_bio_{N}_PiDiversity.qs` | `structure_k/plots/piemap/pi_diversity/piemap_bio_{N}_pi_diversity.qs` |
| `plots/PhenoMap_{trait}.qs` | `phenotype_association/plots/piemap/phenomap_{trait}.qs` |

### Pop Stats Plots

| Old Pattern | New Pattern |
|---|---|
| `plots/MantelTest.qs` | `structure_k/plots/pop_stats/mantel_test.qs` |
| `plots/AMOVA.qs` (or `AMOVA_plot.qs`) | `structure_k/plots/pop_stats/amova_plot.qs` |

### Manhattan Plots (QS)

| Old Pattern | New Pattern |
|---|---|
| `plots/{METHOD}/Manhattan_{trait}_K{k}_{adjust}.qs` | `association/plots/manhattan/{METHOD}/manhattan_{trait}_K{k}_{adjust}.qs` |
| `plots/Manhattan_all_traits_combined_K{k}.qs` | `association/plots/manhattan_combined_K{k}.qs` |
| `plots/association_phenotypes/{METHOD}/Manhattan_{trait}_K{k}_{adjust}.qs` | `phenotype_association/plots/manhattan/{METHOD}/manhattan_{trait}_K{k}_{adjust}.qs` |

### Enrichment Plots (QS)

| Old Pattern | New Pattern |
|---|---|
| `plots/association/enrichment/{trait}/Region_{rid}_{type}.qs` | `association/plots/enrichment/{trait}/region_{rid}_{type}.qs` |
| `plots/association_phenotypes/enrichment/{trait}/Region_{rid}_{type}.qs` | `phenotype_association/plots/enrichment/{trait}/region_{rid}_{type}.qs` |
| `plots/overlapping/enrichment/{trait}/Region_{rid}_{type}.qs` | `overlapping/plots/enrichment/{trait}/region_{rid}_{type}.qs` |

### Maladaptation Plots

| Old Pattern | New Pattern |
|---|---|
| `plots/OverallImportance_{sfx}.qs` | `maladaptation/plots/overall_importance_{sfx}.qs` |
| `plots/CumulativeImportance_{sfx}.qs` | `maladaptation/plots/cumulative_importance_{sfx}.qs` |
| `plots/GeneticOffsetPieMap_{sfx}.qs` | `maladaptation/plots/genetic_offset_piemap_{sfx}.qs` |
| `plots/GeneticOffsetPieMap_{sfx}_TajimaD.qs` | `maladaptation/plots/genetic_offset_piemap_{sfx}_tajima_d.qs` |
| `plots/GeneticOffsetPieMap_{sfx}_PiDiversity.qs` | `maladaptation/plots/genetic_offset_piemap_{sfx}_pi_diversity.qs` |
| `plots/GeneticOffsetPieMap_{sfx}_zoom_{tag}.qs` | `maladaptation/plots/zoom/genetic_offset_piemap_{sfx}_zoom_{tag}.qs` |

### Haplotype Plots

| Old Pattern | New Pattern |
|---|---|
| `plots/haplotype_{tag}/clustree/Region_{rid}_clustree_MG.qs` | `haplotype_scan/{tag}/plots/clustree/region_{rid}_clustree_MG.qs` |
| `plots/haplotype_{tag}/clustree/Region_{rid}_clustree_hap.qs` | `haplotype_scan/{tag}/plots/clustree/region_{rid}_clustree_hap.qs` |
| `plots/haplotype_{tag}/Region_{rid}_crosshap_viz.qs` | `haplotype/{tag}/plots/region_{rid}_crosshap_viz.qs` |
| `plots/haplotype_{tag}/Region_{rid}_boxplot_{trait}.qs` | `haplotype/{tag}/plots/region_{rid}_boxplot_{trait}.qs` |
| `plots/haplotype_{tag}/HaplotypePieMap_Region_{rid}_{trait}.qs` | `haplotype/{tag}/plots/haplotype_piemap_region_{rid}_{trait}.qs` |

### Intermediate QS Files

| Old Pattern | New Pattern |
|---|---|
| `intermediate/haplotype_{tag}/Region_{rid}_HapObject.qs` | `_intermediate/haplotype/{tag}/region_{rid}_hapobject.qs` |
| `intermediate/AMOVA_result.qs` | `_intermediate/amova_result.qs` |

---

## 4. Region ID Format Change

Old format: `chr:start-end` (e.g., `1:50000-60000`)
New format: `chr_start-end` (e.g., `1_50000-60000`)

### Affected Regex Patterns

| Location | Old Regex | New Regex |
|---|---|---|
| Region ID parsing | `([^:]+):(\\d+)-(\\d+)` | `([^_]+)_(\\d+)-(\\d+)` |
| Region dropdown labels | Contains `:` in ID | Contains `_` in ID |
| Enrichment file lookup | `Region_{chr}:{start}-{end}_...` | `region_{chr}_{start}-{end}_...` |
| Haplotype file lookup | `Region_{chr}:{start}-{end}_...` | `region_{chr}_{start}-{end}_...` |

### Region ID Columns in Tables

The `region_id` column in these tables now uses `_` instead of `:`:
- `regions_per_trait.tsv`
- `regions_combined.tsv`
- `genes_per_region.tsv`
- `genes_per_region_collapsed.tsv`
- `selected_regions.tsv` (haplotype)

---

## 5. Config Key Changes

The Shiny app reads config values via `read_project_config()` which uses flat key regex matching.

| Old Key | New Key (nested YAML) | Notes |
|---|---|---|
| `PROJECT_NAME: X` | `project_name: X` | Lowercase |
| `ASSOC_TOP_REGIONS: N` | `association.top_regions` | Nested under `association:` |
| `PHENO_TOP_REGIONS: N` | `phenotype_association.top_regions` | Inherits from `association.top_regions` if not set |
| `OVERLAP_TOP_REGIONS: N` | `overlap.top_regions` | Optional override |

The `read_project_config()` function (line 326) uses `grep(paste0('^', key, ':'), lines)` which only works with flat top-level keys. It must be rewritten to parse nested YAML (use `yaml::read_yaml()` instead).

---

## 6. Discovery Function Changes

### `find_projects()` (line 76)
- Pattern `list.dirs(pipeline_path, ...)` looking for `*_results` — **no change needed**

### `find_k_values()` (line 19)
- Scans `intermediate/` for `*_K{k}*.qs` — needs new path: `structure/plots/` or `_intermediate/`
- Pattern: looks for structure plot filenames

### `find_bio_values()` (line 36)
- Scans `intermediate/` for `PieMap_bio_{N}.qs` — needs:
  - New path: `structure_k/plots/piemap/`
  - New pattern: `piemap_bio_{N}.qs`

### `find_metric_values()` (line 52)
- Scans `intermediate/` for `PieMap_bio_*_{Metric}.qs` — needs:
  - New path: `structure_k/plots/piemap/`
  - New pattern: subdirectories `tajima_d/`, `pi_diversity/`
  - Old suffixes `TajimaD`, `PiDiversity` → `tajima_d`, `pi_diversity`

### `find_gf_suffixes()` (line 96)
- Scans `intermediate/` for `GeneticOffsetPieMap_*.qs` — needs:
  - New path: `maladaptation/plots/`
  - New pattern: `genetic_offset_piemap_*.qs`

### `find_gf_zooms()` (line 108)
- Scans `intermediate/` for `GeneticOffsetPieMap_*_zoom_*.qs` — needs:
  - New path: `maladaptation/plots/zoom/`
  - New pattern: `genetic_offset_piemap_*_zoom_*.qs`

### `find_haplotype_tags()` (line 119)
- Scans `intermediate/` dirs matching `haplotype_*` — needs:
  - New path: `_intermediate/haplotype/` (subdirs are tag names directly, no `haplotype_` prefix)
  - Check for `scan_done.flag` in `haplotype_scan/{tag}/` instead of `intermediate/haplotype_{tag}/`
  - HapObject pattern: `region_*_hapobject.qs` (was `Region_*_HapObject.qs`)

### `find_hap_region_ids()` (line 180)
- Scans `intermediate/haplotype_{tag}/` for `Region_*_HapObject.qs` — needs:
  - New path: `_intermediate/haplotype/{tag}/`
  - New pattern: `region_*_hapobject.qs`

### `find_hap_traits()` (line 188)
- Scans `plots/haplotype_{tag}/` for boxplot files — needs:
  - New path: `haplotype/{tag}/plots/`
  - New pattern: `region_{rid}_boxplot_*.qs`

### `load_plot()` (line 485)
- Searches both `intermediate/` and `plots/` recursively for `.qs` files — needs:
  - Search `_intermediate/` and all `{module}/plots/` directories
  - Or: search entire `{PROJECT}_results/` tree (simpler but slower)

### `load_enrichment_plot()` (line 542)
- Base path: `plots/{mode}/enrichment/{trait}/` — needs:
  - New: `{module}/plots/enrichment/{trait}/`
  - `association_phenotypes` → `phenotype_association`
  - Filename pattern: `Region_{rid}_{type}.qs` → `region_{rid}_{type}.qs`

### `load_assoc_data()` / `load_pheno_assoc_data()` / `load_overlap_data()`
- All base paths change (see Section 1)
- All filenames change (see Section 2)
- `list.dirs()` for method discovery in association: path change only

### `check_haplotype_available()` (line 598)
- Path: `intermediate/haplotype_{tag}/Region_{rid}_HapObject.qs`
- New: `_intermediate/haplotype/{tag}/region_{rid}_hapobject.qs`

---

## 7. `load_plot()` Pattern Changes

The `load_plot()` function matches `.qs` files by regex. All calling code passes patterns that must be updated:

| Old Pattern | New Pattern | Line(s) |
|---|---|---|
| `'^pca\\.qs$'` | `'^pca\\.qs$'` | 1354 (unchanged) |
| `'tracy_widom\\.qs$'` | `'tracy_widom\\.qs$'` | 1358 (unchanged) |
| `'cross_entropy_K[0-9]*-[0-9]*\\.qs$'` | same | 1362 (unchanged) |
| `'pca_structure_K{k}\\.qs$'` | same | 1371 (unchanged) |
| `'pop_diff_K{k}\\.qs$'` | same | 1377 (unchanged) |
| `'^structure_K{k}\\.qs$'` | same | 1383 (unchanged) |
| `'CorrelationHeatmap.*\\.qs$'` | `'correlation_heatmap\\.qs$'` | 1392 |
| `'DensityPlot_bio_{N}.qs'` | `'density_plot_bio_{N}.qs'` | 1405 |
| `'PieMap_bio_{N}.qs'` | `'piemap_bio_{N}.qs'` | 1411 |
| `'PieMap_bio_{N}_{metric}.qs'` | `'piemap_bio_{N}_{metric}.qs'` | 1413 |
| `'PhenoMap_{trait}\\.qs$'` | `'phenomap_{trait}\\.qs$'` | 2178 |
| `'OverallImportance_{sfx}\\.qs$'` | `'overall_importance_{sfx}\\.qs$'` | 1933 |
| `'CumulativeImportance_{sfx}\\.qs$'` | `'cumulative_importance_{sfx}\\.qs$'` | 1938 |
| `'GeneticOffsetPieMap_{sfx}\\.qs$'` | `'genetic_offset_piemap_{sfx}\\.qs$'` | 1944 |
| `'GeneticOffsetPieMap_{sfx}_TajimaD\\.qs$'` | `'genetic_offset_piemap_{sfx}_tajima_d\\.qs$'` | 1949 |
| `'GeneticOffsetPieMap_{sfx}_PiDiversity\\.qs$'` | `'genetic_offset_piemap_{sfx}_pi_diversity\\.qs$'` | 1954 |
| `'GeneticOffsetPieMap_{sfx}_zoom_{zt}\\.qs$'` | `'genetic_offset_piemap_{sfx}_zoom_{zt}\\.qs$'` | 1986 |
| `'haplotype_{tag}/Region_{rid}_crosshap_viz\\.qs$'` | `'region_{rid}_crosshap_viz\\.qs$'` | 1887, 2555, 3136, 3427 |
| `'Region_{rid}_boxplot_{trait}\\.qs$'` | `'region_{rid}_boxplot_{trait}\\.qs$'` | 1903, 2570, 3151, 3445 |
| `'HaplotypePieMap_Region_{rid}_{trait}\\.qs$'` | `'haplotype_piemap_region_{rid}_{trait}\\.qs$'` | 1912, 2579, 3160, 3452 |
| `'haplotype_{tag}/clustree/Region_{rid}_clustree_MG\\.qs$'` | `'region_{rid}_clustree_MG\\.qs$'` | 3401 |
| `'haplotype_{tag}/clustree/Region_{rid}_clustree_hap\\.qs$'` | `'region_{rid}_clustree_hap\\.qs$'` | 3408 |

**Note**: `load_plot()` searches recursively, so if the search base is updated correctly, filename-only patterns (without directory prefixes) may be sufficient.

---

## 8. Haplotype Tables

| Old Path | New Path |
|---|---|
| `tables/haplotype_{tag}/Selected_regions.tsv` | `haplotype_scan/{tag}/tables/selected_regions.tsv` |
| `tables/haplotype_{tag}/scan_status.tsv` | `haplotype_scan/{tag}/tables/scan_status.tsv` |
| `tables/haplotype_{tag}/Region_{rid}_assignments.tsv` | `haplotype/{tag}/tables/region_{rid}_assignments.tsv` |
| `tables/haplotype_{tag}/Region_{rid}_frequencies.tsv` | `haplotype/{tag}/tables/region_{rid}_frequencies.tsv` |

---

## 9. Overlapping Mode Detection

Line 777:
```r
# Old
file.exists(paste0('/pipeline/', proj, '_results/tables/overlapping/Overlap_summary.tsv'))
# New
file.exists(paste0('/pipeline/', proj, '_results/overlapping/tables/overlap_summary.tsv'))
```

Line 2594:
```r
# Old
safe_fread(paste0('/pipeline/', local_proj, '_results/tables/overlapping/Overlap_summary.tsv'))
# New
safe_fread(paste0('/pipeline/', local_proj, '_results/overlapping/tables/overlap_summary.tsv'))
```

---

## 10. Migration Strategy

### Recommended Approach

1. **Replace `read_project_config()`** with `yaml::read_yaml()` to handle nested YAML
2. **Update all base path constants** in helper functions (Section 6)
3. **Update all filename patterns** (Sections 2, 3, 7)
4. **Update region ID regex** from `:` to `_` (Section 4)
5. **Update `load_plot()`** to search module-based directories instead of flat `plots/` + `intermediate/`
6. **Rename mode references**: `association_phenotypes` → `phenotype_association` throughout
7. **Test with SIMDATA** end-to-end

### Backward Compatibility Option

To support both old and new output structures during transition, `load_plot()` could search both old and new paths with a fallback. However, this adds complexity and is not recommended long-term.

---

## 11. Summary of Scale

- **17 helper functions** need path updates
- **~28 QS pattern references** need filename changes
- **~14 TSV file references** need filename + path changes
- **~6 `list.dirs()`/`list.files()` calls** need base path updates
- **3 config key lookups** need nested YAML parsing
- **1 mode rename**: `association_phenotypes` → `phenotype_association`
- **1 region ID format change**: `:` → `_`
