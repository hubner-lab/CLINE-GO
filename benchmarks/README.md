# CLINE-GO Benchmark Datasets

Two published datasets for validating the pipeline against known ground-truth loci.

## Datasets

| Dataset | Accessions | SNPs | Format | Pipeline Test | Climate |
|---------|-----------|------|--------|--------------|---------|
| Arabidopsis 1001G | 1,135 | ~50K (thinned) | VCF | Full (all modes) | yes |
| Populus balsamifera | 336 | ~107K | VCF (ready) | GEA + maladaptation | yes |

## Ground Truth Loci

### Arabidopsis 1001 Genomes (The 1001 Genomes Consortium, 2016)
| Gene | Chromosome | Position (Mb) | Trait | Reference |
|------|-----------|---------------|-------|-----------|
| FRI (FRIGIDA) | 4 | ~0.269 | Flowering time | Johanson et al. 2000 |
| FLC (FLOWERING LOCUS C) | 5 | ~3.17 | Flowering time | Michaels & Amasino 1999 |
| FT (FLOWERING LOCUS T) | 1 | ~24.33 | Flowering time | Kardailsky et al. 1999 |
| CRY2 (CRYPTOCHROME 2) | 1 | ~4.76 | Light-dependent flowering | El-Din El-Assal et al. 2001 |

### Populus balsamifera (Fitzpatrick et al. 2021)
No specific loci — validation is via common garden fitness correlation with GF genetic offset predictions.

## Prerequisites

- `bcftools` (for VCF thinning) — available via `nix shell nixpkgs#bcftools`
- `python3` (for metadata construction)
- `wget` (for data download)
- ~20 GB free disk space (Arabidopsis full VCF cached during thinning)

## Usage

### 1. Prepare datasets

```bash
# Each script downloads data, converts format, builds metadata
bash benchmarks/prepare_arabidopsis.sh  # ~18 GB download + thin to 50K SNPs
bash benchmarks/prepare_populus.sh      # VCF from GitHub, coords from Gougherty et al.

# Use --force to re-download/rebuild
bash benchmarks/prepare_arabidopsis.sh --force
```

### 2. Run pipeline

```bash
# Arabidopsis: full pipeline (GEA + GWAS + maladaptation)
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline cline-go:latest \
  snakemake -c10 -s Snakefile --config mode=processing --configfile config_arabidopsis.yaml --scheduler greedy

# Populus: GEA + maladaptation
docker run ... snakemake -c10 --config mode=processing --configfile config_populus.yaml ...
docker run ... snakemake -c10 --config mode=association --configfile config_populus.yaml ...
docker run ... snakemake -c10 --config mode=maladaptation --configfile config_populus.yaml ...
```

## Data Sources

| Dataset | Source | DOI/URL |
|---------|--------|---------|
| Arabidopsis VCF | 1001genomes.org | [v3.1 release](https://1001genomes.org/data/GMI-MPI/releases/v3.1/) |
| Arabidopsis coords + phenotypes | AraPheno | [REST API](https://arapheno.1001genomes.org/rest/) |
| Arabidopsis GFF | Ensembl Plants | Arabidopsis_thaliana.TAIR10 |
| Populus VCF | GitHub | [stephenrkeller/Fitzpatrick_etal_MER_2021](https://github.com/stephenrkeller/Fitzpatrick_etal_MER_2021) |
| Populus coords | GitHub | [agougher/poplarAdaptiveOffset](https://github.com/agougher/poplarAdaptiveOffset) (popInfo.csv) |
| Populus GFF | Ensembl Plants | Populus_trichocarpa.Pop_tri_v4 |

## Papers

- The 1001 Genomes Consortium. 2016. "1,135 Genomes Reveal the Global Pattern of Polymorphism in Arabidopsis thaliana." *Cell* 166(2):481-491.
- Fitzpatrick et al. 2021. "Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests." *Molecular Ecology Resources* 21(7):2271-2285.
- Gougherty et al. 2020. "Maladaptation, migration and extirpation fuel climate change risk in a forest tree species." *Nature Climate Change* 10:166-171.

## Notes

- Ensembl Plants GFF3 files do NOT include GO terms. Enrichment analysis will not produce results for benchmark datasets.
- The Arabidopsis VCF download is ~18 GB. The script caches it in `benchmarks/.cache/arabidopsis/` and thins to ~50K SNPs for pipeline use.
- Phenotype data (FT10, FT16) has ~1,100 values out of 1,135 accessions. The pipeline's `missing_strategy: "DROP"` handles samples without phenotype values.
- Populus uses P. trichocarpa reference genome for P. balsamifera (standard practice).
- Populus coordinates for 37/42 populations from Gougherty et al.; 5 missing populations (NIC, VER, USDA8, USDA19) use approximate coordinates.
- Downloads are cached in `benchmarks/.cache/` (gitignored).
- `climate: false` mode can be tested by running either dataset with `climate.enabled: false` in config.
