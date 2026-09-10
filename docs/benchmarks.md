# CLINE-GO Benchmark Datasets

Benchmark datasets for validating GWAS, GEA, and maladaptation analysis. Each dataset targets a different sequencing technology, organism, and pipeline mode combination.

## Dataset Overview

| ID | Species | Technology | Samples | SNPs (full) | Chr subset | Pipeline modes | Primary validation |
|----|---------|-----------|---------|-------------|------------|----------------|-------------------|
| RICE | *Oryza sativa* (RDP1) | WGS/array 44K | 413 | 36,901 | chr1 ~3,700 or chr3 ~3,000 | Full (GWAS+GEA+maladaptation) | Known cloned genes |
| SORGHUM | *Sorghum bicolor* (SAP) | GBS 569K | 355 | 569,305 | chr9 ~57K or chr6 ~57K | No-climate (GWAS only) | Known cloned genes |
| POPULUS | *Populus balsamifera* | GBS | 336 | ~10-50K (TBD) | TBD after download | Full (GEA+maladaptation focus) | Common garden fitness |

## Dataset 1: Rice Diversity Panel 1 (RDP1)

### Source
- **Paper**: Zhao K et al. (2011) "Genome-wide association mapping reveals a rich genetic architecture of complex traits in *Oryza sativa*." *Nature Communications* 2:467. DOI: 10.1038/ncomms1467
- **Companion**: McCouch SR et al. (2016) *Nature Communications* 7:10532. DOI: 10.1038/ncomms10532
- **Data**: http://www.ricediversity.org/data/sets/44kgwas/
- **GFF**: Rice MSU7 from rice.uga.edu or Ensembl Plants (IRGSP-1.0)

### Files
| File | Format | Description |
|------|--------|-------------|
| `RiceDiversity.44K.MSU6.Genotypes_PLINK.zip` | PLINK ped/map | 44K SNP genotypes |
| `RiceDiversity_44K_Phenotypes_34traits_PLINK.txt` | PLINK pheno | 34 traits |
| `RiceDiversity.44K.germplasm.csv` | CSV | Passport: lat/lon, country, subpopulation |

### Preparation
```bash
# Convert PLINK to VCF
plink --file RiceDiversity.44K.MSU6 --recode vcf --out rice_rdp1
# Subset to chromosome 1 (contains sd1)
plink --file RiceDiversity.44K.MSU6 --chr 1 --recode vcf --out rice_rdp1_chr1
```

### Phenotypes for benchmarking
| Trait | Expected locus | Chromosome | Gene | Effect |
|-------|---------------|------------|------|--------|
| Plant height | sd1 | chr1 | *SD1* (GA20ox2) | Semi-dwarfism (Green Revolution) |
| Grain length | GS3 | chr3 | *GS3* | Major grain size QTL |
| Amylose content | Waxy | chr6 | *Wx* | Starch synthesis |
| Pericarp color | Rc | chr7 | *Rc* | Proanthocyanidin regulation |

### Ground truth regions
| Locus | Chr | Approximate position | Notes |
|-------|-----|---------------------|-------|
| sd1 | 1 | ~38.3 Mb (MSU7) | Very strong signal in height GWAS |
| GS3 | 3 | ~16.7 Mb | Very strong signal in grain length GWAS |
| Waxy | 6 | ~1.8 Mb | Strong signal in amylose GWAS |
| Rc | 7 | ~6.1 Mb | Strong binary trait signal |
| Hd1 | 6 | ~9.3 Mb | Flowering time |
| Pi-ta | 12 | ~10.6 Mb | Blast resistance |

### Expected pipeline results
- **GWAS**: Manhattan peaks at sd1 (chr1), GS3 (chr3) for height/grain length
- **GEA**: Climate-associated SNPs expected given 82-country, 5-subpopulation distribution
- **Maladaptation**: Gradient Forest should capture latitudinal/temperature gradients across indica/japonica split

---

## Dataset 2: Sorghum Association Panel (SAP)

### Source
- **Paper**: Morris GP et al. (2013) "Population genomic and genome-wide association studies of agroclimatic traits in sorghum." *PNAS* 110(2):453-458. DOI: 10.1073/pnas.1215985110
- **Updated SNPs**: Miao C et al. (2019) Figshare. 569K imputed GBS SNPs.
- **Data**: https://figshare.com/articles/dataset/Untitled_Item/11462469/5
- **GFF**: Sorghum bicolor v3 from Ensembl Plants or Phytozome

### Files
| File | Format | Description |
|------|--------|-------------|
| HapMap genotype files | HapMap (TASSEL) | Per-chromosome, 569K SNPs |
| Phenotype data | TSV | Height components, flowering time (Morris 2013 supplementary) |

### Preparation
```bash
# Convert HapMap to VCF using TASSEL
run_pipeline.pl -Xmx4g -fork1 -h sorghum_chr9.hmp.txt \
    -export sorghum_sap_chr9 -exportType VCF
# Or use custom conversion script
```

### Phenotypes for benchmarking
| Trait | Expected locus | Chromosome | Gene | Effect |
|-------|---------------|------------|------|--------|
| Plant height | Dw1 | chr9 | Sobic.009G229800 | Major dwarfing gene |
| Plant height | Dw2 | chr6 | Sobic.006G067700 | Major dwarfing gene |
| Plant height | Dw3 | chr7 | Sobic.007G163800 | Auxin efflux transporter |
| Flowering time | Ma6 | chr6 | Sobic.006G004400 | Ghd7 ortholog |

### Ground truth regions
| Locus | Chr | Approximate position | Notes |
|-------|-----|---------------------|-------|
| Dw1 | 9 | ~57.4 Mb (v3) | Strong height signal |
| Dw2 | 6 | ~42.9 Mb (v3) | Strong height signal |
| Dw3 | 7 | ~58.6 Mb (v3) | Moderate height signal |
| Ma1 | 6 | ~40.3 Mb | Maturity/flowering |
| Ma6 | 6 | ~3.5 Mb | Maturity/flowering |

### Expected pipeline results
- **GWAS**: Manhattan peaks at Dw1 (chr9) for height, Ma6 (chr6) for flowering time
- **No-climate mode**: This dataset tests GWAS-only workflow (limited geographic coordinates)

---

## Dataset 3: Populus balsamifera

### Source
- **Paper**: Fitzpatrick MC, Chhatre VE, Soolanayakanahally RY, Keller SR (2021) "Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests." *Molecular Ecology Resources*. DOI: 10.1111/1755-0998.13374
- **Data**: https://github.com/fitzLab-AL/geneticOffsetR (VCF + R code)
- **Companion**: https://github.com/stephenrkeller/Fitzpatrick_etal_MER_2021 (data files)

### Files
| File | Format | Description |
|------|--------|-------------|
| `balsam_core336inds_42pops.vcf.gz` | VCF (bgzipped) | GBS genotypes, 336 individuals |
| CSV files | CSV | Genetic offset predictions, outlier SNP lists |
| Common garden data | CSV | Height growth in 2 gardens |

### Preparation
```bash
# VCF is ready to use — subset to chromosome if needed after inspection
bcftools view -r chr1 balsam_core336inds_42pops.vcf.gz -o populus_chr1.vcf
# Inspect SNP count and chromosome naming
bcftools stats balsam_core336inds_42pops.vcf.gz | grep "number of SNPs"
```

### Ground truth
| Validation type | Source | Description |
|----------------|--------|-------------|
| Genetic offset | Common garden | Populations with larger predicted offset performed worse (height growth) |
| GEA outliers | Bayenv2 | SNP lists provided in supplementary |
| GEA outliers | LFMM | SNP lists provided in supplementary |
| Climate adaptation | WorldClim | Range-wide across North America, 42 populations |

### Expected pipeline results
- **GEA**: LFMM should identify climate-associated SNPs; compare overlap with paper's Bayenv2/LFMM outliers
- **Maladaptation**: Gradient Forest genetic offset map should show latitudinal gradient; compare pattern with paper's Fig. 3
- **Validation**: Rank correlation between our genetic offset per population and their common garden fitness measurements

---

## Benchmark Comparison Plan

### For each dataset:
1. Run CLINE-GO pipeline (appropriate modes)
2. Extract significant regions from `regions_combined.tsv`
3. Check if known ground truth loci fall within identified regions
4. For maladaptation: visually compare genetic offset spatial patterns

### Success criteria
| Metric | Method |
|--------|--------|
| GWAS locus recovery | Does a significant region overlap the known gene? (within region_distance) |
| GEA locus overlap | Jaccard or % overlap with published outlier SNP lists |
| Genetic offset pattern | Visual comparison of piemap spatial gradients |
| Population structure | K matches expected subpopulation count |

### Results table template
| Dataset | Trait | Expected locus | Recovered? | Region | Distance to gene | Rank (by SNP count) |
|---------|-------|---------------|------------|--------|-------------------|---------------------|
| RICE | height | sd1 (chr1:38.3Mb) | | | | |
| RICE | grain_length | GS3 (chr3:16.7Mb) | | | | |
| SORGHUM | height | Dw1 (chr9:57.4Mb) | | | | |
| POPULUS | offset | latitudinal gradient | | | | |

---

## No-Climate Mode Plan

Config: `climate: false` (or `climate.enabled: false`)

### What changes:
- **Disabled**: climate download, LFMM, density plots, correlation heatmap, piemaps (ancestry/climate), maladaptation (GF), GEA
- **Enabled**: processing, structure (PCA, sNMF, cross-entropy), structure_K (pop stats without climate plots), GWAS (EMMAX on phenotypes), regions, genes, enrichment, manhattan plots, haplotype scan/analysis

### Implementation approach:
- No separate rules or scripts
- `get_targets()` returns fewer targets when `climate: false`
- Rules that need climate inputs are simply not requested
- Config parsing sets `CLIMATE_ENABLED = config.get('climate', {}).get('enabled', True)` or similar
