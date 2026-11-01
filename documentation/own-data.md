# Run on your own data

This page takes you from your own files to results in the app: the three inputs, a project
folder, a config written by `clinego init`, the settings to change first, and the modes in
order. It assumes the [install](install.md) is done and the [quick check](quickstart.md) passed.

## What you need

VCF
: Biallelic SNPs for every sample, `.vcf` or `.vcf.gz`. Filtering, LD pruning and imputation
  happen inside the pipeline, so give it the SNPs as called. A `chr` prefix on chromosome names
  is removed early, in the VCF and the GFF alike.

Metadata
: A tab-separated file, one row per sample. Columns 1 to 4 must be named `site`, `sample`,
  `latitude` and `longitude`, in that order; every further column is a measured trait.
  - `sample` matches the sample names in the VCF, once each and without spaces.
  - `site` groups samples into sampling sites (populations). Samples of one site should share
    its coordinates: each site is drawn at the mean of its samples' coordinates.
  - `latitude` and `longitude` are decimal degrees. A project without coordinates keeps both
    columns and leaves them empty (see `Regime.mode` below).
  - A trait name starts with a letter and holds only letters, digits and underscores
    (`height`, `flowering_time`, not `100_seed_weight`).

GFF3 (optional)
: A gene annotation for the same assembly as the VCF. With it, associated regions are annotated
  with the genes they hold. Its chromosome names must match the VCF's once the `chr` prefix is
  removed; a GFF that shares no chromosome name with the VCF stops the run.

The first lines of a metadata file:

```{code-block} text
site	sample	latitude	longitude	height	flowering_time
Negev	NEG01	30.854	34.7826	52.1	88
Negev	NEG02	30.854	34.7826	49.7	91
Galilee	GAL01	33.0128	35.4985	61.3	84
```

## Make a project folder

One folder per project holds the inputs, the config and every result. Put the input files in a
subfolder, here `data/`:

```{code-block} text
myproject/
├── data/
│   ├── mydata.vcf.gz
│   ├── mydata_metadata.tsv
│   └── mydata.gff3
└── config_myproject.yaml      written by clinego init
```

From inside `myproject/`, write the config:

```{code-block} text
|init_own|
```

`clinego init myproject` writes `config_myproject.yaml` from the pipeline's template, with
`project_name: myproject`, and refuses to overwrite a config that is already there. Results go
to `myproject_results/` and logs to `myproject_logs/`, both in the same folder. Name each
project after its dataset, and never give two datasets the same project name: the results
folder is keyed by the name alone.

## Set these first

Open `config_myproject.yaml` in a text editor. Every key is described in the
[configuration reference](manual/configuration.md); these are the ones a new project sets:

```{code-block} yaml
Input:
  dir: "data/"        # the input subfolder, with its trailing slash
  vcf: "mydata.vcf.gz"
  metadata: "mydata_metadata.tsv"
  gff: "mydata.gff3"  # or "" without a gene annotation

Regime:
  mode: standard      # standard, or gwas_only for a project without coordinates

Climate:
  # WorldClim bioclimatic variables, bio_1 to bio_19
  predictors: "bio_1,bio_4,bio_12,bio_15"

Map:
  resolution: 2.5     # arc-minutes; the template's 0.5 is a very large download

sNMF:
  k_start: 2          # the range of K (ancestral populations) to test
  k_end: 8
  k_best: null        # set after mode=prestructure, see below

GEA:
  configs:            # which association methods run, and their thresholds
    - method: "LFMM"
      adjust: "bonf"
      threshold: '0.05'
    - method: "EMMAX"
      adjust: "bonf"
      threshold: '0.05'
  snp_sets:           # the candidate SNP set the genomic offset is fitted on
    - name: default
```

- **`Regime.mode`**: `standard` runs everything and needs coordinates. `gwas_only` is for a
  project without coordinates or without climate: it runs `processing`, `prestructure`,
  `structure`, `traits`, `gwas` and `gea_x_gwas`, and refuses the climate modes.
- **`Climate.predictors`** are the variables every GEA method tests and the offset models are
  fitted on. Keep them few and not strongly correlated; `mode=climate` draws their correlations.
- **`GEA.configs`** lists the methods (`LFMM`, `EMMAX`, `RDA`, or a GAPIT model) with a
  significance rule each; `GWAS.configs` does the same for the traits. Keep both lists, even
  in a `gwas_only` project: `mode=gea_x_gwas` needs both to be there.
- **`GEA.snp_sets`** declares the candidate SNP sets for `maladaptation`; `mode=gea` builds
  them. Without one, `maladaptation` stops (see the
  [FAQ](faq.md#maladaptation-stops-with-no-curated-snp-sets)). Leave it out of a `gwas_only`
  project: there it names no predictors, and every mode stops when the run starts.

## Run the modes in order

Each command runs one mode and returns when it is done. Run them one after another, never two
at once for the same project.

### 1. Filter and check the data

```{code-block} text
|run_own_processing|
```

Look at `myproject_results/Processing/tables/filtering_summary.tsv`: how many samples and SNPs
each filter removed. If a filter removed too much, change `Filter.*` and run the mode again.

### 2. Population structure, then choose K

```{code-block} text
|run_own_prestructure|
```

This runs sNMF for every *K* from `sNMF.k_start` to `sNMF.k_end`. Every later mode except
`traits` needs `sNMF.k_best`, so stop here and choose it:

- `PreStructure/tables/cross_entropy_K<start>-<end>.tsv` holds the cross-entropy of every run
  at every *K* (columns `K`, `run`, `cross_entropy`, `best_run`). Lower is better; the usual
  choice is the *K* where the curve stops falling, which
  `PreStructure/plots/cross_entropy_K<start>-<end>.png` shows.
- `PreStructure/plots/K<k>/` has the ancestry bar plot and the PCA coloured by ancestry for
  each *K*.

Write the chosen number into the config as `sNMF.k_best`.

### 3. The rest, in this order

```{code-block} text
|run_own_structure|
|run_own_climate|
|run_own_traits|
|run_own_gea|
|run_own_gwas|
|run_own_gea_x_gwas|
|run_own_maladaptation|
```

- `structure` imputes missing genotypes, extracts the climate at each site (downloading
  WorldClim the first time) and measures LD decay, which sets the distance used to merge
  significant SNPs.
- `climate` characterizes the predictors and builds the spatial variables a spatial Gradient
  Forest needs, so run it before `maladaptation`.
- `traits` and `gwas` need trait columns in the metadata; skip them if it has none.
  `gea_x_gwas` needs both `gea` and `gwas` to have run.
- `maladaptation` computes the genomic offset on the SNP sets `gea` built.
- `pregea` (optional, between `climate` and `gea`) explores LFMM's *K* and the number of
  principal components before the full run.

[The modes](manual/modes.md) says what each mode needs and writes. Two traps:

- `maladaptation` with a spatial Gradient Forest (the template's `spatial_correction: with`)
  needs more sampling sites than `Climate.predictors` plus one: with ten predictors, at least
  twelve sites. With fewer sites, list fewer `Climate.predictors` (and run `climate` again), or
  set `Maladaptation.methods.gradient_forest.spatial_correction: without`.
- EMMAX does not run on a Mac with Apple silicon (see the
  [FAQ](faq.md#emmax-stops-at-once-on-my-mac)).

## Look at the results

Start the app in the project folder; it serves on <http://localhost:3838>:

```{code-block} text
|app|
```

Pick `myproject` in the project list. Each mode has its tab, and the GEA and GWAS tabs let you
change thresholds, combine methods and save SNP sets without running the pipeline again.
