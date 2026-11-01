# Populus tutorial: sources

Every fact the tutorial data rests on, with where it was read and when. conf.py excludes this
directory from the page build; this file is a record, not a page.

Status 2026-10-31: **incomplete.** Phase 1 stopped at its hard stop and was set aside by the
user; Phase 2 went ahead on what is sourced. `sites.tsv` and `populus_metadata.tsv` hold 38 of
the 42 populations (325 of 336 trees); `gardens.tsv` and the predictor list are not written yet.
See "Open" at the end.

## Paper

- Fitzpatrick, M.C., Chhatre, V.E., Soolanayakanahally, R.Y. & Keller, S.R. (2021). Experimental
  support for genomic prediction of climate maladaptation using the machine learning approach
  Gradient Forests. Mol Ecol Resour 21(8):2749-2765. doi:10.1111/1755-0998.13374 (paywalled).
- Open preprint: Authorea, doi:10.22541/au.159863198.86187354, posted 2020-08-28 (Crossref record,
  read 2026-10-31). Crossref marks it `is-preprint-of` the MER article.
  - Not readable from the lab on 2026-10-31: `authorea.com/doi/full/...` and `.../doi/pdf/...`
    answer HTTP 403 (bot protection), to both curl and WebFetch.
  - The NSF Public Access Repository record `par.nsf.gov/servlets/purl/10210915` (found by web
    search, probably the accepted manuscript) refused the connection.
  - Not in Europe PMC as a preprint; OpenAlex lists the Authorea PDF as the only open copy.

## Authors' data repository

`github.com/stephenrkeller/Fitzpatrick_etal_MER_2021`, commit
`7b39b0034e9224da5a62ee3f63f03e197c5735df` (2021-03-03), cloned 2026-10-31. **No license:** read,
never copied into this repo; readers download it themselves.

| file | what it holds (README, checked against the file) |
|---|---|
| `balsam_core336inds_42pops.vcf.gz` | 57,482,278 bytes; 336 trees, 42 populations ("core ancestry deme") |
| `core336.gf.allSNPs.gardenOffset.csv` | GF offset per population on all SNPs; columns `pop_code, VT, IH`; 42 rows |
| `core336.gf.stdFreqs.gardenOffset.csv` | GF-X offsets (XTX-standardised frequencies); the paper's main offset; extra columns to the right are a pasted copy and residuals |
| `core336.gf.rawFreqs`, `.allOls`, `.lfmm`, `.bayenv`, `.Bayenv-lfmm` `.gardenOffset.csv` | offsets on each outlier set |
| `core336.500_RANDOM.gardenOffset.csv` | 999 replicates (`monteC` 1-999) x 42 populations = 41,958 rows; offsets on 500 random SNPs each. The authors' script averages them per population. |
| `core336.mahal.gardenOffset.csv` | Mahalanobis climate distance home -> garden per population, `pop_code, VT, IH`. **How it was computed is not in the repo** (needs the paper). |
| `core336.bf.pc{1,2}.emp.txt`, `manhattan_data_core_pc{1,2}.txt` | Bayenv2 and LFMM outliers on climate PC1 / PC2 |
| `snpR2.csv` | GF R^2 of modelled SNPs per outlier test |
| `gardenQuery.RData` | Vermont garden data (see below) |
| `IH_phenology_2015.csv` | Indian Head garden data (see below) |
| `GF_GBS_predict_gardens.R` | the authors' growth ~ offset analysis (last modified 2021-03-02) |

### Growth data and its model (from `GF_GBS_predict_gardens.R`, read 2026-10-31)

- Vermont (VT): `gardenQuery.RData` objects `HI_VT` (rep_code, pop_code, ind_code,
  height_increment) merged with `joiner` (row_num, tree_num, rep_code, garden_id). Zero
  `height_increment` is set to NA. The script counts 1,461 non-missing observations before the
  garden filter; restricted to `garden_id == "VT"`: 875 observations, 58 populations, **41 of the
  42 genotyped** (all but NIC). Other objects: `HI2_VT`, `DAM_VT` (winter damage, unused), and
  `drv`, a PostgreSQL driver object that does not load without RPostgreSQL (harmless; skip it).
- Indian Head (IH): `IH_phenology_2015.csv` (field, pop_code, ind, ind_code, rep, rep_100,
  rep_code, height_increment); fields `New_CG` 936, `Old_CG` 360, `P1` 45, `S2` 90 rows. P1 and S2
  are excluded: 592 non-missing observations, 19 populations, **9 of the 42 genotyped** (CLK CPL
  GAM KAP LON NBY OUT SKN TUR).
- Per-garden mixed models (lmerTest), BLUP = `coef(model)$pop_code` intercept:
  - VT: `height_increment ~ tree_num*row_num + (1|pop_code) + (1|ind_code)`
  - IH: `height_increment ~ field + (1|pop_code) + (1|ind_code)`
- Offset models, both gardens stacked: `lm(HtBlups ~ offset + I(offset^2) + Garden)` for each GF
  offset set; the climate-only model is linear, `lm(HtBlups ~ offset + Garden)` with the
  Mahalanobis distance as `offset`. Compared by `MuMIn::AICc`.
- So the chapter-7 fits run on at most 41 (VT) + 9 (IH) = 50 population-garden points.
- **Bug in the authors' script, recorded, not to be copied:** the IH "All outliers" merge takes
  `GF_predictAllOutliers[,c(1,2)]`, which is the VT column (the GF-X, random and climate merges
  take column 3 for IH). Any published AICc for the all-outliers model may carry it.

## Population coordinates

- `sites.tsv`: 38 of 42 populations, full precision, from Gougherty, A.V., Keller, S.R. &
  Fitzpatrick, M.C. (2021). Maladaptation, migration and extirpation fuel climate change risk in a
  forest tree species. Nat Clim Change 11:166-171, doi:10.1038/s41558-020-00968-6; data file
  `data/popInfo.csv` (81 populations; columns code, lat, long) in
  `github.com/agougher/poplarAdaptiveOffset`, commit `b4f82e432599607e40af7e9dda5a7a56fab2c90a`
  (2022-02-02), cloned 2026-10-31. Gougherty writes `USDA_N`; the VCF and the authors' tables
  write `USDAN`, so `USDA_1` -> `USDA1` and so on.
- The March metadata (`/mnt/data/eugene/ADAPTOGENE/data/Populus_balsamifera_metadata.tsv`) agrees
  with `popInfo.csv` for these 38 to its 3-decimal rounding.
- **Missing: 4 populations, not 5** (`benchmarks/README.md:81` says "5" and names these four): NIC,
  VER, USDA8, USDA19. Not in `popInfo.csv`. The March values for them (NIC 46.5/-81.0, VER
  44.48/-73.21, USDA8 48.17/-93.75, USDA19 45.42/-86.5) have no recorded source and are not used.
- Candidate source, not reachable from the lab: Dryad doi:10.5061/dryad.gp78p (Keller et al. 2018
  J Hered 109(1):47, 90 populations; its data statement says "sampling locations and SNP
  genotypes"), file `Keller_etal_2017_JoH_SNPdata_DRYAD.xlsx`. The Dryad API answers 401 without a
  bearer token.

## Data script (`get_populus_data.sh`, checked 2026-10-31)

- VCF: `raw.githubusercontent.com/stephenrkeller/Fitzpatrick_etal_MER_2021/7b39b00.../balsam_core336inds_42pops.vcf.gz`,
  sha256 `053c486247aa9ea78c9fb9d447c9e727273f35e3c81647b4aa68d99a960109f6` (same bytes as the
  clone; the repo uses no Git LFS). VCFv4.0, no `##contig` lines, FORMAT GT:AD:DP:GQ:PL (so the
  pipeline's depth QC runs), 336 samples named `<site>_<nn>`, 107,309 SNPs on chromosomes 1-19
  (ids `S<chr>_<pos>`). bcftools warns `PL should be declared as Number=G` and reads it.
- Kept gzipped: every rule that reads the input VCF (`workflow/rules/processing.smk`: plink
  `--vcf`, bcftools, vcftools with its `--gzvcf` branch) takes `.vcf.gz`. Checked in the image:
  `plink --vcf ...vcf.gz --const-fid --keep` loads 107,309 variants and keeps all 325 metadata
  trees with their underscore names intact (FID 0, IID `CHL_01`).
- Sorted copy, `balsam_core336inds_42pops_sorted.vcf.gz`, the file `populus.yaml` reads: the
  same header and the same 107,309 records, body sorted by chromosome then position
  (`LC_ALL=C sort -k1,1n -k2,2n -s`, `gzip -n`; the script checks the record count). The
  authors' file lists chromosomes 10-19 before 1-9, and `scripts/ld_decay_pairs.R` (mode
  structure) pairs PLINK's `.bim`, which PLINK writes in numeric order, with the VCF's rows by
  position: on the unsorted file the rule stops with `PLINK merged distinct VCF chromosome
  names` (run 2, 2026-10-31). Filed as a pipeline bug (lane Handoff); the user chose the sort
  as the tutorial's workaround.
- GFF: Ensembl Plants release 52,
  `ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-52/gff3/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v3.52.chr.gff3.gz`
  (the `.chr` file: chromosomes only), sha256
  `3152ef21bb160998efc1e31367a867ed28df29be730b239f27642a5e9e2046b5`, matching Ensembl's
  `CHECKSUMS` sum `20747 10486`. 39,514 `gene` rows; `##sequence-region` 1 is 50,495,391 bp (v3).
  Gene rows carry `ID=gene:POPTR_...v3`, `biotype`, `description`, `gene_id`; no `Name` and no
  GO terms. Unpacked, because `scripts/normalize_gff.py` opens the GFF as plain text.
- `populus_metadata.tsv`: one row per VCF sample, in VCF header order, site = the sample name
  before its last `_`, coordinates from `sites.tsv`. **325 trees of 38 sites**: the 11 trees of
  NIC (2), VER (3), USDA8 (3) and USDA19 (3) are left out until their coordinates have a source.
  The pipeline keeps only metadata samples (plink `--keep`), so those 11 trees are not analysed.
  When the four sites are sourced: add them to `sites.tsv` and rebuild this file the same way.

## Decisions made on the first run (Phase 3, 2026-10-31)

- **Predictors: provisional.** The paper's list is not sourced, so `populus.yaml` carries the
  schema default `bio_1,2,3,4,8,9,10,12,14,15` (the user's choice). When the paper's list is
  sourced, `climate`, `gea` and `maladaptation` are rerun; `processing`, `prestructure` and
  `structure` do not depend on it.
- **Clones removed, relatedness cutoff 0.975 (the user's choice).** Run 1
  (`populus-run1`, `relatedness_action: keep`) showed, among the 311 trees left after
  missingness filtering (`Processing/tables/relatedness_pairs.tsv`, plink IBS):
  - all-pair IBS median 0.943, 99th percentile 0.951; 48 pairs >= 0.99, 43 of them between
    trees of two different sites, mostly USDA sites;
  - pairs between 0.955 and 0.99: 10; above the highest same-or-cross-site non-clone pair
    (0.969, HST_09-HST_12) the next is 0.986, and three cross-site pairs sit at 0.989
    (SSR_12-USDA10_10, USDA13_01-USDA17_01, USDA18_01-USDA1_08), below the 0.99 default;
  - sNMF on these trees gives a 6-tree clone group (SSR_07, USDA1_04, USDA3_02, USDA7_07,
    USDA14_09, USDA14_10; IBS 0.992-0.995) a cluster of its own at every K from 2 to 8, and a
    4-tree group (USDA5_05, USDA9_07, USDA18_08, USDA18_09) one from K = 3: the clones, not
    ancestry, shape the cross-entropy curve.
  - The pipeline's rule (`relatedness_action: remove`) drops, pair by pair, the copy with more
    missing calls. A simulation of it on run 1's tables gives 281 trees at 0.99 (what the
    pipeline's own preview list says) and 278 at 0.975. USDA1 (4 trees) and USDA6 (1) lose all
    trees, USDA5 and USDA7 keep 1, USDA3, USDA9, USDA12, USDA14 and USDA21 keep 2.
  - SSR_07 survives as its group's copy (lowest missingness, 0.035): a genotype shared with
    four Great Lakes sites stays at a Wyoming site. Kept as the pipeline decides; filed as a
    gap (cross-site clones), not hand-filtered.
  - The authors' 336 trees include these clones. The chapter-7 check compares population-level
    offsets with the authors' garden BLUPs, which are not affected.

- **`sNMF.k_best: 2`**, by the rule "minimum mean cross-entropy" on
  `PreStructure/tables/cross_entropy_K2-8.tsv` (run 3, 278 trees, 10 runs per K): mean 0.69861
  at K = 2, rising at every K after (0.69982, 0.70254, ... 0.71314 at K = 8). The PCA agrees:
  one axis stands out (variance fraction 2.08%, then 0.77% and a plateau near 0.6%). Tracy-Widom
  p-values are not used: every axis sits at the 8e-09 floor. With the clones kept (run 1) the
  curve fell through K = 8 with no minimum.
- **`Map.resolution: 2.5`** arc-minutes: the present-day WorldClim 2.1 bioclim archive at 2.5'
  is 658 MB (`wc2.1_2.5m_bio.zip`, 697 MB unpacked); `download_climate_present` took 3 min on
  the lab (run 3, 2026-10-31). The finer 30'' set is about 10 GB (CLAUDE.md, SIMDATA's cache), too
  much for a tutorial; 2.5' (0.042 degrees) resolves the 36 sites, the closest two of which
  (USDA12, USDA13) lie 0.5 degrees apart.
- **RDA: `axes: auto`, `permutations: 99` (the user's choice; default 999).** At 999, mode=gea
  was killed twice at the documented `--memory=20g` and needed about 84 min for the by-axis
  `anova.cca` alone: measured in the image, the fit takes 3 s, by-axis 504 s per 99
  permutations. The permutations set only the ANOVA p-values (and through by-axis, K); the
  per-SNP p-values come from rdadapt. The constrained eigenvalues are flat (0.124 to 0.085 of
  constrained inertia over 10 axes) and adjusted R^2 after Condition(PC1 + PC2) is 0.005.
  The stopped 999 run peaked at 19.2 GB in its first 7 minutes; at 99 the documented
  `--memory=20g` with no BLAS flags is enough (run 4).
- **Runs.** Run 1: clones kept. Run 2: unsorted VCF, stopped in mode=structure (LD-decay bug).
  Run 3: sorted VCF, clones removed, empty WorldClim cache; gea and maladaptation rerun with
  64 GB and BLAS caps while chasing the RDA kill. **Run 4 is the recorded run:** a fresh folder
  and cache, all six modes in the reader's documented command, 2026-10-31; its tables and
  figures are byte-identical to run 3's.

## What run 3 showed (the plan's watch list)

- **Against March.** March (old code, 336 trees, clones kept) kept 320 trees and 11,699 SNPs.
  Run 3 keeps 278 trees and 11,819 SNPs: 11 trees are out because their four sites have no
  coordinate source yet, 14 fail the sample-missingness filter, 33 are clone copies. The SNP
  count is within 1%; the tree count differs by the clones and the four sites, not by a
  pipeline change.
- **EMMAX** ran on x86 through `scripts/emmax_run.sh`; the plink `--double-id` / `--const-fid`
  steps kept the underscore sample names (`CHL_01`) intact; lambda_GC 0.97-1.01 over the 10
  predictors.
- **Clumping distance.** `auto_genome_wide` resolved without stopping: the genome-wide LD-decay
  row is `reliable` (short-range r^2 0.448, critical r^2 of unlinked markers 0.0158), so the
  distance is its `decay_critical_bp`, 184,933 bp (half-decay would be 3,072 bp). Regions are
  therefore about 370 kb wide and hold 7-74 genes each (`verified/ch5_top_regions.tsv`).
- **GEA calls (Bonferroni):** LFMM 27 calls over 7 predictors, EMMAX 4, RDA 8 (0.01, K = 5
  axes chosen by the by-axis test): 26 SNPs. Only one SNP is called by two methods within the
  clumping distance (chr 2, EMMAX bio_10 + RDA), so the `agreement` SNP set has 1 SNP and
  mode=maladaptation skips it (below the 3-SNP floor); the offset uses `union` (26 SNPs).
  Chapter 7 compares `union` with random SNP sets and the climate distance (the user's
  choice, 2026-10-31).
- **Offset (GF, union, spatial, 2061-2080 SSP5-8.5, MPI-ESM1-2-HR):** site offsets 0.0083
  (USDA21) to 0.026 (CHL); 96.8% of future cells lie outside the training climate range
  (`gradient_forest_diagnostics.tsv`), so the map is mostly extrapolation (`extrap: true`).
- **Figures to judge:** `ch4_varpart_venn` (title and subtitle overprint) and
  `ch6_cumulative_importance` (axis titles and ticks collide) carry pipeline plot defects,
  filed in the lane Handoff; the combined Manhattan draws one threshold line (4.25e-6), not
  RDA's 8.5e-7.

## Open (blocks Phase 1)

- Coordinates of NIC, VER, USDA8, USDA19.
- Coordinates of the Vermont and Indian Head common gardens.
- The climate predictors of the GF models (and so the `bio_N` mapping).
- How the Mahalanobis climate distance was computed.

All four need the paper text (MER or the Authorea preprint), or for the coordinates the Dryad
gp78p sheet. The repos above do not hold them. Gougherty's `code/gfOffset.R` uses WorldClim 1.4
for a different paper and is not a stand-in for Fitzpatrick 2021's predictors.
