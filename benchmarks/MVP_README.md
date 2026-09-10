# MVP benchmark — 14 converted replicates

Simulated GEA benchmark replacing the retired Láruson set (`_archive/laruson/`).
Source: **Lotterhos 2023**, PNAS 120(12):e2220313120, doi:10.1073/pnas.2220313120.
Data: **BCO-DMO doi:10.26008/1912/bco-dmo.889769.1**.
Selection rationale and the paper argument: `docs/gea_simulation_benchmarks.md` §8.2–8.4 — do not re-derive it.

## What is here

| Path | What |
|---|---|
| `benchmarks/mvp_select_seeds.R` | Builds the 14-seed manifest from the deposit's own summary table |
| `benchmarks/mvp_seeds.tsv` | The manifest: 12 primary + 2 degenerate controls, with each seed's published statistics |
| `benchmarks/fetch_mvp.sh` | Streams the three tars, extracting only the wanted seeds' members |
| `benchmarks/convert_mvp.R` | One replicate → VCF + metadata + env table + 3 truth tables |
| `benchmarks/convert_mvp_all.sh` | Runs the converter over the manifest, rolls up provenance |
| `benchmarks/mvp_write_configs.R` | Emits the 14 `config_MVP{seed}.yaml` |
| `data/mvp/selection/` | The two selection inputs, persisted (summary CSV + params file) |
| `data/mvp/raw/` | 60 MB of extracted per-seed members |
| `data/mvp/MVP{seed}/` | Pipeline-ready inputs + truth tables, one dir per replicate |
| `data/mvp/conversion_summary.tsv` | Per-seed provenance roll-up |

Regenerate everything from scratch:
```bash
benchmarks/fetch_mvp.sh          # ~40 min, streams 23.5 GB, keeps 60 MB
benchmarks/convert_mvp_all.sh 7  # ~2 min
docker … Rscript /pipeline/benchmarks/mvp_write_configs.R
```

## The scoring contract — the one thing that is easy to get wrong

`causal_temp` and `causal_sal` are **independent** per-axis classifications: a locus can be
causal for temperature and neutral for Env2. `lib_detection.R:load_truth()` takes a single
`category` column and `call_by_threshold()` unions calls across traits, so pairing the wrong
truth table with the wrong predictor silently converts real detections into false positives.

`load_pvalues()` treats **every** column that is not `SNPID/chr/pos/key/n_snps/mean_maf` as a
trait, and the trait column names differ by method family — so the correct pairing is not the
same for all three:

**Univariate GEA (LFMM, EMMAX)** — one column per predictor, named `bio_1` / `bio_2`:
```
--truth=truth_temp.tsv  --traits=bio_1     # confounded axis,    R²(PC1~temp) 0.31–0.57
--truth=truth_sal.tsv   --traits=bio_2     # orthogonal control, R²(PC1~Env2) ≈ 0.00
--truth=truth_any.tsv   --traits=all       # union
```

**RDA** — multivariate over both predictors at once. `rda.R:689` writes a **single**
`climate_multivariate` column, so there is no per-axis p-value and `--traits=bio_1` would match
nothing. RDA has exactly one valid pairing:
```
--truth=truth_any.tsv   --traits=all       # the ONLY correct pairing for RDA
```

**GWAS** — trait columns are the metadata phenotype names, not `bio_*`:
```
--truth=truth_temp.tsv  --traits=phen_temp
--truth=truth_sal.tsv   --traits=phen_sal
--truth=truth_any.tsv   --traits=fitness   # fitness ONLY against truth_any (see below)
```
`fitness` is home-deme fitness, which is polygenic across **both** environmental axes — scoring
it against a single-axis truth table is a category error, so it pairs with `truth_any` only.

Never mix. Any other combination is a scoring error, not a result. Before scoring a table you
have not scored before, run `head -1` on it and confirm its trait columns.

**Headline precision is `precision_strict`.** Only `background_neutral` loci (linkage groups
11–20) count as false positives — which is the source paper's own convention: *"the performance
statistics were calculated by only including truly neutral loci unaffected by selection on
linkage groups 11 to 20 and the QTNs"* (Lotterhos 2023, Methods). Report `precision_all`
alongside for transparency; roughly half the genome is `linked_neutral`, so the two differ a lot.

For the four 1-trait seeds, `truth_sal.tsv` has **zero** causal loci. That is not a defect — it
makes bio_2 a pure false-positive control on those replicates.

## Design decisions worth defending

**`Filter.maf: 0.01`, not the 0.05 real-data default.** The deposit's genotype matrix and its
mutations/truth table are both already filtered at MAF > 0.01 in the same 1000 individuals, so
0.01 is the no-op that keeps every causal locus testable. Verified on seed 1231288: 11,194 →
11,194 SNPs, zero loss. This is exactly the defect that disqualified Láruson (truth defined at
one MAF, VCF at another → recall capped at ~27% by construction). `loci_freq.tsv` carries
`a_freq_subset`, so any other denominator is a column filter rather than a re-run.

**`sNMF.k_best` comes from the deposit's per-seed `K` column**, floored at 2. sNMF cross-entropy
has no interior minimum on this landscape — 100 demes on a continuous stepping-stone grid, so
there is no discrete K to find. Measured on seed 1231288, cross-entropy falls monotonically
across the whole tested range (K=2 → 10: 0.5416, 0.5204, 0.5061, 0.4987, 0.4943, 0.4889, 0.4848,
0.4807, 0.4770), so "pick the minimum" would just return `k_end`. Using the authors' published K
follows the same principle the seed selection rests on. `k_best` is in any case only the centre
of the PreGEA sweeps; the benchmark exists to find the best setting by scoring against truth.

**Genome-scale rescalings.** The simulated genome is 20 linkage groups × 50,000 sites = 1 Mb, so
several pipeline defaults sized for real genomes were reduced: `LD.window` 100 → 10 kb (the
default spans an entire LG), `snp_clumping_distance` 100,000 → 5,000 (the default collapses a
whole LG into one region), `LDdecay.max_distance` 500 → 50 kb.

**No `future_table`, no `mode=maladaptation`.** MVP ships home-deme fitness only, never transplant
fitness, so there is no honest offset target. A synthetic "future" shift would be a fabricated
target. The offset arm belongs to Capblancq et al. 2026 (Dryad doi:10.5061/dryad.jdfn2z3m6).

**No GFF, no WZA.** Simulated genome — annotation and GO enrichment are inert. ~10k SNPs over
20 LGs gives order 10¹–10² windows, an order of magnitude short of what window statistics need.

## Coordinate conventions (established empirically from the data)

- `pos_pyslim` is **1-based** in the 1 Mb genome; linkage group *k* spans
  `(k-1)*50000+1 .. k*50000`, and a locus exactly on a boundary belongs to the **lower** group.
  Confirmed by the deposit: `5-250000` has LG 5, `13-650000` has LG 13.
- VCF `CHROM` = LG (1–20), `POS` = `pos_pyslim - (LG-1)*50000` ∈ 1..50000, `ID` = `mutname`.
- `mutname` is `<LG>-<pos>` with an optional `_<n>` suffix marking a recurrent mutation at the
  same map position. Those are the duplicate `(chr,pos)` rows the converter collapses, keeping
  the causal row; the counts match one-for-one (e.g. seed 1231288: 93 suffixed names, 93 collapsed).
- Demes sit on a 10×10 grid mapped to lon 34.05–35.85 / lat 31.05–32.85, 0.2° apart, offset by
  half a raster cell so they land at cell **centres** — `stage_custom_climate.R:160` errors if two
  demes resolve to the same cell. The box is arbitrary; env values come from `temp_opt`/`sal_opt`
  and never from a raster lookup.

## Alignment guarantees

The converter does not trust row order. Gmat row names are `mutname` and column names are
`indID`, so both joins are by key. On top of that it re-derives each locus's allele frequency
from the genotype matrix and compares it to the deposit's `a_freq_subset`:

**All 14 replicates converted with `af_max_deviation = 0`** — exact agreement, so row *i* of the
matrix provably is locus *i* of the mutations table.

Independently verified against the deposit's own summary table: SNP counts 14/14, causal counts
14/14, 100 demes × 1000 individuals 14/14.

## Verified pipeline behaviour (seed 1231288)

- 1000/1000 samples retained, zero coordinate drops, zero missing genotypes
  (`vcf2lfmm.R`'s zero-missingness guard handles the SLiM case).
- LD pruning 11,194 → 3,648 SNPs.
- **PC1 = 5.85 %** of genetic variance (Láruson: 0.38 %) — structure correction is load-bearing here.
- **R²(PC1 ~ bio_1) = 0.495**, inside the 0.30–0.60 partial-confounding target band.
- **R²(PC1 ~ bio_2) = 0.006** — Env2 is a genuine orthogonal control.

Those last three are the whole reason this dataset replaces Láruson, measured through our own
pipeline rather than taken from the paper.

## Pilot baseline (seed 1232548, `k_best=5`) — the full chain is proven

`processing` → `prestructure` → `structure` → `gea` all exit 0, and the scoring harness runs
against the result. Truth-join gate: **10,325 / 10,325 p-value keys intersect the truth table**,
zero drift either direction, all 30 temperature-axis causal loci present.

`benchmarks/eval_detection.R --mode=rank`, 10,325 ranked SNPs:

| Method | Truth / trait | testable causal | AUC-PR | best F1 |
|---|---|---:|---:|---:|
| LFMM | `truth_temp` / `bio_1` (confounded) | 30 | **0.222** | 0.300 |
| LFMM | `truth_sal` / `bio_2` (orthogonal control) | 40 | **0.396** | 0.492 |
| EMMAX | `truth_temp` / `bio_1` (confounded) | 30 | **0.031** | 0.091 |
| EMMAX | `truth_sal` / `bio_2` (orthogonal control) | 40 | **0.104** | 0.204 |
| RDA | `truth_any` / all (multivariate) | 70 | 0.222 | 0.302 |

Two things to read off this, both of which are the benchmark working as designed:

1. **The orthogonal control axis outscores the confounded axis by ~1.8× for both univariate
   methods** (LFMM 0.396 vs 0.222; EMMAX 0.104 vs 0.031). Same genotypes, same run, same
   number of causal loci in the same order of magnitude — the only difference is how much
   population structure aligns with the predictor. That within-dataset positive/negative
   control is precisely what the seed selection was for, and no other candidate dataset in
   `docs/gea_simulation_benchmarks.md` offers it.
2. **EMMAX at 5 PCs collapses on the confounded axis** (0.031 vs LFMM's 0.222). This is the
   structure-correction cost/benefit tradeoff the benchmark exists to measure — and the thing
   Láruson structurally could not exercise. Do not read it as "EMMAX is bad": `n_pcs` was not
   tuned here, and finding the right rung is the job of the `mode=pregea` ladders plus
   `benchmarks/score_ladders.sh`. These numbers are an untuned baseline, not a result.

## Pipeline bug found and fixed during this pilot

`scripts/find_genes_around_regions.R` crashed whenever regions existed but no GFF was supplied.
`Input.gff` is optional and `normalize_gff` writes a **zero-byte** file when it is unset, but
`read_gff()` runs `grep -v '#'` on it, grep exits 1 for "no lines matched", and `fread` surfaces
that as a misleading *"disk is full in the temporary directory"* error. Every pre-existing
project ships a GFF, so this path had never run. Fixed by adding an empty-GFF guard alongside
the existing empty-regions guard, emitting the same empty gene tables. Gene annotation is
correctly inert on all 14 MVP replicates.

**`MVP1231288_results/` is stale.** It was built before the linkage-group-boundary fix, so its
`_work/` tree carries POS one higher than the current VCF. The genotypes are byte-identical, so
the three numbers above stand unchanged — but rerun it from `mode=processing` before scoring
anything from it. Snakemake will rebuild automatically (the input mtime changed); this note
exists only so the directory is not mistaken for current output.

## Published baseline — compare every future number against the original study

**`benchmarks/mvp_published_baseline.tsv`** — one row per benchmark seed, 36 columns, holding
Lotterhos 2023's *own* published detection metrics for LFMM, RDA, Kendall correlation and
phenotype GWAS. Regenerate with `benchmarks/mvp_published_baseline.R`. This is the fast-compare
path: a join on `seed`, not a re-derivation.

```r
pub <- fread("benchmarks/mvp_published_baseline.tsv")
pub[seed == 1232548, .(pub_lfmm_aucpr_sal_neut, pub_rda_aucpr_neut, pub_lfmm_tpr_temp)]
```

### Which column to compare against

**Always `*_neutSNPs`, never `*_allSNPs`.** Their `*_neutSNPs` is defined as *"including only
causal loci and neutral loci not affected by selection (any non-causal loci that arose on the
half of the genome affected by selection was excluded)"* — identical to our convention, where
`lib_detection.R` counts only `background_neutral` as false positives and ranks but skips
`linked_neutral` in the precision denominator. Verified numerically equal on seed 1232548:
deleting linked-neutral rows entirely and skipping them in the denominator both give AP =
0.2027 (temp) / 0.3938 (Env2).

| Our output | Compare against |
|---|---|
| LFMM AUC-PR, `truth_temp` + `bio_1` | `pub_lfmm_aucpr_temp_neut` **(see defect below)** |
| LFMM AUC-PR, `truth_sal` + `bio_2` | `pub_lfmm_aucpr_sal_neut` |
| RDA AUC-PR, `truth_any` + all | `pub_rda_aucpr_neut_corr` — **but see the RDA parameter section; neither published arm is an exact match** |
| LFMM TPR/FDR at `--adjust=qval --value=0.05` | `pub_lfmm_tpr_*` / `pub_lfmm_fdr_*_neut` |
| GWAS TPR/FDR (`phen_temp`/`phen_sal`) | `pub_gwas_tpr_*` / `pub_gwas_fdr_*` |
| EMMAX (any axis) | **no counterpart** — they ran no univariate mixed model on the environment |

### Known defect in the published table

`pub_lfmm_aucpr_temp_neut` is **not reproducible** from the deposit's own temperature p-value
column. Re-scoring `LEA3.2_lfmm2_mlog10P_tempenv` gives **0.2027** for seed 1232548 against a
published **0.0783**; the salinity column gives **0.0791**, i.e. the published "temp" figure
matches the *salinity* p-values within the same ±0.001–0.007 offset every other published cell
shows. Consistent with a column swap in that one metric. Not confirmed against the authors'
`MVP-NonClinalAF` source — **verify before citing**, and for comparison use the re-scored value
from `benchmarks/mvp_score_published.R` rather than this cell. Flagged in the table as
`lfmm_temp_aucpr_suspect`.

### Re-scoring their p-values (the cleanest comparison)

`benchmarks/mvp_score_published.R --seed=N` rebuilds the authors' per-SNP p-values
(`LEA3.2_lfmm2_mlog10P_*`, `RDA_mlog10P`, `RDA_mlog10P_corr`) as harness-format wide tables
under `benchmarks/mvp_eval/published/MVP{seed}/`, using the *same* coordinate derivation and
duplicate collapse as `convert_mvp.R` so the keys join to the same truth tables. Scoring those
through `eval_detection.R` removes every scoring-convention difference, so any remaining gap is
a difference in how the association was run — not in how detection was counted.

Validation on seed 1232548 (their p-values → our harness vs their published value):

| Metric | Published | Re-scored | Δ |
|---|---:|---:|---:|
| LFMM Env2 AUC-PR | 0.3924 | 0.3939 | +0.0015 ✓ |
| RDA AUC-PR (no corr.) | 0.1373 | 0.1396 | +0.0023 ✓ |
| RDA AUC-PR (struct. corr.) | 0.0738 | 0.0805 | +0.0067 ✓ |
| LFMM temp AUC-PR | 0.0783 | 0.2027 | +0.124 ✗ (the defect above) |

### Methodological differences between our run and theirs

These are the deltas that could legitimately move a number. State them whenever a comparison
is reported.

| # | Aspect | Lotterhos 2023 | CLINE-GO | Effect on comparability |
|---|---|---|---|---|
| 1 | SNP set | 10,394 loci (seed 1232548) | 10,325 | We collapse recurrent mutations sharing a map position (69 here), keeping the causal row; they keep both. <0.7% of loci |
| 2 | MAF | MAF > 0.01, applied by them before deposit | `Filter.maf: 0.01` — a verified no-op | **None** — identical effective threshold |
| 3 | LFMM training set | latent factors on the full SNP set | latent factors on the **LD-pruned** set (4,066 of 10,325), tested on the full set | Different structure model. Plausible cause of our +10% AUC-PR |
| 4 | LFMM K | their `K` column (5 for this seed) | `sNMF.k_best`, sourced from that same column | **None by construction** |
| 5 | RDA structure correction | reports both arms: none, and `Condition(PC1+PC2)` = **2** PCs | `Condition(PC1..PC5)` = **5** PCs (`condition_pcs` defaults to `@k_best`) | see the RDA section below — neither published arm matches ours exactly |
| 6 | RDA thresholding | q < 0.05 → TPR 0.214 at FDR 0.857 | 0 calls at q < 0.05 | **explained** — our `pmax` intersection rule, see below. Not a calibration defect |
| 7 | EMMAX on environment | not run | run | No baseline exists; our EMMAX numbers are novel, not a reproduction |
| 8 | Significance threshold | q < 0.05 throughout | config default is `bonf 0.05` | For comparison always pass `--adjust=qval --value=0.05` |
| 9 | Truth categories | `causal` / `neutral-linked` / `neutral` | same, renamed to `causal` / `linked_neutral` / `background_neutral` | **None** — counts verified identical (30 / 40 / 70) |
| 10 | Coordinates | genome-global `pos_pyslim` | per-linkage-group `POS` | **None** — both our VCF and our truth use the same derivation |

### RDA: exact parameter comparison, and why our scan calls nothing

Their code (`ModelValidationProgram/MVP-NonClinalAF`, `src/c-AnalyzeSimOutput.R`) versus ours
(`scripts/rda.R` + `scripts/R/lib/rdadapt.R`). The `rdadapt` function itself is **identical** in
both — `covRob(..., estim="pairwiseGK")` → `lambda <- median(resmaha)/qchisq(0.5, df=K)` →
`pchisq(resmaha/lambda, K, lower.tail=FALSE)`. Everything that differs is around it:

| Parameter | Lotterhos 2023 | CLINE-GO | Consequence |
|---|---|---|---|
| Ordination call | `rda(t(G) ~ sal + temp)` — no `scale=` → **`scale=FALSE`** (covariance-based) | `rda(..., scale=TRUE)` (correlation-based) | SNP columns standardised to unit variance → different loadings → different Mahalanobis distances |
| Structure correction | two separate arms: none, and `Condition(PC1 + PC2)` — **2 PCs** | `Condition(PC1..PC5)` — **5 PCs**, because `condition_pcs` defaults to the `@k_best` sentinel = `sNMF.k_best` | we remove more structure-aligned variance; closest published arm is `_corr`, but not equal |
| Axes K for rdadapt | **fixed `rdadapt(rdaout, 2)`** | `axes: auto` → axes with `anova.cca(by="axis")` p < 0.05, capped at rank, floored at 2 → **K = 2 here** | **same value**, arrived at differently. Cannot diverge with 2 predictors (rank ≤ 2) |
| p-value combination | **none** — each fit reported separately | **`pmax(p_partial, p_unconstrained)`** (`rda.R:478`, design note B6) — a SNP survives only if *both* the corrected and uncorrected fits call it | the whole thresholding difference; see below |
| Significance rule | `qvalue(...)$qvalues < 0.05` | `bonf 0.05` in these configs (`gea.py` default is `bonf 0.01`, per Capblancq & Forester 2021) | ours far stricter even before `pmax` |
| Marker set for the fit | full set | full set (`fit_mode: full`, 10,325 markers) | none |

**Why our scan calls 0 SNPs — it is the `pmax`, not miscalibration.** Measured p-value
quantiles on seed 1232548:

| | q=0.001 | q=0.01 | q=0.05 | median | λ_GC(df=1) |
|---|---:|---:|---:|---:|---:|
| ours (`pmax` of 2 fits) | 0.0092 | 0.0750 | 0.2169 | **0.7060** | 0.313 |
| theirs, no correction | 0.0000 | 0.0000 | 0.0031 | 0.4998 | 1.001 |
| theirs, `Condition(PC1+PC2)` | 0.0000 | 0.0000 | 0.0036 | 0.5004 | 0.998 |
| uniform null | 0.0010 | 0.0100 | 0.0500 | 0.5000 | 1.000 |

For two ~uniform p-values the median of their maximum is `0.5^(1/2) = 0.7071`. We observe
**0.7060**. That is the `pmax` rule reproducing itself exactly — each fit is individually
well-calibrated (the GIF step forces median p = 0.5 by construction), and taking the maximum
squares the tail, so our smallest p across 10,325 SNPs is 3.7e-4 and nothing survives any
multiple-testing rule. **λ_GC = 0.313 is a property of the intersection rule, not evidence of a
broken test.**

**Superseded in part by the sweep below.** The `pmax` arithmetic above is correct, but `pmax` is
not the *cause* of the zero calls: cell C (`scale=TRUE`, 2 PCs, **no** `pmax`) already calls 0,
while cell B — identical but `scale=FALSE` — calls 408. `scale=TRUE` is the driver; `pmax`
compounds it. And because `pmax` combines *two* p-value vectors it is **not** a monotone
transform of one, so it does change the ranking — measurably for the better (+139% AUC-PR).

Consequence for interpretation: **our RDA is a ranking method on this dataset, not a thresholding
one.** AUC-PR is 0.222 versus their 0.137 uncorrected / 0.074 corrected. Any comparison of RDA
*call counts* against theirs is meaningless unless the `pmax` rule is disabled or `condition_pcs`
is set to match. To reproduce their arms directly:

- their uncorrected arm → `condition_pcs: 0` (this also skips the second fit, so `pmax` is a no-op)
- their corrected arm → `condition_pcs: 2`, and note the `pmax` and `scale=TRUE` deltas remain

### RDA parameter sweep — which settings actually detect better (one seed, NOT yet applied)

`benchmarks/rda_param_sweep.R --seed=N` isolates the three deltas by re-fitting the same data
with each combination. It is a **benchmark probe, not pipeline code**: it replicates `rda.R`'s
data loading and re-uses the pipeline's own `rdadapt()` and `load_pca_covariates()` verbatim, so
only the three parameters vary. Writes nothing into any `{PROJECT}_results/` tree. Results:
`benchmarks/mvp_eval/rda_sweep/{PROJECT}/sweep_scored.tsv` plus per-cell p-value tables.

Seed 1232548, truth `truth_any` (70 causal / 10,325 SNPs), `rdadapt` K = 2 throughout:

| cell | `scale` | cond PCs | `pmax` | **AUC-PR** | best F1 | q<0.05 calls | TP | FP | precision |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|
| **D — current default** | TRUE | 5 | ✓ | **0.2223** | **0.302** | 0 | 0 | 0 | — |
| G | TRUE | 2 | ✓ | 0.1820 | 0.236 | 0 | 0 | 0 | — |
| A — *their uncorrected arm* | FALSE | 0 | — | 0.1353 | 0.205 | 416 | 15 | 86 | 0.149 |
| C | TRUE | 2 | — | 0.1285 | 0.170 | 0 | 0 | 0 | — |
| H | FALSE | 2 | ✓ | 0.1065 | 0.175 | 71 | 5 | 10 | **0.333** |
| E | TRUE | 5 | — | 0.0927 | 0.108 | 1 | 1 | 0 | 1.000 |
| B — *their corrected arm* | FALSE | 2 | — | 0.0813 | 0.120 | 408 | 10 | 122 | 0.076 |
| F | FALSE | 5 | — | 0.0666 | 0.115 | 394 | 12 | 169 | 0.066 |

Probe fidelity: cell A reproduces their published `RDA_AUCPR_neutSNPs` 0.1373 (we get 0.1353) and
cell B their `_corr` 0.0738 (0.0813) — same ±0.001–0.008 tolerance as every other cross-check.

**Each parameter isolated, holding the other two fixed:**

- **`pmax` helps in every context, and is the largest single effect.** E→D **+139%**
  (0.093→0.222), C→G **+41%**, B→H **+31%**. It also lifts thresholded precision sharply where
  thresholding fires at all: H = 0.333 versus their 0.076–0.149.
- **`scale=TRUE` helps in every context.** B→C **+58%**, H→G **+70%**, F→E **+39%**.
- **`condition_pcs` 2→5 is mixed and only pays off together with `pmax`.** G→D +22%, but
  C→E **−28%** and B→F **−17%**. This is the one genuinely open knob.

**Provisional recommendation — do NOT apply until replicated across the other 13 seeds:**

1. **Keep `scale=TRUE` and `pmax`.** Both are positive in every cell tested, and the current
   default D is the best of all eight on both AUC-PR and best-achievable F1.
2. **The real defect is calibration, not the model.** D has the best ranking of any cell yet
   reaches none of it through a standard threshold — its own best F1 sits at n=133 calls while
   q<0.05 yields 0. `max` of two ~uniform p-values is **not** uniform, so the combined value is a
   valid *ordering* but not a valid *p-value*. Options, in preference order:
   a. rank-based candidate selection for RDA (top-N, or an FDR against the empirical null)
      instead of q/Bonferroni on the `pmax` statistic;
   b. apply the GIF step to the *combined* statistic rather than to each fit separately;
   c. treat RDA as ranking-only and take candidates from its intersection with EMMAX/LFMM —
      which is Forester et al. 2018 §4.4's own recommendation for high-structure systems.
3. **Leave `condition_pcs` alone** until the 2-vs-5 question is tested on more seeds.

**Do not read their call counts as superiority.** Their arms reach TPR only at precision
0.076–0.149, matching their published FDR of 0.857/0.885 — 416 calls means 86 false positives
for 15 true ones. On this dataset the honest summary is that RDA ranks usefully and thresholds
badly, in both implementations.

### Where we currently stand (seed 1232548, K=5)

| Method | Axis | Published | Their p-values, our harness | **Our pipeline** |
|---|---|---:|---:|---:|
| LFMM | temp `bio_1` | 0.0783 ⚠ | 0.203 | **0.222** |
| LFMM | Env2 `bio_2` | 0.392 | 0.394 | **0.396** |
| RDA (no corr., 0 PCs) | multivariate | 0.137 | 0.140 | not run |
| RDA (struct. corr., 2 PCs) | multivariate | 0.074 | 0.081 | not run |
| RDA (ours: 5 PCs + `pmax`) | multivariate | — | — | **0.222** |
| Kendall AF~env | temp / Env2 | 0.106 / 0.117 | — | — |

At their own operating point (q < 0.05), **our LFMM TPR matches theirs exactly on both axes**:
temp **4/30 = 0.133** vs published **0.1333**; Env2 **8/40 = 0.200** vs published **0.2000**.
Our FDR is 0 on both (they report 0.200 / 0.000 — a one-SNP difference on the temperature axis).

Their phenotype GWAS is the reference for our `mode=gwas` arm and is a cautionary number:
**TPR 0.767 (temp) / 0.925 (Env2)** at **FDR 0.991 / 0.987** — nearly every causal locus found,
buried in ~99% false discoveries.

## Next: run and score, one seed at a time

```bash
D="nix shell nixpkgs#docker-client -c docker run --user $(id -u):$(id -g) --rm -e USER=pipeline \
   --cpus=64 --memory=300g -v /mnt/data/eugene/CLINE-GO:/pipeline cline-go:latest"
for m in processing prestructure structure climate pregea gea gwas; do
  $D snakemake -c16 -s Snakefile --config mode=$m --configfile config_MVP1232548.yaml --scheduler greedy
done
```

Never run two modes concurrently for one project — they share a Snakemake lock. Note `.snakemake/`
lives at the repo root and is shared across projects; if a run is killed, clear it with
`snakemake --unlock`, never by deleting the lock files.

Then score with the existing dataset-agnostic harness:

```bash
Rscript benchmarks/eval_detection.R --mode=rank \
  --pvalues=MVP1232548_results/GEA/tables/methods/LFMM/LFMM_pvalues_K5.tsv \
  --truth=data/mvp/MVP1232548/truth_temp.tsv --traits=bio_1 \
  --label=MVP1232548_LFMM --out=benchmarks/mvp_eval/MVP1232548.tsv
```

`benchmarks/score_ladders.sh` and `sweep_rda.sh` still default to the retired `LARUSON1K`
project and `data/laruson_1k/causal_loci.tsv`; retarget them at an `MVP{seed}` project and its
`truth_temp.tsv` before use.

The two `Est-Clines` controls (**1232568**, **1231578**) exist to be run, not just argued about:
with R²(PC1~temp) ≈ 0.92–0.94 they should reproduce the Láruson collapse *on purpose* once
structure correction is switched on. If they do not, the benchmark is not measuring what §8.4
claims.
