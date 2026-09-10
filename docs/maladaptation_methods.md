# Genetic Offset Methods — Research Reference

Candidate methods for expanding the maladaptation module beyond the current Gradient Forest baseline. This document is **literature and comparison only** — implementation details (Snakemake rules, R scripts, registry entries) are deferred to a separate session.

> Confidence markers: claims backed by adversarial multi-source verification from the deep-research harness are marked **[verified]**. Claims drawn from the primary literature without harness verification are marked **[literature]**. Both types cite primary sources.

---

## Background

Genetic (genomic) offset quantifies how much a population's current genotype-environment association deviates from the genotype-environment association expected under a projected future climate — a proxy for maladaptation risk. The GF-based offset already in the pipeline (Fitzpatrick & Keller 2015; Ellis et al. 2012) is the most widely used method but is not the only option. Running multiple complementary methods serves two goals:

1. **Robustness** — if two or more methods with different statistical assumptions converge on the same populations as high-risk, that convergence is more credible than any single method alone.
2. **Comparison / ensemble** — offset estimates from different methods can be rank-correlated, averaged, or used to bound uncertainty across methods.

All current genomic offset methods share a core assumption that cannot be fully tested: that the genotype-environment relationship observed today holds beyond the range of present climatic values. **[verified]** (Rellstab et al. 2021; Lind & Lotterhos 2025; Gain et al. 2023) This assumption is unverified for truly novel future climates; any offset estimate should be treated as a relative risk indicator, not an absolute prediction.

---

## Method Catalog

### Gradient Forest (GF) — existing baseline

**Principle.** Aggregated Random Forest regression of allele frequencies on environmental predictors across all loci. Turnover functions (importance-weighted cumulative R²) transform each environmental predictor into "genetic space." Offset = Euclidean distance in transformed genetic space between present and projected future environment. **[literature]** (Ellis et al. 2012; Fitzpatrick & Keller 2015)

**Inputs required.** Genotype matrix (samples × SNPs, 0/1/2 dosage); per-site present climate table; present and future climate raster (all cells). Can operate on whole-genome SNPs or GEA-significant SNPs.

**Extrapolation to novel climates.** Non-parametric (tree-based), so extrapolation uses "flat" turnover functions beyond the training range — better behaved than linear extrapolation but still unvalidated. **[verified]** (Rellstab et al. 2021)

**Simulation benchmarks.** Broadly correlated with fitness under both monogenic and polygenic architectures in SLiM simulations. **[verified, high]** (Láruson et al. 2022, DOI: 10.1111/eva.13354) Consistently outperforms F~ST~-based differentiation as a fitness predictor. **[verified, high]** Three key confounders identified:
- Neutral demographic effects — GF offset is negatively correlated with deme size because genetic drift in small demes generates spurious allele-frequency turnover that GF mistakes for adaptive signal. **[verified, high]**
- Nonlinear environments — performance weakest when causal environmental predictors are not the highest-ranked predictors across all loci. **[verified, high]**
- General confounders — neutral demography, genomic architecture, and adaptive environment all confound offset-fitness relationships. **[verified, high]** (Láruson et al. 2022)

**Common garden validation.** GF offset had a significant negative effect on red spruce (*Picea rubens*) seedling height growth (p < 0.001) across all tested climate scenarios in three common gardens (Vermont, Maryland, North Carolina). **[verified, high]** (Lachmuth et al. 2023, DOI: 10.3389/fevo.2023.1155783) This is partial empirical validation; independent multi-source validation is required before conservation applications. **[verified, high]** (Archambeau et al. 2026, DOI: 10.1086/739111)

**Adaptive vs. whole-genome SNP sets.** Using GEA-identified "adaptive" markers provides only marginal advantages: median performance increase < 3%, advantages in only 67–68% of comparisons, and advantages are least prevalent precisely under climate novelty conditions. **[verified, high]** (Lind & Lotterhos 2025, DOI: 10.1111/1755-0998.14008)

**R package.** `gradientForest` (Ellis et al.); `gradientForestOffset` wrapper available.

**Status in pipeline.** Implemented. Runs on curated SNP sets from the Shiny GEA tab. Supports spatial correction (PCNM) and optional neutral random model.

---

### GDM Offset — Generalized Dissimilarity Modelling

**Principle.** Generalized Dissimilarity Modelling uses I-spline-weighted regression of pairwise genetic dissimilarity on environmental and geographic distances between pairs of sites. **[literature]** (Fitzpatrick & Keller 2015, DOI: 10.1111/ele.12376) Offset is computed by transforming present-climate site values through GDM's I-spline functions, then taking Euclidean distance in transformed space between present and projected future environments. Conceptually similar to GF but operates on pairwise dissimilarities rather than individual locus×sample genotypes.

**Inputs required.** Population-level allele frequencies (site × allele) or individual genotype matrix aggregated to sites; pairwise geographic distances or lat/long; per-site climate table (present and future). **[literature]**

**Extrapolation.** I-splines are piecewise monotone — extrapolation beyond training range uses the terminal slope of the spline, which is more constrained than GF tree-flattening but still unvalidated at true climate novelty. **[literature]**

**Simulation benchmarks.** GDM offset was established alongside GF in Fitzpatrick & Keller 2015 (balsam poplar). **[verified, medium]** One unverified-source quote from the harness (source access failed) suggested "GDM and GF produced the strongest and most consistent associations with empirical mortality data" in maritime pine (Archambeau et al. 2024, preprint). Direct head-to-head simulation benchmark vs. GF and geometric GO is not yet available from verified sources.

**Pros.** Handles nonlinear environment-genotype relationships via I-splines. Naturally operates at population/site level (appropriate when individual assignment is uncertain). Mature package.

**Cons.** Requires per-site (population-level) allele frequencies — **not currently produced by the CLINE-GO pipeline**; would need a new rule to average dosages within sites. Pairwise structure scales as O(n²) samples. Less widely benchmarked than GF for offset specifically.

**R package.** `gdm` (CRAN); key functions: `gdm::gdm()`, `gdm::gdm.transform()`, `gdm::predict.gdm()`.

**Key citation.** Fitzpatrick MC & Keller SR (2015) *Ecology Letters* 18:1–11. DOI: 10.1111/ele.12376

---

### RDA Genetic Offset — Redundancy Analysis

**Status: implemented** as a batch maladaptation method (`workflow/methods/maladaptation.py`, `scripts/rda_offset.R`). Full decision table, gotchas, and citations: `docs/rda_research.md` Part B (row IDs B1–B22 referenced below). This section is a condensed summary reconciled against that primary spec — read Part B before modifying `rda_offset.R`.

**Principle.** A *second* constrained ordination (RDA), independent of the GEA-scan RDA (`scripts/rda.R`) — "compute once, reuse downstream" does not apply here (B0): the two fits differ on both the SNP set (candidate-only vs. all) and the conditioning term (none vs. partial) simultaneously. The offset RDA fits `rda(candidate_SNPs ~ climate)` with **no** `Condition()` by default (B7 — canonical construction). Present and future RDA scores are projected via `method="loadings"` (`Σ_v z_v × CCA$biplot[v,i]`, B4) — **never** `predict(type="lc")`, which is broken in both of Capblancq & Forester's own shipped source files (gotcha G1). Offset = eigenvalue-weighted Euclidean distance between predicted present and future scores, weights applied to the scores *before* squaring in `dist()` (B2, B3). (Capblancq & Forester 2021, DOI: 10.1111/2041-210X.13722)

Two variants exist in the literature — CLINE-GO implements the uncorrected form by default, with the corrected form available as a documented deviation:
- **RDA-uncorrected** (default, `condition_pcs: 0`) — canonical construction per B7: conditioning removes exactly the climate-correlated structure variance the adaptive index needs to project.
- **RDA-corrected (partial RDA)** — `condition_pcs > 0`; a labeled deviation from the canonical method, retained as a config option because Lind et al. 2025 treat corrected/uncorrected as two distinct implementations with no declared winner overall (B.6 open disagreement #2).

Lind & Lotterhos 2025 (*Mol Ecol Resour* 25(4):e14008, >4.8M simulation evaluations) find performance driven mainly by degree of local adaptation, not method choice — adaptive (candidate) marker sets give only a **minimal** advantage over whole-genome panels (median <3% gain), least prevalent under climate novelty, the exact condition where offset is applied. This directly informs the candidate-set-only design here (B5) and the general caution against over-interpreting marker-set choice.

Gain et al. 2023 prove RDA offset ≡ geometric offset **when RDA includes latent (structure) predictors alongside environment**. Since CLINE-GO's offset RDA omits latent predictors by default (B7), the two methods are expected to differ on real data — by how much is an empirically testable, currently unquantified question (B15). The pipeline's `compare_offsets.R` / Shiny cross-model comparison tab (Spearman ρ, Jaccard top-K overlap, ExDet novelty stratification) measures this directly, with no extra code required once a method is registered.

Gain et al. 2023 also found linear methods (geometric GO, RDA) achieve a better bias-variance trade-off than machine-learning methods such as GF, an explanation for cases where linear methods outperform GF at limited sample sizes. **[verified, medium — paper uses hedged language "An explanation may be..."]**

**Inputs required.** Full imputed genotype matrix (`W['lfmm_imp_full_climate']`), subset internally to the curated `{run_label}` candidate-SNP set (`snp_set_file`, same mechanism as `geometric_offset`); present/future climate at sample sites and raster cells (`O['climate_site']`, `O['climate_future_site']`, `O['climate_all']`, `O['climate_future_all']`). All already produced by the pipeline — **no pipeline gap**.

**Extrapolation.** Linear extrapolation of RDA loadings to novel climates — same failure class as the geometric offset (unbounded, no plateauing the way GF's tree-based construction has). **[literature]** (B16)

**Simulation benchmarks.** Gain et al. 2023 compared geometric GO and RDA favourably against GF and squared Euclidean distance under Gaussian stabilizing selection. Lind & Lotterhos 2025's large-simulation study covered GF, RDA, LFMM2, and RONA; the adaptive-marker finding above is its most directly applicable result to this implementation. Archambeau et al. 2026 (*Am Nat* 207(3):389–414, maritime pine, 5 common gardens) found substantial cross-method variability — best-validated methods gave *differing* geographic patterns, motivating multi-method concordance (B19) over trust in any single map.

**Interactivity, resolved.** The original plan deferred RDA offset to a fully interactive, Shiny-driven flow because axis selection and the structure-conditioning choice need human judgment. As implemented, that judgment is exercised **once, upstream**, through the existing curated-SNP-set mechanism (the GEA tab's "Save SNP set" action, shared with `geometric_offset`) plus config-level `axes`/`condition_pcs` — not per-run in the maladaptation module itself. Axis count defaults to `"auto"` (anova by-axis at `axis_alpha`, floored at 2 constrained axes — a 1-axis offset degenerates to a single scaled climate distance, the "genomics added nothing" case B.5 warns against).

**Pros.** Linear, interpretable, theoretically grounded (Gain et al. 2023 unification). Fast — a second fit on a small candidate-SNP matrix, not the full-marker-set scale problem the GEA scan has (A5/A6). Directly reuses inputs already in the pipeline. Registered as a standard batch maladaptation method — same fan-out over SNP sets as GF/GeoOff, `nospatial`-only (B7).

**Cons.** Linear extrapolation assumption (B16). Candidate-set choice (B5, B6) is a documented deviation from the canonical partial∩unconstrained intersection — CLINE-GO uses the curated `run_label` set instead, justified by Lind et al. 2025's marker-set-barely-matters finding. Offset should be reported as a **relative site ranking only**, never an absolute magnitude (B20) — eigenvalue weighting is an implementation convention, not derived from an explicit fitness model (contrast the geometric offset's Gaussian-selection derivation). Unvalidated absent common-garden data (B21).

**R package / functions.** `vegan::rda()`, `vegan::anova.cca()`; offset projection ported from Capblancq & Forester 2021's `genomic_offset.R` (`scripts/rda_offset.R`, `terra`-based).

**Key citations.**
- Capblancq T & Forester BR (2021) Redundancy analysis: A Swiss Army Knife for landscape genomics. *Methods in Ecology and Evolution* 12(12):2298–2309. DOI: 10.1111/2041-210X.13722
- Capblancq T, Fitzpatrick MC, Bay RA, Exposito-Alonso M & Keller SR (2020) Genomic prediction of (mal)adaptation across current and future climatic landscapes. *Annual Review of Ecology, Evolution, and Systematics* 51: 245–269. DOI: 10.1146/annurev-ecolsys-020720-042553
- Gain C, Rhoné B, Cubry P, Salazar I, Forbes F, Vigouroux Y, Jay F & François O (2023) *Molecular Biology and Evolution* 40(6):msad140. DOI: 10.1093/molbev/msad140
- Lind BM & Lotterhos KE (2025) *Molecular Ecology Resources* 25(4):e14008. DOI: 10.1111/1755-0998.14008
- Archambeau et al. (2026) *American Naturalist* 207(3):389–414. DOI: 10.1086/739045

---

### Geometric Genetic Offset / Genetic Gap (LEA)

**Principle.** Theoretically grounded in Gaussian stabilizing selection. Uses LFMM2 regression effect-size coefficients (β matrix, loci × predictors) as a genetic-effect-weighted covariance matrix. Offset between present environment **x** and future environment **x*** is a Mahalanobis-type quadratic distance in environmental space: G² ∝ (x − x*) Σ_β (x − x*)^T, where Σ_β is derived from the LFMM2 β matrix. **[literature]** (Gain et al. 2023) The exact formula claim was killed by the adversarial harness (1-2 vote on precision), so the formula above is presented as conceptual rather than verified-verbatim.

Key theoretical results from Gain et al. 2023 (MBE): geometric GO achieves r² ≈ 78% correlation with true fitness offset in forward simulations under Gaussian stabilizing selection, compared to r² ≈ 45% for squared Euclidean environmental distance. **[verified, medium]** The r² = 97% fit of the offset to a quadratic function of selection intensity provides evidence for the method's theoretical grounding. **[verified, medium]** In pearl millet common garden experiments, geometric GO achieved r² = 61% — lower than the simulation ideal but still best-performing among methods tested. **[literature]** (Gain et al. 2023)

**Inputs required.** Genotype matrix (LFMM format, samples × SNPs); per-site or per-sample present climate table; future climate table. LFMM2 is run internally by `LEA::genetic.gap()`. **[literature]** All these inputs are **already produced by the CLINE-GO pipeline** (`W['lfmm_full']`, `O['climate_site']`, future climate tables) — this method has the cleanest data contract fit of any new candidate.

**Extrapolation.** Linear projection of LFMM effect sizes to future environment. Same linearity caveat as RDA offset. **[literature]**

**Simulation benchmarks.** r² ≈ 78% vs. fitness under Gaussian stabilizing selection; outperforms Euclidean environmental distance. **[verified, medium]** Linear methods including geometric GO achieve better bias-variance tradeoff than GF in data-limited conditions. **[verified, medium, hedged]** (Gain et al. 2023)

**Pros.** Strong theoretical grounding (unified framework with RDA). Easiest data contract for CLINE-GO (LFMM inputs already exist). Single R function call. Accounts for polygenic architecture through the full β matrix.

**Cons.** Assumes Gaussian stabilizing selection at equilibrium. Performance in real non-equilibrium populations with drift likely lower than simulation ideals. LFMM2 runtime scales with SNP count (can be slow at WGS density without filtering).

**R package.** `LEA` (Bioconductor); key function: `LEA::genetic.gap(input, env, pred.env, K)` where `env` is the present-climate matrix, `pred.env` is the future-climate matrix, and K is the number of latent factors (same K as used in LFMM2 GEA). **[verified]**

**Key citation.** Gain C, Rhoné B, Cubry P, Salazar I, Forbes F, Vigouroux Y, Jay F & François O (2023) *Molecular Biology and Evolution* 40(6): msad140. DOI: 10.1093/molbev/msad140

---

### RONA — Deferred

**Why deferred.** RONA (Risk of Nonadaptedness; Rellstab et al. 2016, DOI: 10.1111/mec.13809) estimates the expected allele-frequency change per unit climate shift via a linear regression of allele frequency on each environmental predictor, then sums over a predefined set of putatively adaptive loci. The key extrapolation problem is mechanistic: RONA assumes the linear regression holds outside the range of present climate values observed across sampled populations. When future climates are truly novel (outside this range), linear extrapolation becomes arbitrary — the regression is used as a prediction function far beyond its training domain. Additionally, RONA is inherently a candidate-locus method (cannot aggregate across a genome-wide SNP set in a principled way) and does not account for correlations among predictors or among loci. **[literature; direct harness verification of specific mechanism claims was unsuccessful — all RONA-specific claims submitted to adversarial harness received 0-3 votes, likely due to source access failures rather than claim incorrectness]**

RONA will not be added to the pipeline. It is documented here for completeness.

**Key citation.** Rellstab C et al. (2016) *Molecular Ecology* 25:5479–5495. DOI: 10.1111/mec.13809

---

### Modern / Post-2023 Survey

**algatr R package.** A unified landscape genomics analysis pipeline implementing six methods: TESS, wingen, MMRR, GDM, RDA, and LFMM — covering population structure, diversity, IBD/IBE, and GEA. **[verified, high]** (Chambers, Bishop & Wang 2025, Mol Ecol Resour, DOI: 10.1111/1755-0998.13884, PMC12142729) Importantly, **algatr does not implement offset extensions** — the GEA and dissimilarity methods are present, but the maladaptation-projection layer (GF offset, GDM offset, RDA offset, geometric GO, RONA) is absent. **[verified, high]** algatr is thus useful for the GEA step (upstream of maladaptation) but not for the offset computation itself.

**Mahalanobis-type offsets.** Methods computing offset as Mahalanobis distance in environmental space (weighting by observed genetic covariance across environments) have been proposed but the deep-research harness did not locate verified primary-source benchmark claims for a distinct "Mahalanobis offset" method separate from the geometric GO. The geometric GO (Gain et al. 2023) is itself Mahalanobis-type, and represents the current best-supported implementation of this class.

**Deep-learning offsets.** No published, benchmarked deep-learning genomic offset method was identified in the harness's search through June 2026. Deep learning has been applied to population structure and ancestry inference but not yet to genetic offset in a peer-reviewed benchmark.

**Ensemble approaches.** No verified claim for a formal ensemble offset method (e.g., rank-aggregation or weighted average of GF + geometric GO + RDA + GDM outputs) was recovered. However, the Archambeau et al. 2026 study (*American Naturalist*, DOI: 10.1086/739111) demonstrates that "validation of GO predictions with a range of independent data sources" is necessary before conservation use **[verified, high]** — which implicitly endorses multi-method comparison as a validation strategy. Rank correlation among methods (Spearman ρ across sites) is a natural diagnostic: high concordance strengthens inference; disagreement flags populations where method assumptions diverge.

---

## Benchmark Synthesis

Summary of evidence from simulation studies and common garden experiments. Blank cells = not tested or source not accessible.

| Method | Simulation benchmark | Fitness correlation | Outperforms FST | Key confounders | Novel-climate performance | Common garden | Package | Maturity |
|--------|---------------------|--------------------|-----------------|-----------------|-----------------------------|---------------|---------|---------|
| **Gradient Forest** | Láruson et al. 2022 (SLiM) | Broadly correlated [V] | Yes [V] | Drift (deme size), nonlinear env, demography [V] | Degrades with novelty [V] | Red spruce p<0.001 [V] | `gradientForest` | Production |
| **Geometric GO (LEA)** | Gain et al. 2023 (forward sim) | r²≈78% vs. fitness [V] | — | Equilibrium / Gaussian stab-sel assumption [L] | Linear extrapolation [L] | Pearl millet r²≈61% [L] | `LEA::genetic.gap()` | Mature |
| **RDA offset** | Gain et al. 2023, Lind & Lotterhos 2025 (4.85M evals) | — | — | Linearity; structure correction choice [L] | Linear extrapolation [L] | Not available from verified sources | `vegan` + custom (`scripts/rda_offset.R`) | **Implemented** |
| **GDM offset** | Fitzpatrick & Keller 2015 [V-med] | — | — | Requires per-site allele freq [L] | I-spline extrapolation [L] | Not available from verified sources | `gdm` | Mature |
| **RONA** | Lind & Lotterhos 2025 (4.85M evals) | Lower than GF/RDA [UV] | No [UV] | Single-locus; linear out-of-range [L] | Breaks at climate novelty [L] | — | Custom | Deferred |

**Confidence key:** [V] = adversarially verified (≥2/3 votes), [L] = literature knowledge (pre-Jan 2026 training), [UV] = unverified (source not accessible by harness — treat as indicative only).

**Notable benchmark result (Lind & Lotterhos 2025, ~4.85 million simulation evaluations):** Adaptive (GEA-identified) SNP sets provide **minimal** performance advantages over whole-genome panels — median gain < 3%, observable in only 67–68% of comparisons, and **least prevalent under climate novelty** — the exact conditions where offset is applied. **[verified, high]** This finding applies across GF, RDA, LFMM2, and RONA and weakens the rationale for restricting any offset method to GEA-significant SNPs only.

---

## Recommendation and Next Steps

### Priority order for implementation

1. **Geometric Genetic Offset (`LEA::genetic.gap()`)** — **implemented** (`scripts/geometric_offset.R`). Cleanest data contract: LFMM input format (`W['lfmm_full']`), `O['climate_site']`, and future climate tables are all already produced. Single function call. Strong theoretical grounding (Gain et al. 2023).

2. **RDA Offset** — **implemented** (`scripts/rda_offset.R`, `workflow/methods/maladaptation.py`), after RDA GEA (`scripts/rda.R`) landed. A second, candidate-SNP-only, uncorrected-by-default RDA fit — see the RDA Genetic Offset section above and `docs/rda_research.md` Part B for the full design. Registered as a standard batch method (fan-out over curated SNP sets, `nospatial`-only), not a fully bespoke interactive flow — the axis-selection / structure-conditioning judgment calls are exercised once via config (`axes`, `condition_pcs`), reusing the same curated-SNP-set mechanism as `geometric_offset`.

3. **GDM Offset** — lower priority, not yet implemented. Requires a new pipeline intermediate: per-site allele-frequency table (currently not produced). GDM is already in algatr's GEA layer, so the upstream step exists conceptually; the offset projection would be new. Would complete the benchmark triangle (nonparametric GDM vs. nonparametric GF vs. linear geometric GO vs. linear RDA).

### Comparing and combining offset estimates

- **Rank correlation (Spearman ρ)** across sites between pairs of methods is the primary diagnostic. High concordance (ρ > 0.7–0.8) strengthens interpretation. Disagreement between methods identifies populations where method-specific assumptions are most influential.
- **Rank-based ensemble** (average ranks across methods, then re-rank) is robust to differences in scale and distribution between methods. Appropriate when all implemented methods are considered equally credible.
- **Uncertainty bounds** — display per-site offset range (min, max across methods) alongside the ensemble mean. Populations consistently high across all methods are the highest-confidence priority areas.
- **No single method dominates** across all simulation conditions. **[literature; unverified-source confirmation]** Multi-method consensus is a more defensible basis for conservation recommendations than any single estimate.

### Data contract compatibility (CLINE-GO pipeline)

| Method | Present climate | Future climate | Genotype input | Per-site allele freq | Pipeline gap |
|--------|----------------|---------------|----------------|---------------------|--------------|
| GF (existing) | `O['climate_all']` | `O['climate_future_all']` | `W['lfmm_full']` + `W['vcfsnp_full']` | Not needed | None |
| Geometric GO | `O['climate_site']` | `O['climate_future_site']` | `W['lfmm_full']` | Not needed | **None** — ready |
| RDA offset | `O['climate_site']` | `O['climate_future_site']` | `W['lfmm_imp_full_climate']` (subset to candidate SNPs) | Not needed | **None** — ready (implemented) |
| GDM offset | `O['climate_site']` | `O['climate_future_site']` | Derived from `W['lfmm_full']` | **Required — not produced** | New rule: allele-freq-per-site |

---

## Open Questions (from deep-research harness)

These questions were not answered by verified primary sources and require further targeted literature search:

1. **RONA extrapolation mechanism (precise)** — every harness claim about RONA's linear-regression-per-locus assumption and its failure under novel climates received 0–3 votes (likely source access failures, not claim incorrectness). The primary Rellstab et al. 2016 and 2021 papers should be read directly to document the mechanism precisely.

2. **RDA offset vs. GDM offset vs. GF vs. geometric GO — direct simulation comparison** — no single benchmark comparing all four under common conditions survived adversarial verification. The Lind & Lotterhos 2025 paper (DOI: 10.1111/1755-0998.14008) likely contains this comparison but source access failed for that specific result. Read directly.

3. **Ensemble performance** — no quantitative result on whether rank-aggregated multi-method ensembles outperform individual methods was verified. Archambeau et al. 2026 and Lotterhos 2024 (*Evolutionary Letters* 8:331) emphasize multi-source validation but not ensemble statistics specifically.

4. **Post-2023 methods** — no published, benchmarked deep-learning or novel Mahalanobis-type offset method was identified beyond geometric GO. The field is active; a re-search in 2026–2027 is warranted.

---

## References

**Verified primary sources (adversarially verified by deep-research harness):**

- Fitzpatrick MC & Keller SR (2015) Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. *Ecology Letters* 18: 1–11. DOI: 10.1111/ele.12376

- Láruson AJ, Haller BC, Keller SR, Fitzpatrick MC & Lotterhos KE (2022) Seeing the forest for the trees: assessing genetic offset predictions from Gradient Forest. *Evolutionary Applications* 15: 403–416. DOI: 10.1111/eva.13354 PMC8965365

- Gain C, Rhoné B, Cubry P, Salazar I, Forbes F, Vigouroux Y, Jay F & François O (2023) A quantitative theory for genomic offset statistics. *Molecular Biology and Evolution* 40(6): msad140. DOI: 10.1093/molbev/msad140

- Lind BM & Lotterhos KE (2025) The accuracy of predicting maladaptation to new environments with genomic data. *Molecular Ecology Resources* 25: e14008. DOI: 10.1111/1755-0998.14008 PMC11969643

- Lachmuth S, Capblancq T, Keller SR & Fitzpatrick MC (2023) Novel genomic tools for the assessment of climate-driven maladaptation of red spruce (Picea rubens). *Frontiers in Ecology and Evolution*. DOI: 10.3389/fevo.2023.1155783

- Archambeau J et al. (2026) Reciprocal evaluation of genomic offset predictions of climate maladaptation with independent empirical datasets. *American Naturalist*. DOI: 10.1086/739111 Preprint: bioRxiv 2024.05.17.594631

- Rellstab C et al. (2021) Genomics helps to predict maladaptation. *Evolutionary Applications*. PMC8127717

- Chambers EA, Bishop AP & Wang IJ (2025) Individual-based landscape genomics for conservation: an analysis pipeline (algatr). *Molecular Ecology Resources* 25: e13884. DOI: 10.1111/1755-0998.13884 PMC12142729

**Literature sources (not harness-verified, cited from primary knowledge):**

- Ellis N, Smith SJ & Pitcher CR (2012) Gradient forests: calculating importance gradients on physical predictors. *Ecology* 93: 156–168. DOI: 10.1890/11-0252.1

- Capblancq T & Forester BR (2021) Redundancy analysis: a Swiss army knife for landscape genomics. *Methods in Ecology and Evolution* 12: 2298–2309. DOI: 10.1111/2041-210X.13722

- Capblancq T, Fitzpatrick MC, Bay RA, Exposito-Alonso M & Keller SR (2020) Genomic prediction of (mal)adaptation across current and future climatic landscapes. *Annual Review of Ecology, Evolution, and Systematics* 51: 245–269. DOI: 10.1146/annurev-ecolsys-020720-042553

- Rellstab C et al. (2016) Signatures of local adaptation in candidate genes of oaks (*Quercus* spp.) in comparison to their neutral landscape genomic structure. *Molecular Ecology* 25: 5444–5457. DOI: 10.1111/mec.13809 [RONA original]

- Fitzpatrick MC & Keller SR (2015) — same as above; GDM also introduced here.
