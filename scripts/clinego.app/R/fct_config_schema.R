#' Config parameter schema for all pipeline tabs
#'
#' Each entry describes one config parameter: where it lives in the YAML,
#' how to display it, which tab it belongs to, and validation rules.
#'
#' Fields:
#'   key        dot-path into config (e.g. "sNMF.k_start")
#'   label      human-readable label
#'   tab        tab name: home | processing | prestructure | structure | gea |
#'                         gwas | gea_x_gwas | haplotype | maladaptation
#'   section    grouping label within the tab (subsection when `group` is set)
#'   group      optional top-level grouping label. When ANY entry on a tab sets
#'              `group`, that tab renders as flat named sections (one per group,
#'              `section` becomes a subsection label within it) instead of the
#'              default mandatory-core / Advanced-accordion split.
#'   type       numeric | text | select | checkbox | checkbox_invert | textarea |
#'                       method_table
#'              checkbox_invert displays/writes the negation of the stored value
#'              (e.g. a "Disable X" box backed by an `X.enabled` key)
#'   mandatory  TRUE = shown in core section; FALSE = shown in Advanced accordion.
#'              Honored on `group` tabs too — the Advanced accordion nests
#'              inside each group rather than at the tab's top level.
#'   show_if    optional dot-key of a checkbox entry; this field only renders
#'              while that checkbox is checked (client-side conditionalPanel)
#'   help       tooltip / hint text shown below input
#'   placeholder hint text for empty text inputs
#'   min/max/step  for numeric type
#'   choices    for select type (character vector)
#'
#' @noRd
config_schema <- function() {
    # helper to build one entry compactly
    s <- function(key, label, tab, section, type, mandatory, ..., group = NULL, show_if = NULL) {
        args <- list(...)
        list(
            key         = key,
            label       = label,
            tab         = tab,
            section     = section,
            group       = group,
            type        = type,
            mandatory   = mandatory,
            show_if     = show_if,
            help        = args$help,
            placeholder = args$placeholder,
            min         = args$min,
            max         = args$max,
            step        = args$step,
            choices     = args$choices
        )
    }

    list(
        # ── HOME TAB (project input files only) ───────────────────────────────
        s("Input.dir",      "Input directory",    "home", "Input Files",
          "text",    TRUE,
          help = "Directory containing input files (relative to /pipeline/)",
          placeholder = "data/"),
        s("Input.vcf",      "VCF file",           "home", "Input Files",
          "text",    TRUE,
          help = "VCF filename inside input dir (.vcf or .vcf.gz)",
          placeholder = "samples.vcf.gz"),
        s("Input.metadata", "Metadata file",      "home", "Input Files",
          "text",    TRUE,
          help = "Tab-separated: site, sample, lat, lon, [phenotype columns…]",
          placeholder = "metadata.tsv"),
        s("Input.gff",      "GFF3 annotation",    "home", "Input Files",
          "text",    TRUE,
          help = "Optional. Enables gene annotation and GO enrichment.",
          placeholder = "annotation.gff3"),

        # ── PROCESSING TAB ────────────────────────────────────────────────────
        # Two top-level groups (SNP Filtering / Sample Filtering), no Advanced
        # accordion — every processing filter is a first-look decision, not a
        # rarely-touched knob. LD Pruning is a subsection within SNP Filtering.
        s("Filter.maf",         "Minor allele freq",  "processing", NULL,
          "numeric", TRUE, group = "SNP Filtering",
          min = 0, max = 0.5, step = 0.01,
          help = "Minimum MAF threshold (e.g. 0.05)"),
        s("Filter.snp_miss",    "SNP missingness",    "processing", NULL,
          "numeric", TRUE, group = "SNP Filtering",
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per SNP"),
        s("LD.window", "LD window (kb)",  "processing", "LD Pruning",
          "numeric", TRUE, group = "SNP Filtering",
          min = 1, step = 1,
          help = "Sliding window size for LD pruning in kilobases"),
        s("LD.step",   "LD step size",   "processing", "LD Pruning",
          "numeric", TRUE, group = "SNP Filtering",
          min = 1, step = 1,
          help = "Number of SNPs to slide per step"),
        s("LD.r2",     "LD r\u00b2 threshold", "processing", "LD Pruning",
          "numeric", TRUE, group = "SNP Filtering",
          min = 0, max = 1, step = 0.05,
          help = "Prune SNP pairs in LD above this r\u00b2 threshold"),
        s("Filter.sample_miss", "Sample missingness", "processing", NULL,
          "numeric", TRUE, group = "Sample Filtering",
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per sample (default 0.5)"),
        s("Filter.relatedness", "Max relatedness (IBS)", "processing", NULL,
          "numeric", FALSE, group = "Sample Filtering",
          min = 0, max = 1, step = 0.05,
          help = "Optional. Colors/counts related pairs (plink IBS allele-sharing) above this in the Processing relatedness histogram (blank = skip). Does NOT remove samples by itself \u2014 see 'Relatedness action' below. IBS is model-free \u2014 works for selfers and outcrossers alike; duplicates/clones cluster near 1.0 \u2014 set from the histogram, not a fixed outbred default."),
        s("Filter.relatedness_action", "Relatedness action", "processing", NULL,
          "select",  FALSE, group = "Sample Filtering",
          choices = c("keep", "remove"),
          help = "keep (default): only visualize/count related pairs, remove nothing. remove: drop the higher-missingness member of each pair above the threshold. Recommended: set the threshold, inspect the relatedness histogram, then switch to remove and re-run Processing."),

        # ── PRESTRUCTURE TAB ──────────────────────────────────────────────────
        s("sNMF.k_start", "K range start",  "prestructure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Minimum K (ancestral populations) to test in cross-validation"),
        s("sNMF.k_end",   "K range end",    "prestructure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Maximum K to test in cross-validation"),
        s("sNMF.repeats", "Repeats per K",  "prestructure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 100, step = 1,
          help = "Number of independent sNMF runs per K for robust cross-entropy estimation"),

        # ── STRUCTURE TAB ─────────────────────────────────────────────────────
        s("sNMF.k_best",           "Best K",            "structure", "Structure",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Selected K for all downstream analysis. Review cross-entropy plot first."),
        # Climate.enabled is no longer user-editable here -- it is derived from
        # Regime.mode on every write (fct_config_writer.R). A project declares its
        # regime at creation; see config_regime() in fct_config.R.
        s("Climate.predictors",    "Climate predictors", "gea", "Climate",
          "bio_chips", TRUE,
          help = "Click to toggle variables for GEA. Review piemaps in Structure first to drop collinear or low-variance predictors. Same list as the Climate tab's copy of this field — editing either updates both."),
        # Optional — map
        s("Map.climate_extent", "Map extent",        "structure", "Map",
          "text",    TRUE,
          help = "Bounding box [min_lon,max_lon,min_lat,max_lat] or 'auto'",
          placeholder = "auto"),
        s("Map.gap",            "Auto extent gap",   "structure", "Map",
          "numeric", FALSE,
          min = 0, max = 10, step = 0.5,
          help = "Degrees of buffer around sample coordinates (auto extent only)"),
        s("Map.resolution",     "Resolution (min)",  "structure", "Map",
          "numeric", TRUE,
          min = 0.5, max = 10, step = 0.5,
          help = "WorldClim resolution in arc-minutes (lower = finer but slower)"),
        s("Map.zoom_extent",    "Zoom extent",       "structure", "Map",
          "text",    FALSE,
          help = "Optional zoom region for piemaps: xmin,xmax,ymin,ymax",
          placeholder = "NULL"),
        # Population stats — core (not Advanced): dependents reveal only once
        # calc_stats is checked, so the section stays compact when unused.
        s("Population.calc_stats",       "Calculate pop stats", "structure", "Pop Stats",
          "checkbox", TRUE,
          help = "Compute Tajima's D, Pi diversity, AMOVA, and IBD analysis"),
        s("Population.window_size",      "Window size (bp)",   "structure", "Pop Stats",
          "numeric",  TRUE,
          min = 1000, step = 1000,
          show_if = "Population.calc_stats",
          help = "Sliding window size in base pairs for Tajima's D and Pi diversity"),
        s("Population.custom_trait_file","Custom trait file",  "structure", "Pop Stats",
          "text",     TRUE,
          show_if = "Population.calc_stats",
          help = "Custom trait file for piemap sizing (NULL = use metadata phenotype columns)",
          placeholder = "NULL"),
        # Optional — piemap. Order: size/appearance first, then the pie->points
        # style switch, then label controls (label size gated on show_labels).
        s("Piemap.pie_scale",   "Pie size scale",    "structure", "Piemap",
          "numeric", FALSE,
          min = 0.1, max = 5, step = 0.1,
          help = "Multiplier for pie chart diameter (1.0 = default)"),
        s("Piemap.alpha",       "Pie transparency",  "structure", "Piemap",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Transparency of pie slices (0 = transparent, 1 = opaque)"),
        s("Piemap.use_points",  "Replace pies with points", "structure", "Piemap",
          "checkbox", FALSE,
          help = "Draw sample sites as simple dots instead of pie charts — clearer when dense sampling makes pies overlap. Affects on-demand haplotype-viz maps; the main Structure/Maladaptation/GWAS maps always emit both and switch live via the Points toggle in the app."),
        s("Piemap.show_labels", "Show pop labels",   "structure", "Piemap",
          "checkbox", FALSE,
          help = "Overlay population name labels on piemaps"),
        s("Piemap.label_size",  "Label size",        "structure", "Piemap",
          "numeric",  FALSE,
          min = 4, max = 20, step = 1,
          show_if = "Piemap.show_labels",
          help = "Font size for population labels"),
        # Optional — LD decay
        s("LDdecay.group_by",    "Group by",        "structure", "LD Decay",
          "select", TRUE,
          choices = c("site", "cluster"),
          help = "Variable to group samples by for LD decay curves"),
        s("LDdecay.min_samples", "Min samples",     "structure", "LD Decay",
          "numeric", TRUE,
          min = 2, step = 1,
          help = "Minimum group size; smaller groups are excluded"),
        s("LDdecay.max_distance","Max dist (kb)",   "structure", "LD Decay",
          "numeric", TRUE,
          min = 10, step = 50,
          help = "Maximum pairwise distance in kb for LD calculations"),
        s("LDdecay.scope",       "Scope",           "structure", "LD Decay",
          "select", TRUE,
          choices = c("both", "genome_wide", "per_chromosome"),
          help = "Compute LD decay genome-wide, per chromosome, or both"),

        # ── CLIMATE TAB (predictor characterization + spatial/varpart) ────────
        # Split out of PreGEA (module split, see CLAUDE.md "Phase" notes):
        # correlation/density/invariant-predictor display now lives in
        # mod_climate.R (producers unchanged, in structure.smk); this tab's
        # config covers the variance-partitioning sub-block that used to live
        # under PreGEA.Varpart, plus (2026-08-03) its own copy of the
        # Climate.predictors chip picker — SAME config key as the GEA tab's
        # entry below, both write into the one shared config_state$working
        # path, so editing on either tab is instantly reflected on the other
        # (Filter(e$tab==tab_name, config_schema()) just scopes which tabs a
        # given key renders on; project/invariant-drop plumbing in
        # render_bio_chips is already tab-agnostic). Added because varpart
        # itself consumes Climate.predictors directly (X_clim) and this tab
        # has no other way to fix an invariant predictor (e.g. bio_14)
        # silently NaN-cascading climate_varpart to an empty
        # variance_partition.tsv — see
        # knowledge/ADAPTOGENE/55-rda-scale-gate-and-invariant-predictor-cascade.
        #
        # Model fitted: Y ~ X_clim + X_struct + X_geo. Y (response, this
        # section) is genomic PCs by default — a documented scalability
        # adaptation vs. Lasky et al. 2012's raw-SNP response (WGS-scale
        # vegan::rda on millions of columns is impractical); 'snps' recovers
        # literal Lasky methodology. X_clim/X_struct/X_geo (explanatory
        # tables, next section) are NOT response-matrix choices: X_clim is
        # Climate.predictors (chip picker below), X_struct is this section's
        # structure_table, X_geo is always the site-level dbMEM eigenvectors
        # (no config — see scripts/pregea_dbmem.R). The climate/geography
        # confounding diagnostic is always computed, not switchable. dbMEM
        # itself has no user-facing knobs: always one point per site
        # (coordinate centroid, down to n=1 sample), skip-with-warning below
        # 3 sites rather than a min-sites config.
        s("Climate.predictors", "Climate predictors", "climate", "Climate predictors",
          "bio_chips", TRUE, group = "Predictors",
          help = "Same list as the GEA tab's Climate predictors — editing here updates there too. Invariant predictors (flagged by check_climate_variance, e.g. zero-variance across sites) are auto-blocked: including one would z-score to an all-NaN column and silently empty out the variance-partition table below."),
        s("Climate.Varpart.response", "Response matrix (Y)", "climate", "Response matrix (Y)",
          "select", TRUE, group = "Variance partitioning",
          choices = c("pcs", "snps"),
          help = "'pcs' = genomic PCA scores (default, scalability adaptation). 'snps' = raw LD-pruned genotypes (literal Lasky et al. 2012 methodology, impractical at WGS scale)."),
        s("Climate.Varpart.response_var_cutoff", "Response variance cutoff", "climate", "Response matrix (Y)",
          "numeric", TRUE, group = "Variance partitioning",
          min = 0, max = 1, step = 0.05,
          help = "Keep the leading genomic PCs that together explain at least this fraction of genetic variance (0.8 = 80%, the standard choice). Applies only when response = 'pcs'."),
        s("Climate.Varpart.structure_table", "Structure covariate (X)", "climate", "Explanatory tables (X)",
          "select", TRUE, group = "Variance partitioning",
          choices = c("qmatrix", "none"),
          help = "Population-structure explanatory table: 'qmatrix' = sNMF Q at sNMF.k_best, last column dropped (Q rows sum to 1, so all K columns are collinear). Uses Q rather than PCs here because Y is already PCs — an ancestry-PC covariate would make X = Y, a tautology. 'none' drops the structure fraction (2-table climate-vs-geography varpart only). Climate (X_clim) is Climate.predictors on the GEA tab; geography (X_geo) is always the site-level dbMEM eigenvectors, no config."),
        s("Climate.Varpart.response_max_pcs", "Response PCs (max)", "climate", "Response matrix (Y)",
          "numeric", FALSE, group = "Variance partitioning",
          min = 1, step = 1,
          help = "Cap on retained PCs regardless of variance cutoff"),
        s("Climate.Varpart.response_min_pcs", "Response PCs (min)", "climate", "Response matrix (Y)",
          "numeric", FALSE, group = "Variance partitioning",
          min = 1, step = 1,
          help = "Floor so varpart always has a multivariate Y"),
        s("Climate.Varpart.permutations", "Varpart permutations", "climate", "Explanatory tables (X)",
          "numeric", FALSE, group = "Variance partitioning",
          min = 99, step = 100,
          help = "anova.cca permutation count for each testable variance fraction (999 = literature/production default; lower values run faster)"),

        # ── TRAITS TAB (Phenotypic Factors, mode=traits) ──────────────────────
        # Deliberately one knob: which figures get built is decided by the data
        # (traits present) and the regime (climate correlogram only when there
        # is climate), not by config.
        s("Traits.pairs_max_factors", "Max traits in pairs plot", "traits", "Pairs plot",
          "numeric", TRUE,
          min = 2, step = 1,
          help = "Above this many traits the pairs grid becomes unreadable (and quadratically expensive), so the pipeline writes a placeholder instead of the plot. Default 8."),

        # ── PREGEA TAB ────────────────────────────────────────────────────────
        # Focused hyperparameter-choice module: "how many K / how many PCs".
        # Predictor characterization (correlations, densities, varpart, dbMEM)
        # moved to the Climate tab above. Uses the mandatory/Advanced split
        # (no `group`, unlike the old per-block layout) so Advanced/TransferGuard
        # params collapse behind an accordion instead of always showing.
        #
        # Fixed constants no longer exposed here at all (see README "Fixed
        # constants"): seed=42, LFMM/EMMAX/RDA always run together (no
        # per-block enabled switches), kinship=BN always (GAPIT uses the same
        # matrix), genomic_control=FALSE always (the GEA tab's own
        # LFMM.genomic_control param, above, is the only place this toggles),
        # FDR=0.1, RDA.min_predictors=2, lambda_tol=0.15, deflation_floor=0.90,
        # ordiR2step's Pin/R2permutations (Blanchet double-stopping-rule standard).
        s("PreGEA.predictors", "Predictors", "pregea", "Predictors",
          "bio_chips", TRUE,
          help = "Which climate predictors the LFMM-K/EMMAX-#PC/RDA ladders sweep. Falls back to Climate.predictors when left at its seeded default."),
        s("PreGEA.k_offset", "K offset", "pregea", "Ladder range",
          "numeric", TRUE,
          min = 1, max = 10, step = 1,
          help = "LFMM K sweep = sNMF.k_best +/- this (e.g. offset=2 with k_best=3 sweeps K in 1..5)"),
        s("PreGEA.n_pcs_max", "#PCs / Condition-PCs (max)", "pregea", "Ladder range",
          "numeric", TRUE,
          min = 1, max = 30, step = 1,
          help = "Shared ceiling for BOTH the EMMAX/GAPIT #PC ladder and RDA's Condition()-PC ladder (0..this). One range, one shared recommendation — most of the time the recommended #PCs and condition_pcs are equal."),
        s("PreGEA.Advanced.collinearity_r", "Collinearity |r| pre-screen", "pregea", "Advanced",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Predictor pairs above this |r| are dropped before the RDA fit (pairwise pre-screen). Complementary to VIF below, which catches multi-way collinearity a pairwise |r| check cannot see."),
        s("PreGEA.Advanced.vif_max", "VIF max (post-fit)", "pregea", "Advanced",
          "numeric", FALSE,
          min = 1, step = 0.5,
          help = "vif.cca() cutoff applied after the |r| pre-screen — a Condition()-PC value whose max VIF reaches this is excluded from the recommendation (a real gate, not just a plot line)."),
        s("PreGEA.Advanced.axis_alpha", "RDA axis significance alpha", "pregea", "Advanced",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.01,
          help = "anova.cca per-axis significance threshold"),
        s("PreGEA.Advanced.permutations", "RDA permutations", "pregea", "Advanced",
          "numeric", FALSE,
          min = 99, step = 100,
          help = "anova.cca permutation count (production default 999; lower values run faster for testing)"),
        s("PreGEA.TransferGuard.enabled", "Run transfer guard", "pregea", "Transfer guard (opt-in)",
          "checkbox", FALSE,
          help = "Re-estimates LFMM/EMMAX lambda on the FULL marker set at the recommended hyperparameters, before committing to the full GEA run. OFF by default — pulls the full-set imputation chain into this otherwise LD-pruned-only mode."),
        s("PreGEA.TransferGuard.lfmm_k", "LFMM K override", "pregea", "Transfer guard (opt-in)",
          "text", FALSE, show_if = "PreGEA.TransferGuard.enabled",
          help = "'auto' resolves from pregea_recommendations.tsv, or enter an explicit K",
          placeholder = "auto"),
        s("PreGEA.TransferGuard.emmax_n_pcs", "EMMAX #PCs override", "pregea", "Transfer guard (opt-in)",
          "text", FALSE, show_if = "PreGEA.TransferGuard.enabled",
          help = "'auto' resolves from pregea_recommendations.tsv, or enter an explicit #PCs",
          placeholder = "auto"),

        # ── GEA TAB ───────────────────────────────────────────────────────────
        s("GEA.configs", "GEA methods", "gea", "Methods",
          "method_table", TRUE,
          help = "Methods to run, with optional per-method hyperparameters (expand a row's Params to edit). Method list and hyperparameter defaults come from the pipeline's method registry (workflow/methods/gea.py) — read the preGEA ladders before overriding a default."),
        # NOTE: SNP clumping distance is no longer config-editable \u2014 it is a purely
        # interactive parameter in the GEA/GWAS filter bar (default 100 kb), kept in
        # sync with the region table/rectangles. See mod_gea.R / mod_gwas.R.
        s("GEA.wza.window_size", "WZA window size",    "gea", "WZA",
          "text",    FALSE,
          help = "Window size for Weighted-Z Analysis (shown in the WZA regime view). auto_genome_wide (default) uses genome-wide LD r\u00b2=0.2; auto_per_chromosome uses per-chr; or enter a fixed bp integer. Changing this requires re-running the GEA module.",
          placeholder = "auto_genome_wide"),
        s("GEA.wza.fallback_window_bp", "WZA fallback window (bp)", "gea", "WZA",
          "numeric", FALSE,
          min = 1000, step = 1000,
          help = "Window size (bp) used when LD decay data is unavailable (default: 10000 bp). Applies only in auto mode. Changing this requires re-running the GEA module."),
        # Gene annotation
        s("GFF.feature",   "GFF feature type", "gea", "Gene Annotation",
          "text",    FALSE,
          help = "Feature type to extract from GFF3 for gene overlaps",
          placeholder = "mRNA"),
        s("GFF.gene_name", "Gene name field",  "gea", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute used as gene display name in output tables",
          placeholder = "description"),
        s("GFF.biotype",   "Biotype field",    "gea", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute for gene biotype (used in filtering)",
          placeholder = "biotype"),
        s("GFF.go_field", "GO field in GFF",  "gea", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute containing GO term IDs. NULL disables enrichment analysis.",
          placeholder = "NULL"),
        s("GEA.promoter_length",  "Promoter length (bp)","gea","Gene Annotation",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Upstream distance from gene start counted as promoter region"),
        # Enrichment
        s("Enrichment.top_terms",        "Top GO terms",      "gea", "Enrichment",
          "numeric", FALSE,
          min = 5, max = 100, step = 5,
          help = "Number of enriched GO terms shown in dot/emap/cnetplot"),
        s("Enrichment.cnet_label",       "Cnet gene label",   "gea", "Enrichment",
          "select",  FALSE,
          choices = c("gene_id", "Name", "description"),
          help = "Which gene attribute to show as labels in the concept network plot"),

        # ── GWAS TAB ─────────────────────────────────────────────────────────
        s("GWAS.traits",           "Phenotype traits",    "gwas", "Traits",
          "pheno_chips", TRUE,
          help = "Click to toggle which phenotype traits to include in GWAS. Discovered from metadata columns 5+. Empty = all traits."),
        s("GWAS.configs",          "GWAS methods",        "gwas", "Methods",
          "method_table", TRUE,
          help = "Association methods for phenotype GWAS. Same options as GEA methods."),
        s("GWAS.missing_strategy", "Missing data strategy","gwas", "Methods",
          "select",       TRUE,
          choices = c("DROP", "MEAN", "MEDIAN"),
          help = "How to handle samples missing phenotype data. DROP removes them per trait."),
        # NOTE: SNP clumping distance is no longer config-editable \u2014 see GEA note above.
        s("GWAS.wza.window_size", "WZA window size",    "gwas", "WZA",
          "text",    FALSE,
          help = "Window size for Weighted-Z Analysis (shown in the WZA regime view). auto_genome_wide (default) uses genome-wide LD r\u00b2=0.2; auto_per_chromosome uses per-chr; or enter a fixed bp integer. Changing this requires re-running the GWAS module.",
          placeholder = "auto_genome_wide"),
        s("GWAS.wza.fallback_window_bp", "WZA fallback window (bp)", "gwas", "WZA",
          "numeric", FALSE,
          min = 1000, step = 1000,
          help = "Window size (bp) used when LD decay data is unavailable (default: 10000 bp). Applies only in auto mode. Changing this requires re-running the GWAS module."),
        s("GWAS.promoter_length",  "Promoter length (bp)","gwas", "Gene Annotation",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Promoter length (bp) for gene annotation around significant regions."),

        # ── GEAxGWAS TAB ──────────────────────────────────────────────────────
        s("GEAxGWAS.pairwise.window_size", "Pairwise window (bp)","gea_x_gwas", "Pairwise Analysis",
          "numeric", FALSE,
          min = 0, step = 100000,
          help = "Window size for scoring GEA vs GWAS pairwise region overlap"),
        s("GEAxGWAS.pairwise.min_snps",    "Min pairwise SNPs",   "gea_x_gwas", "Pairwise Analysis",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Minimum shared significant SNPs to call a pairwise overlap"),

        # ── HAPLOTYPE TAB ──────────────────────────────────────────────────────
        s("haplotype.scan.regions_source", "Regions source",    "haplotype", "Scan Setup",
          "select",  TRUE,
          choices = c("gea", "gwas", "gea_x_gwas", "all", "custom"),
          help = "Which pipeline output to pull candidate regions from"),
        s("haplotype.epsilon_selected",    "Selected epsilon",  "haplotype", "Scan Setup",
          "text",    TRUE,
          help = "Epsilon chosen after reviewing clustree. Required before running haplotype mode.",
          placeholder = "NULL"),
        s("haplotype.scan.top_regions",       "Top regions",       "haplotype", "Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Number of top regions to scan per source"),
        s("haplotype.scan.min_snps",          "Min SNPs per region","haplotype", "Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Skip regions with fewer SNPs than this threshold"),
        s("haplotype.scan.min_group_size",    "Min group size (MGmin)","haplotype","Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "DBSCAN MinPts parameter controlling cluster density"),
        s("haplotype.scan.min_haplotype_size","Min haplotype size","haplotype", "Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Minimum number of samples to call a valid haplotype group"),
        s("haplotype.scan.epsilon_range",     "Epsilon range",     "haplotype", "Scan Parameters",
          "textarea", FALSE,
          help = "Comma-separated epsilon values to scan (e.g. 0.1,0.2,0.4,0.6,0.8)",
          placeholder = "0.01,0.1,0.2,0.4,0.6,0.8,1.0"),
        s("haplotype.scan.metadata_type",     "Metadata type",     "haplotype", "Scan Parameters",
          "select",  FALSE,
          choices = c("site", "cluster", "both"),
          help = "Metadata variable to use for haplotype group annotation"),
        s("haplotype.scan.regions_file",      "Custom regions file","haplotype","Scan Parameters",
          "text",    FALSE,
          help = "Path to custom regions TSV — only used when source = 'custom'",
          placeholder = "NULL"),

        # ── MALADAPTATION TAB ──────────────────────────────────────────────────
        s("Future.ssp",    "Climate scenario (SSP)","maladaptation","Future Climate",
          "select",  TRUE,
          choices = c("126", "245", "370", "585"),
          help = "Shared Socioeconomic Pathway (126 = low emissions, 585 = high emissions)"),
        s("Future.year",   "Projection period",    "maladaptation","Future Climate",
          "select",  TRUE,
          choices = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
          help = "Decade range for future climate projections"),
        s("Future.models", "Climate models",       "maladaptation","Future Climate",
          "text",    FALSE,
          help = "Comma-separated CMIP6 model names to average. Blank = use all available.",
          placeholder = "(all models)"),
        s("Maladaptation.snp_sets",
          "SNP sets to run",    "maladaptation","SNP Sets",
          "snp_set_picker", TRUE,
          help = "Saved sets from the GEA tab. Selected sets define which adaptive-SNP sets GF runs on."),
        s("Maladaptation.methods.gradient_forest.ntree",
          "Number of trees",    "maladaptation","Gradient Forest",
          "numeric", FALSE,
          min = 100, step = 100,
          help = "Number of trees per Gradient Forest run (more = more stable, slower)"),
        s("Maladaptation.methods.gradient_forest.cor_threshold",
          "R\u00b2 threshold",  "maladaptation","Gradient Forest",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Minimum R\u00b2 to include a climate predictor in the model"),
        s("Maladaptation.methods.gradient_forest.spatial_correction",
          "Spatial correction", "maladaptation","Gradient Forest",
          "select",  FALSE,
          choices = c("with", "without", "both"),
          help = "Include forward-selected dbMEM spatial vectors as covariates (requires mode=climate to have run first). 'both' runs each SNP set as spatial AND nospatial."),
        s("Maladaptation.methods.gradient_forest.random_model",
          "Random model",       "maladaptation","Gradient Forest",
          "checkbox", FALSE,
          help = "Also build a random-SNP neutral model for comparison with adaptive model"),
        s("Maladaptation.methods.geometric_offset.k",
          "Latent factors (K)", "maladaptation","Geometric Offset",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Number of LFMM2 latent factors. Leave blank to use sNMF k_best."),
        s("Maladaptation.methods.geometric_offset.scale",
          "Scale environment",  "maladaptation","Geometric Offset",
          "checkbox", FALSE,
          help = "Standardise climate predictors before fitting LFMM2 (recommended)."),
        s("Maladaptation.methods.rda_offset.axes",
          "RDA axes (K)",       "maladaptation","RDA Offset",
          "text", FALSE,
          help = "'auto' (anova by-axis at axis_alpha, floored at 2) or an integer >= 2."),
        s("Maladaptation.methods.rda_offset.axis_alpha",
          "Axis alpha",         "maladaptation","RDA Offset",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.01,
          help = "Significance threshold for anova.cca(by=\"axis\") when axes='auto'."),
        s("Maladaptation.methods.rda_offset.permutations",
          "Permutations",       "maladaptation","RDA Offset",
          "numeric", FALSE,
          min = 99, step = 100,
          help = "Permutations for the by-axis anova.cca significance test."),
        s("Maladaptation.methods.rda_offset.condition_pcs",
          "Condition PCs",      "maladaptation","RDA Offset",
          "numeric", FALSE,
          min = 0, step = 1,
          help = "0 = canonical offset RDA (no structure conditioning, docs/rda_research.md B7). >0 is a documented deviation."),
        s("Maladaptation.methods.rda_offset.seed",
          "Random seed",        "maladaptation","RDA Offset",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Seed for the by-axis permutation test.")
    )
}

#' Get config value by dot-path key
#'
#' @param config config list from read_project_config()
#' @param key_path dot-separated path (e.g. "snmf.k_start")
#' @return the value, or NULL if not found
#' @noRd
config_get_by_path <- function(config, key_path) {
    keys <- strsplit(key_path, "\\.")[[1]]
    val <- config
    for (k in keys) {
        if (is.null(val) || !is.list(val) || !k %in% names(val)) return(NULL)
        val <- val[[k]]
    }
    val
}

#' Set config value by dot-path key (returns modified config list)
#'
#' @param config config list
#' @param key_path dot-separated path (e.g. "snmf.k_start")
#' @param value value to set
#' @return modified config list
#' @noRd
config_set_by_path <- function(config, key_path, value) {
    keys <- strsplit(key_path, "\\.")[[1]]
    if (length(keys) == 1) {
        config[[keys]] <- value
        return(config)
    }
    # Recurse: ensure intermediate list exists
    top <- keys[1]
    rest <- paste(keys[-1], collapse = ".")
    if (is.null(config[[top]]) || !is.list(config[[top]])) config[[top]] <- list()
    config[[top]] <- config_set_by_path(config[[top]], rest, value)
    config
}

#' Coerce a sidebar input value back to the appropriate R type for YAML
#'
#' @param raw_value the raw value from input[[id]]
#' @param type schema type: numeric | text | select | checkbox | textarea | method_table
#' @return R value suitable for yaml::write_yaml()
#' @noRd
input_to_config_value <- function(raw_value, type) {
    switch(type,
        "numeric" = {
            n <- suppressWarnings(as.numeric(raw_value))
            if (is.na(n) || length(n) == 0) NULL else n
        },
        "checkbox" = {
            isTRUE(raw_value)
        },
        "checkbox_invert" = {
            !isTRUE(raw_value)
        },
        "textarea" = {
            # Comma-separated string -> list of trimmed values
            if (is.null(raw_value) || !nzchar(trimws(raw_value))) return(NULL)
            parts <- trimws(strsplit(as.character(raw_value), ",")[[1]])
            parts <- parts[nzchar(parts)]
            if (length(parts) == 0) return(NULL)
            # Try numeric coercion; fall back to character list
            nums <- suppressWarnings(as.numeric(parts))
            if (!anyNA(nums)) as.list(nums) else as.list(parts)
        },
        "method_table" = raw_value,  # read-only, pass through
        {
            # text / select / default
            if (is.null(raw_value)) return(NULL)
            val <- as.character(raw_value)
            if (!nzchar(trimws(val))) NULL else val
        }
    )
}

#' Get schema entries for a specific tab
#' @noRd
schema_for_tab <- function(tab_name) {
    Filter(function(e) e$tab == tab_name, config_schema())
}

#' Split schema entries for a tab into mandatory and optional
#' @noRd
schema_split <- function(entries) {
    list(
        mandatory = Filter(function(e) isTRUE(e$mandatory), entries),
        optional  = Filter(function(e) !isTRUE(e$mandatory), entries)
    )
}
