#=============================================================================
# MODULE 4C: OVERLAPPING REGIONS (GEA + GWAS COMBINED)
#=============================================================================
if GWAS_CONFIGS and GEA_RUNS:

    rule miami_plot:
      """Static Miami plot combining GEA (top) and GWAS (bottom)."""
      input:
          gea_tables = [assoc_pvalues(m) for m in GEA_CONFIGS],
          gwas_tables = [pheno_pvalues(m) for m in GWAS_CONFIGS]
      output:
          png = O['overlap_miami'],
          svg = O['overlap_miami_svg'],
          background = f"{MOD_GEAXGWAS}plots/miami_combined_K{K_BEST}_background.png",
          coords_json = f"{MOD_GEAXGWAS}plots/miami_combined_K{K_BEST}_coords.json"
      params:
          gea_str = ','.join([
              f"{method}:{adjust}:{assoc_pvalues(method)}"
              for method, adjust in GEA_CONFIGS.items()
          ]),
          gwas_str = ','.join([
              f"{method}:{adjust}:{pheno_pvalues(method)}"
              for method, adjust in GWAS_CONFIGS.items()
          ]),
          gea_preds = PREDICTORS_SELECTED,
          gwas_preds = PHENO_PREDICTORS,
          kbest = K_BEST,
          plot_dir = f"{MOD_GEAXGWAS}plots/"
      log: f"{LOGDIR}gea_x_gwas/miami_plot.log"
      # The mode's one long job, alone: scattermore's blur runs on the job's threads
      # (scattermore_threads.R, de7512), 560 s on one thread, 268 s on four (SIMDATA). The combined
      # Manhattans stay at one: at four they held up gea's and gwas's per-trait plots (p5b).
      threads: CPU
      shell:
          """
          Rscript /pipeline/scripts/plot_miami.R \
              "{params.gea_str}" "{params.gwas_str}" \
              {params.gea_preds} {params.gwas_preds} \
              {params.kbest} {params.plot_dir} > {log} 2>&1
          """

# Pairwise trait overlap — works with GEA-only, GWAS-only, or both
if GEA_RUNS or GWAS_CONFIGS:

    rule compute_pairwise_overlaps:
      """Compute all pairwise trait overlaps from method-collapsed sig SNPs."""
      input:
          selected = [f for f in [
              O.get('selected_snps')      if GEA_RUNS     else None,
              O.get('pheno_selected_snps') if GWAS_CONFIGS else None,
          ] if f is not None]
      output:
          collapsed = O['pairwise_collapsed_snps'],
          pairwise  = O['pairwise_overlap_table']
      params:
          gea_snps  = O.get('selected_snps',       'NULL') if GEA_RUNS     else 'NULL',
          gwas_snps = O.get('pheno_selected_snps',  'NULL') if GWAS_CONFIGS else 'NULL',
          window    = PAIRWISE_WINDOW_SIZE,
          min_snps  = PAIRWISE_MIN_SNPS,
          # A single-method combine_method gives SNPs their traits from that method alone.
          gea_combine  = SIGSNPS_METHOD       if GEA_RUNS     else 'All',
          gwas_combine = PHENO_COMBINE_METHOD if GWAS_CONFIGS else 'All'
      log: f"{LOGDIR}gea_x_gwas/compute_pairwise_overlaps.log"
      shell:
          """
          Rscript /pipeline/scripts/compute_pairwise_overlaps.R \
              {params.gea_snps} {params.gwas_snps} \
              {params.window} {params.min_snps} \
              {output.collapsed} {output.pairwise} \
              "{params.gea_combine}" "{params.gwas_combine}" > {log} 2>&1
          """
