# Quick check

Runs the first two modes on the small simulated dataset that ships in `test_data/`
({{qc_samples_raw}} samples, {{qc_snps_raw}} SNPs), to show that your installation works. It needs
no climate download. Run time on the lab: about {{qc_minutes_total}} minutes.

:::{margin} Tips
- **Change:** nothing; the test config is meant to run as it is.
- **Trap:** run from the repository root, after the [install](install.md); the commands mount
  the current directory.
- **In Shiny:** start the [app](install.md#3-run-a-mode) and pick the `testdata` project.
:::

## Run it

```{code-block} text
:class: clinego-run

|run_testdata_processing|
|run_testdata_prestructure|
```

## What should come out

Results go to `testdata_results/`:

- `Processing/tables/filtering_summary.tsv`: samples and SNPs after each filter
- `PreStructure/plots/pca.png` and `PreStructure/plots/cross_entropy_K2-7.png`
- `PreStructure/tables/cross_entropy_K2-7.tsv`: cross-entropy for each *K* and run

What a run on the lab gave:

| | Samples | SNPs |
|---|---|---|
| Raw VCF | {{qc_samples_raw}} | {{qc_snps_raw}} |
| After the filters | {{qc_samples_kept}} | {{qc_snps_filtered}} |
| After LD pruning | {{qc_samples_kept}} | {{qc_snps_ld_pruned}} |

- The filters drop {{qc_samples_removed_missingness}} sample for missing data and
  {{qc_samples_removed_relatedness}} as close relatives (the test data carries planted duplicates).
- Cross-entropy is lowest at *K* = {{qc_k_min_cross_entropy}}; the config keeps
  *K* = {{qc_k_best_config}}, the number of populations the data was simulated from.
- `processing` took {{qc_seconds_processing}} s and `prestructure` {{qc_seconds_prestructure}} s.

If your counts match, the pipeline works. Go on to the [tutorial](tutorial/index.md).
