# Quick check

Runs the first two modes on the small simulated dataset that ships in the image
({{qc_samples_raw}} samples, {{qc_snps_raw}} SNPs), to show that your installation works. It needs
no climate download. Run time on the lab: about {{qc_minutes_total}} minutes.

Do the [install](install.md) first: at its end, `clinego version` says the image is built.

:::{margin} Tips
- **Change:** nothing; the test config is meant to run as it is.
- **Trap:** run every command below from the same folder; `clinego` reads and writes only the
  current folder.
- **In Shiny:** the last command below starts the app; pick the `testdata` project.
:::

## Make a folder for it

Any empty folder outside `CLINE-GO` works, so that the results stay out of the code. The install
leaves you inside `CLINE-GO`, so make one in your home folder:

```{code-block} text
mkdir ~/clinego-quick && cd ~/clinego-quick
```

## Run it

Copy the test dataset into the folder, then run the two modes:

```{code-block} text
:class: clinego-run

|example|
|run_testdata_processing|
|run_testdata_prestructure|
```

`clinego example` writes `test_data/` (the VCF, metadata and gene annotation) and
`config_testdata.yaml`, and refuses to overwrite them if they are already there.

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

## Look at it in the app

Start the Shiny viewer in the same folder; it serves on <http://localhost:3838>:

```{code-block} text
:class: clinego-run

|app|
```

Stop it with Ctrl-C.

If your counts match, the pipeline works. Go on to [run it on your own data](own-data.md).
