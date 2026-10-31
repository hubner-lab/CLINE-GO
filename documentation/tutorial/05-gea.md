# 5. Which SNPs follow climate, and do the methods agree?

Tests every SNP for association with each climate predictor using LFMM, EMMAX and RDA on the same
imputed genotypes, keeps the SNPs the methods agree on, and groups them into regions with their
genes.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** `GEA.configs` sets each method's p-value adjustment and threshold.
- **Trap:** EMMAX does not run on an Apple-silicon Mac (see the warning below).
- **In Shiny:** the GEA tab; click a SNP in the Manhattan plot to open its region.
:::

:::{warning}
On a Mac with Apple silicon, EMMAX refuses to start and this mode stops. Remove EMMAX from
`GEA.configs` there, or run this chapter on a Linux machine. See the [FAQ](../faq.md#emmax-stops-at-once-on-my-mac).
:::

## Run it

```{code-block} text
:class: clinego-run

|run_populus_gea|
```

## What should come out

- `PopulusTutorial_results/GEA/plots/manhattan/combined/manhattan_combined_K*.png`
- `PopulusTutorial_results/GEA/tables/regions_combined.tsv`: associated regions
- `PopulusTutorial_results/GEA/tables/genes_combined.tsv`: genes in those regions

:::{slot} Expected numbers
Significant SNPs per method, SNPs the methods agree on, and regions.
:::

## Figures

:::{slot} Figure: ch5_manhattan_combined
All methods on one Manhattan plot, with the agreed SNPs marked.
:::

:::{slot} Table: top regions
The top regions and their genes.
:::

## How to read it

:::{slot} How to read it
What agreement between methods adds, and how to read a region.
:::

## Not run here

:::{tip}
- **Choosing LFMM's *K* and EMMAX's number of PCs**: `mode=pregea` sweeps them on the LD-pruned
  SNPs and recommends one value per method before the full run.
- **GO enrichment**: when your GFF carries GO terms (`GFF.go_field`), the Shiny app tests each
  region for enriched functions. The poplar GFF used here has none.
- **Traits instead of climate**: when your metadata has measured traits, `mode=gwas` runs the same
  methods against them.
:::

**Next:** [](06-offset.md)
