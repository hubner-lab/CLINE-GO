# 2. How much of the data survives filtering?

Filters samples and SNPs (missingness, minor allele frequency, relatedness, heterozygosity), prunes
SNPs in linkage disequilibrium, and records every step's count.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** `Filter.maf`, `Filter.snp_miss` and `Filter.sample_miss` set the thresholds.
- **Trap:** after editing a threshold, rerun this mode; later modes read its output.
- **In Shiny:** the Processing tab.
:::

## Run it

```{code-block} text
:class: clinego-run

|run_populus_processing|
```

## What should come out

- `PopulusTutorial_results/Processing/tables/filtering_summary.tsv`: one row per filter step
- `PopulusTutorial_results/Processing/tables/metadata.tsv`

:::{slot} Expected numbers
Trees and SNPs before and after filtering.
:::

## Figures

:::{slot} Figure: ch2_filtering_attrition
How many samples and SNPs each filter removes.
:::

## How to read it

:::{slot} How to read it
Which filter removes most, and whether that is expected for this dataset.
:::

**Next:** [](03-structure.md)
