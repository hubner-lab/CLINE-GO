# 4. How much genetic variation does climate explain?

Characterizes the climate predictors (how they correlate, which do not vary between sites) and
partitions genetic variation between climate, population structure and geography.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** `Climate.predictors` picks the WorldClim variables; drop one of each highly
  correlated pair.
- **Trap:** this mode reads the climate that `structure` extracted; run chapter 3 first.
- **In Shiny:** the Climate tab.
:::

## Run it

```{code-block} text
:class: clinego-run

|run_populus_climate|
```

## What should come out

- `PopulusTutorial_results/climate/plots/correlation_heatmap.png`
- `PopulusTutorial_results/climate/plots/varpart/varpart_venn.png`
- `PopulusTutorial_results/climate/tables/varpart/variance_partition.tsv`

:::{slot} Expected numbers
The shares of genetic variation explained by climate, structure and geography, alone and shared.
:::

## Figures

:::{slot} Figure: ch4_correlation_heatmap
Correlations between the climate predictors.
:::

:::{slot} Figure: ch4_varpart_venn
Variance partitioning between climate, structure and geography.
:::

## How to read it

:::{slot} How to read it
How much of the climate signal is confounded with structure, and what that means for GEA.
:::

**Next:** [](05-gea.md)
