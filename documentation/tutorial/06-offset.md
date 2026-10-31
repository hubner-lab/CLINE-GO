# 6. Which populations will be most out of step with future climate?

Fits Gradient Forest, geometric and RDA offset models on the associated SNPs, projects them onto
a future climate scenario, and repeats the fit on random SNP sets for comparison.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** `Future.ssp`, `Future.year` and `Future.models` choose the climate scenario.
- **Trap:** keep `Maladaptation.predictors` well below the number of sites; with more predictors
  than sites minus one, the RDA offset drops the surplus.
- **In Shiny:** the Maladaptation tab.
:::

## Run it

```{code-block} text
:class: clinego-run

|run_populus_maladaptation|
```

## What should come out

- `PopulusTutorial_results/Maladaptation/plots/gradient_forest/`: `cumulative_importance.png` and
  the offset maps
- `PopulusTutorial_results/Maladaptation/tables/gradient_forest/`: `genetic_offset_site.tsv`,
  offset per population

:::{slot} Expected numbers
The range of offsets, and the populations with the highest and lowest offset.
:::

## Figures

:::{slot} Figure: ch6_cumulative_importance
Which predictors carry the turnover in allele frequencies, and at what values.
:::

:::{slot} Figure: ch6_offset_map
Genomic offset per population under the future scenario.
:::

## How to read it

:::{slot} How to read it
What a high offset means, and how the adaptive set compares with random SNPs.
:::

**Next:** [](07-gardens.md)
