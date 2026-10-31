# 3. Is there population structure in balsam poplar?

Runs PCA and sNMF ancestry over a range of *K*, then maps the ancestry of each population. The
`structure` mode also downloads the WorldClim climate rasters.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** `sNMF.k_best`, chosen from the cross-entropy table that `prestructure` writes; set it
  before running `structure`.
- **Trap:** the first `structure` run downloads WorldClim; it takes a while and is reused after.
- **In Shiny:** the PreStructure and Structure tabs.
:::

## Run it

```{code-block} text
:class: clinego-run

|run_populus_prestructure|
|run_populus_structure|
```

## What should come out

- `PopulusTutorial_results/PreStructure/plots/pca.png`
- `PopulusTutorial_results/PreStructure/tables/cross_entropy_K*.tsv`: cross-entropy per *K* and run
- `PopulusTutorial_results/Structure/plots/piemap/`: ancestry pie charts on climate maps

:::{slot} Expected numbers
The chosen *K* and the variance explained by the first principal components.
:::

## Figures

:::{slot} Figure: ch3_pca_structure
PCA of the trees, coloured by ancestry at the chosen *K*.
:::

:::{slot} Figure: ch3_ancestry_piemap
Ancestry proportions per population on the map.
:::

## How to read it

:::{slot} How to read it
What the clusters are, and why structure must be corrected for in chapter 5.
:::

**Next:** [](04-climate.md)
