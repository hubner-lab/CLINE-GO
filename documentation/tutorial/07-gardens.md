# 7. Does genomic offset predict how trees grow?

Fitzpatrick et al. (2021) planted these populations in two common gardens, Vermont and Indian
Head. Treating each garden's climate as a "future" climate, a population's
offset to a garden should predict how poorly it grows there. This chapter repeats that test with
CLINE-GO's offsets.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** any table of site climates can serve as a "future" through `Climate.custom`.
- **Trap:** this is a second project (`PopulusGardens`) with its own config; it reads the same
  `data/populus/`.
- **In Shiny:** the Maladaptation tab of the `PopulusGardens` project.
:::

:::{in-development}
A Shiny form for custom climates is being designed; today a script builds the garden tables.
:::

## Run it

% The two garden scripts come from lane populus; the mode chain below must match the chain
% lane populus actually ran for PopulusGardens (Handoff in the site lane plan).

:::{slot} Build the garden climate tables
The command that writes the sites' present climate and each garden's climate as "future" tables.
It runs first: the chain below reads them.
:::

```{code-block} text
:class: clinego-run

|run_gardens_processing|
|run_gardens_prestructure|
|run_gardens_structure|
|run_gardens_climate|
|run_gardens_gea|
|run_gardens_maladaptation|
```

:::{slot} Compare offset with growth
The command that fits growth against offset in both gardens and compares it with the published
offsets.
:::

## What should come out

:::{slot} Expected numbers
Per garden: the association of offset with growth, offset against climate distance by AICc, and
the Spearman correlation with the published offsets.
:::

## Figures

:::{slot} Figure: ch7_offset_vs_growth
Growth against offset, one panel per garden.
:::

## How to read it

:::{slot} How to read it
Whether CLINE-GO reproduces the paper's conclusion, and what random SNP sets show.
:::

**Next:** [](08-shiny.md)
