# 1. Where were the trees sampled?

Downloads the tutorial data (the authors' VCF, the poplar gene annotation, and a metadata table
of samples and site coordinates) into `data/populus/`.

Run time on the lab: *to be measured in the tutorial run.*

:::{margin} Tips
- **Change:** nothing; the script pins the exact files and checks their checksums.
- **Trap:** `data/populus/` also holds the WorldClim download later on. Keep it between runs.
- **In Shiny:** the Home tab shows the input files a project reads.
:::

## Run it

```{code-block} text
:class: clinego-run

bash documentation/tutorial/populus/get_populus_data.sh data/populus/
```

## What should come out

- `data/populus/`: the VCF, the GFF3 and `populus_metadata.tsv`

:::{slot} Expected numbers
Samples, populations and SNPs in the downloaded VCF.
:::

## Figures

:::{slot} Figure: ch1_sites_map
The 42 sampling sites on the map.
:::

## How to read it

:::{slot} How to read it
What the sampling covers, and what the metadata columns mean.
:::

**Next:** [](02-filtering.md)
