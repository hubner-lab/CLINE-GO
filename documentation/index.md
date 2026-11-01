# CLINE-GO

CLINE-GO is a Dockerized Snakemake pipeline, with a Shiny viewer, that goes from a VCF and the
coordinates of your sampling sites to environment-associated loci and the genomic offset of each
population under future climate.

Local adaptation leaves allele frequencies that track the environment, but in most samples that
signal is entangled with population structure, and genotype-environment association (GEA) methods
that correct for structure in different ways often disagree on which loci matter. CLINE-GO runs
LFMM, EMMAX and RDA on one imputed genotype matrix and keeps the loci they agree on, within the
distance at which linkage disequilibrium decays. Those loci are grouped into regions and annotated
with genes. Genomic offset, the mismatch between a population's adaptive genotype and its future
climate, is computed three ways (Gradient Forest, geometric and RDA offset), and each is checked
against offsets from random SNP sets.

:::{figure} ../assets/clinego-metromap.png
:alt: CLINE-GO pipeline overview, drawn as a metro map
:width: 100%
:target: _images/clinego-metromap.png

The pipeline at a glance: each coloured line follows one kind of data (genotypes, climate,
candidate loci, phenotypes) through the numbered stages. Click it for the full-size image.
:::

## What you provide

::::{grid} 1 1 3 3
:gutter: 3

:::{grid-item-card} VCF
Biallelic SNPs for every sample. Filtering, LD pruning and imputation happen inside the
pipeline.
:::

:::{grid-item-card} Metadata
A TSV with the columns `site`, `sample`, `latitude`, `longitude`, in that order, plus any
measured traits.
:::

:::{grid-item-card} GFF3 (optional)
A gene annotation for the same assembly as the VCF, used to name genes near associated regions.
:::
::::

Everything else runs in one Docker image that you build once from the repository; present-day
and future climate rasters are downloaded for you.

## Where to start

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Quick check
:link: quickstart
:link-type: doc

Run the bundled test dataset in a few minutes and see that your installation works.
:::

:::{grid-item-card} Run on your own data
:link: own-data
:link-type: doc

Set up a project folder for your VCF and metadata, write its config and run the modes in order.
:::
::::

```{toctree}
:hidden:
:maxdepth: 2

install
quickstart
own-data
manual/index
faq
citing
```
