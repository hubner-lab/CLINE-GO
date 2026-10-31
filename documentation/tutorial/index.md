# Tutorial: balsam poplar

This tutorial runs CLINE-GO on public data from balsam poplar (*Populus balsamifera*): 336 trees
from 42 populations, genotyped and published by Fitzpatrick et al. (2021). You go from the VCF to
a map of genomic offset, and in chapter 7 you test that offset the way the paper did, against how
the same populations grew in two common gardens.

Each chapter answers one question with one command, shows what should come out, and explains how
to read the figure. Run the chapters in order, from the repository root, after the
[install](../install.md).

| Chapter | Question | Mode |
|---|---|---|
| [1](01-data.md) | Where were the trees sampled? | data download |
| [2](02-filtering.md) | How much of the data survives filtering? | `processing` |
| [3](03-structure.md) | Is there population structure in balsam poplar? | `prestructure`, `structure` |
| [4](04-climate.md) | How much genetic variation does climate explain? | `climate` |
| [5](05-gea.md) | Which SNPs follow climate, and do the methods agree? | `gea` |
| [6](06-offset.md) | Which populations will be most out of step with future climate? | `maladaptation` |
| [7](07-gardens.md) | Does genomic offset predict how trees grow? | the gardens chain |
| [8](08-shiny.md) | How do I explore the results? | the Shiny app |

:::{slot} Run times
Each chapter's run time on the lab, from the tutorial run.
:::

```{toctree}
:hidden:
:maxdepth: 1

01-data
02-filtering
03-structure
04-climate
05-gea
06-offset
07-gardens
08-shiny
```
