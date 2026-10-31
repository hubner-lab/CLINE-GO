% Tutorial chapter template (design step 9). Excluded from the build (conf.py exclude_patterns).
% Copy it, keep the order, keep it short: about 150 words of prose per chapter in all.
%
% - Title: a question the chapter answers, e.g. "Is there population structure in balsam poplar?"
% - First line: what the chapter does, and its run time on the lab (a tut_chN_minutes number once known).
% - Margin: up to three tips, one each of: a parameter worth changing, a trap, where it is in Shiny.
% - Run it: ONE code block, `text` lexer, class clinego-run, substitutions only (|run_populus_<mode>|).
% - What should come out: the key files (paths under PopulusTutorial_results/), then 2-3 numbers as
%   {{tut_...}} substitutions.
% - Figures: one or two {figure} blocks from tutorial/fig/<name>.png, each caption saying what to look at.
% - How to read it: 2-4 sentences.
% - A {warning} only when something breaks the run.
% - Last line: a link to the next chapter.

# <Question?>

<What this chapter does.>

Run time on the lab: <time>.

:::{margin} Tips
- **Change:** <a parameter worth changing>
- **Trap:** <a trap>
- **In Shiny:** <tab and card>
:::

## Run it

```{code-block} text
:class: clinego-run

|run_populus_<mode>|
```

## What should come out

- `PopulusTutorial_results/<Module>/<file>`

:::{slot} Expected numbers
<2-3 numbers from tutorial/verified/numbers.tsv>
:::

## Figures

:::{slot} Figure: <name>
<what the figure shows>
:::

## How to read it

:::{slot} How to read it
<2-4 sentences>
:::

**Next:** [](<next chapter>.md)
