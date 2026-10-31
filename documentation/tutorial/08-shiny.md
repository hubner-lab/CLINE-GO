# 8. How do I explore the results?

Starts the Shiny viewer over every result of the tutorial, at <http://localhost:3838>.

:::{margin} Tips
- **Change:** pick the project in the selector at the top; both tutorial projects are listed.
- **Trap:** stop the app with `docker stop clinego-app` when you are done; it holds port 3838.
- **In Shiny:** haplotypes: in a region of the GEA tab, *Run Haplotype Scan* groups the trees by
  their haplotypes there.
:::

## Run it

```{code-block} text
:class: clinego-run

|app|
```

Then open <http://localhost:3838> in a browser.

## What you will see

:::{slot} Screenshot: the GEA tab
A region selected from the Manhattan plot, with its genes and regional plot.
:::

:::{slot} Screenshot: the Maladaptation tab
The offset map with its controls.
:::

## How to read it

:::{slot} How to read it
Where each chapter's figures live in the app, and what the app computes on demand.
:::

**Back to:** [](index.md)
