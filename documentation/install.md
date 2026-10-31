# Install

CLINE-GO runs in Docker. Today the image holds the software but not the pipeline code: you clone
the repository, build the image once, and every command mounts the repository into the container
as `/pipeline`. Run all commands from the repository root.

:::{margin} Tip
A packaged image with the code inside and a single `clinego` command is planned. When it ships,
these pages switch to it.
:::

## 1. Get the code

```{code-block} text
git clone https://github.com/hubner-lab/CLINE-GO.git
cd CLINE-GO
```

## 2. Build the image

```{code-block} text
docker build -t cline-go .
```

The build compiles R and Bioconductor packages and takes a while. Build again only after the
`Dockerfile` changes; a change to the pipeline code needs no rebuild, because the code is mounted.

## 3. Run a mode

Every pipeline step is one mode, run like this (here `processing` on the bundled test data):

```{code-block} text
|run_testdata_processing|
```

- `--user $(id -u):$(id -g)` makes the output files yours, not root's.
- `-e USER=pipeline` is needed with `--user`; see the [FAQ](faq.md#the-run-stops-at-once-with-keyerror-getpwuid).
- `--name` lets you stop the run by name; see the [FAQ](faq.md#how-do-i-stop-a-run).
- `--memory=20g` caps the container; lower it on a small machine and see the
  [FAQ](faq.md#the-container-runs-out-of-memory).

The Shiny viewer starts the same way and serves on <http://localhost:3838>:

```{code-block} text
|app|
```

## Linux and Mac

The image is built for x86-64 (amd64) only.

- **Linux on x86-64**: everything runs natively.
- **Mac with Apple silicon (M1 and later)**: Docker runs the image under emulation. Most modes
  work, more slowly. EMMAX does not: it refuses to start under emulation, so a GEA or GWAS run
  that includes EMMAX stops there. Use LFMM and RDA on a Mac, or run EMMAX on a Linux machine.
