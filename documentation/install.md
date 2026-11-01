# Install

CLINE-GO runs in Docker. You clone the repository, build the image once with the `clinego`
launcher that ships in it, and from then on run every mode and the Shiny app from your own
project folder with short `clinego` commands. The pipeline code is built into the image, so
nothing from the repository is mounted at run time.

## Requirements

- **Docker** on **x86-64 (amd64) Linux**. The image is built for amd64 only; there everything
  runs natively.
- **Mac with Apple silicon (M1 and later)**: Docker runs the image under emulation. The build is
  slow and can stop with a compiler crash (`internal compiler error: Segmentation fault`); run
  `./clinego build` again, it resumes from the packages already built. Most modes then work,
  more slowly. EMMAX does not: it refuses to start under emulation, so a GEA or GWAS run that
  includes EMMAX stops there (see the [FAQ](faq.md#emmax-stops-at-once-on-my-mac)).
- **git** and **bash**, to get the code and run the launcher.
- Disk and memory to match your data; see the
  [FAQ](faq.md#the-container-runs-out-of-memory) for a small machine.

## 1. Get the code

```{code-block} text
|clone|
cd CLINE-GO
```

## 2. Build the image

```{code-block} text
|build|
```

The first build compiles R and Bioconductor packages from source and takes a long time, longer
still under emulation. The pipeline code is copied into the image last, so a
rebuild after [updating the code](faq.md#how-do-i-update) reuses every package layer and takes a
minute or two.

The image is tagged `cline-go:<version>`, after the version of the code it was built from (`1.0`
for the release). Each clone's launcher runs the image built from that clone's code.

## 3. Put `clinego` on your PATH

From the repository folder:

```{code-block} text
|link|
```

or add the repository folder to your `PATH`. Then check that the launcher finds its image:

```{code-block} text
|version|
```

It prints the code version and whether the image for it exists. If it says there is no image,
run `./clinego build` again.

## 4. Next

The [quick check](quickstart.md) runs the bundled test dataset in a fresh folder and starts the
Shiny viewer, to show that everything works. Every `clinego` command runs in the current folder:
inputs, configs and results all live there.
