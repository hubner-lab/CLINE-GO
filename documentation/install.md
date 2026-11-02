# Install

CLINE-GO runs in Docker. You clone the repository, build the image once with the `clinego`
launcher that ships in it, and from then on run every mode and the Shiny app from your own
project folder with short `clinego` commands. The pipeline code is built into the image, so
nothing from the repository is mounted at run time.

## Requirements

- **Docker** on **x86-64 (amd64) Linux**. The image is built for amd64 only; there everything
  runs natively.
- **git** and **bash**, to get the code and run the launcher.
- Disk and memory to match your data; see the
  [FAQ](faq.md#the-container-runs-out-of-memory) for a small machine.

## 1. Get the code

```{code-block} text
|clone|
cd CLINE-GO
```

This makes a `CLINE-GO` folder with the code and the `clinego` launcher, and moves you into it.
Steps 2 and 3 run there.

## 2. Build the image

```{code-block} text
|build|
```

The first build compiles R and Bioconductor packages from source and takes a long time. The
pipeline code is copied into the image last, so a rebuild after
[updating the code](faq.md#how-do-i-update) reuses every package layer and takes a minute or two.

The image is tagged `cline-go:<version>`, after the version of the code it was built from: the
commit of a clone of `main` (such as `5c1cce1`), or `1.0` for the release. Each clone's launcher
runs the image built from that clone's code.

## 3. Put `clinego` on your PATH

`PATH` is the list of folders your shell searches when you type a command. Add the `CLINE-GO`
folder to it, so that `clinego` works from any folder:

```{code-block} text
|path|
```

This lasts until you close the terminal. To keep it for every new terminal, run this once, still
in the `CLINE-GO` folder; it adds the same line, with the folder's full path, to your shell's
startup file:

```{code-block} text
|path_keep|
```

With zsh (the default shell on macOS), write `~/.zshrc` instead of `~/.bashrc`. Then check that
the launcher finds its image:

```{code-block} text
|version|
```

It prints the code version and whether the image for it is built:

```{code-block} text
clinego <version> (code at /home/you/CLINE-GO)
image   cline-go:<version>: built (sha256:…)
```

If it says `not built`, run `./clinego build` again in the `CLINE-GO` folder.

## 4. Next

The [quick check](quickstart.md) runs the bundled test dataset in a fresh folder and starts the
Shiny viewer, to show that everything works. Every `clinego` command runs in the current folder:
inputs, configs and results all live there.
