# FAQ and troubleshooting

## EMMAX stops at once on my Mac

The image is x86-64 only. On Apple silicon Docker emulates it, and the EMMAX binaries are not
reliable under emulation (they used to hang forever), so the pipeline refuses to start them
there and says why. LFMM, RDA and the GAPIT models are not affected. Remove EMMAX from
`GEA.configs` (or `GWAS.configs`) on a Mac, or run that mode on a Linux machine.

## The container runs out of memory

A container killed for memory exits with status 137. Two settings help:

- `--memory=<size>` on `docker run` sets the container's limit; keep it below your machine's
  memory (on Docker Desktop, below the memory given to its VM).
- `--workflow-profile workflow/profiles/lowmem` in place of `-c4` limits how many jobs run at
  once by a memory budget rather than by core count. Raise `resources.mem_mb` in that profile to
  match your `--memory` on a bigger machine.

## The climate download is slow; can I keep it?

Yes, and you should. WorldClim rasters are downloaded once, into the config's `Input.dir`, and
reused by every later run that reads the same directory. Keep that directory between runs and
across projects; deleting it means downloading again.

## Snakemake says the directory is locked

A run that was killed (out of memory, Ctrl-C, a closed terminal) leaves a lock on
`<project>_results/`. Make sure no run of that project is still going (`docker ps`), then run the
same command with `--unlock` added after `snakemake`, and start the mode again.

## The run stops at once with `KeyError: 'getpwuid()'`

`--user $(id -u):$(id -g)` runs the container as a user it has no name for, and Snakemake looks
the name up at startup. `-e USER=pipeline` gives it one; every `docker run` command on this site
already carries it.

## How do I stop a run?

By its container name: `docker stop clinego-processing` (the name is the `--name` in the command).
Pressing Ctrl-C or closing the terminal stops the Docker client, not always the container, which
can keep running and keep its lock. Check with `docker ps`, and stop it by name.

## Can I run two modes at once?

Not for the same project: modes share the project's results directory and its lock. Run them one
after another. Two different projects can run side by side.
