# FAQ and troubleshooting

## What does `clinego run` actually run?

One `docker run` of the image, with your current folder mounted as the container's working
directory, and inside it one `snakemake` call for the mode you named. Add `--show` to print the
full docker command without running it:

```{code-block} text
clinego run processing --config config_myproject.yaml --show
```

Anything after `--` goes to Snakemake unchanged, for example `-- -n` for a dry run that lists
the jobs a mode would run.

## The container runs out of memory

A container killed for memory exits with status 137. Three settings help:

- `--memory <size>` (for example `--memory 16g`) caps the container. Keep it below your
  machine's memory (on Docker Desktop, below the memory given to its VM).
- `--cores <N>` (default 4) is how many jobs Snakemake runs at once, and so how many are in
  memory together. Lower it on a small machine.
- `-- --workflow-profile /pipeline/workflow/profiles/lowmem` schedules jobs by a memory budget
  (10 GB in all, 5 GB per job) instead of by core count.

## The climate download is slow; can I keep it?

Yes, and you should. WorldClim rasters are downloaded once, into the config's `Input.dir` inside
your project folder, and reused by every later run that reads the same directory. Keep that
directory between runs; deleting it means downloading again. Projects whose configs sit in the
same folder share the download when their `Input.dir` names the same subfolder. `clinego` sees
only the current folder, so a project in another folder downloads again.

## Snakemake says the directory is locked

A run that was killed (out of memory, a closed terminal) leaves a lock on
`<project>_results/`. Make sure no run of that project is still going (`docker ps`), then remove
the lock and start the mode again:

```{code-block} text
clinego unlock --config config_myproject.yaml
```

## The run stops at once with `KeyError: 'getpwuid()'`

Only a hand-written `docker run` command meets this; `clinego` never does. `--user $(id -u):$(id -g)`
runs the container as a user it has no name for, and Snakemake looks the name up at startup.
Images built from this repository set `USER=pipeline` themselves; with an older image, add
`-e USER=pipeline` to your command.

## How do I stop a run?

Press Ctrl-C in the terminal that runs it: the signal reaches Snakemake, which stops its jobs
and removes its lock, and the container goes away. From another terminal, stop the container by
name; `docker ps` lists it as `clinego-<config>-<mode>`, and `--show` prints the name a command
would use:

```{code-block} text
docker ps
docker stop <name>
```

## Can I run two modes at once?

Not for the same project: modes share the project's results directory and its lock. Run them one
after another. Two different projects can run side by side.

## maladaptation stops with "No curated SNP sets"

Genomic offset is fitted on a set of candidate SNPs, and `mode=maladaptation` stops with
"No curated SNP sets with at least 3 SNPs found" when there is none. Make one in either of two
ways:

- Declare it in the config, under `GEA`, then run `mode=gea` again; it builds every declared set
  (here one named `default`, each method's own threshold, all of them joined):

  ```{code-block} yaml
  GEA:
    snp_sets:
      - name: default
  ```

- Or open the GEA tab of the app, choose the SNPs, and click "Save SNP set for maladaptation".

A set declared but not built yet stops with "Declared SNP set(s) ... not built yet": run
`mode=gea` first. `Maladaptation.snp_sets: all` uses every set with at least 3 SNPs.

## How do I update?

Get the new code, then build its image:

```{code-block} text
cd CLINE-GO
git fetch --tags
git checkout <new release tag>
./clinego build
```

(On a clone of the `main` branch, `git pull` does the first two steps.) The new image gets the
new version's tag, so the old one stays until you remove it. Your project folders and results
are not touched. They are not recomputed either: Snakemake reruns a step when its inputs or
settings change, not when the code does. To redo one step with the new code, force its rule with
`-- -R <rule>`, which also reruns every step that reads its output; the rule names are those of
the log files under `<project>_logs/`, and `-- -n` lists the jobs a mode would run. `-- -F`
forces everything: the mode and every earlier step it reads from, climate extraction included.
