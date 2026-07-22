# Full-node Enamine feature extraction

## What changes

Each Slurm job handles one `M/H*` or `S/H*` folder.

Inside that single allocation, the new `node` mode:

1. creates or reuses the H-folder manifest;
2. identifies incomplete manifest tasks;
3. launches several independent feature-extraction subprocesses;
4. keeps at most `--local-workers N` subprocesses active;
5. writes each task to its original `batch_XXXXXX` directory.

With a QOS limit of four running jobs, four H-folder jobs can run concurrently.
If each uses eight local workers, the practical extraction concurrency is up to
32 independent batches.

## Files

Place these together in `65b_fragment_processing`:

```text
enamine_recipe_feature_pipeline_with_query.py
enamine_recipe_feature_pipeline_with_query_slurm_node.py
```

The second file imports the first and adds node-local concurrency.

## Submit all H folders

Edit the configuration at the top of:

```text
submit_enamine_H_full_node_jobs.sh
```

Then run:

```bash
chmod +x submit_enamine_H_full_node_jobs.sh
./submit_enamine_H_full_node_jobs.sh
```

All H-folder jobs may be queued. The QOS should allow only four to run at once.

## Recommended initial worker count

Start with:

```bash
LOCAL_WORKERS=8
CPUS_PER_TASK=32
```

The worker count, not `CPUS_PER_TASK`, controls the number of simultaneous
feature-extraction processes.

Each current extraction process is single-threaded, so eight workers use roughly
eight CPUs. Requesting 32 CPUs reserves capacity for later tuning and requests
a substantial node allocation, but whether it grants the full physical node
depends on local Slurm configuration. `#SBATCH --exclusive` is the setting that
requests exclusive use of the allocated node.

## Memory tuning

`#SBATCH --mem=0` commonly requests all memory on the node. Verify this cluster
accepts that convention.

After the first jobs finish, inspect memory:

```bash
sacct -j JOBID --format=JobID,JobName,Elapsed,AllocCPUS,MaxRSS,State
```

Because the parent launches child processes, also inspect steps:

```bash
sacct -j JOBID --units=G \
    --format=JobID,JobName,Elapsed,AllocCPUS,MaxRSS,AveRSS,State
```

Increase `LOCAL_WORKERS` gradually only if memory and I/O remain safe.

## Restart behavior

The launcher passes:

```text
--reuse-manifest
--skip-completed
```

A task is considered complete when its `batch_XXXXXX.stats.json` exists and can
be parsed. Requeued or rerun H-folder jobs skip completed batches.

## Important caveat

This increases compute parallelism within each of the four scheduler jobs, but
all workers still read compressed input and write output through the shared
filesystem. Beyond some worker count, storage I/O may become the bottleneck.
