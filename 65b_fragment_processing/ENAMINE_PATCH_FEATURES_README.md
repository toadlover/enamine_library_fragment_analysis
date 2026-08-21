# Enamine 90B additive feature patcher

## What changed from the early delta-filter prototype

`enamine_patch_features.cpp` now:

- defines molecule identity from **column 1 SMILES/CXSMILES**, not ligand name/ID;
- recursively scans every `.bz2` file under the selected old and new shard roots;
- parallelizes old/new partitioning across source files with `--threads`;
- parallelizes membership comparison across disk hash buckets with `--threads`;
- deduplicates repeated novel SMILES in the new release by default (`--keep-new-duplicates` disables that);
- never rewrites completed 64.9B feature batches;
- chooses patch batch IDs starting at `max(existing batch_XXXXXX)+1`;
- creates a persistent `patch_manifest.tsv` pointing directly into checkpointed novel bucket `.cxsmiles.bz2` files;
- invokes the existing Python/RDKit `worker` mode to generate the actual feature files;
- parallelizes feature generation with `--feature-workers`;
- checkpoints source-file partitioning, comparison buckets, the patch manifest, and completed feature batches;
- safely resumes with the identical command plus `--resume`.

## Important identity assumption

The current implementation compares the literal first-column SMILES/CXSMILES text. This is fast and exact at the string level, but it assumes Enamine serializes an unchanged molecule identically in the 64.9B and 90B releases.

Before scaling, sample known old compounds and confirm that their first-column strings are unchanged in the new release. If not, add a chemistry-aware canonicalization stage before using this tool; raw text comparison would otherwise call equivalent molecules novel.

## Build

```bash
g++ -O3 -std=c++17 -pthread \
    enamine_patch_features.cpp \
    -lbz2 \
    -o enamine_patch_features
```

If bzip2 was installed in Conda:

```bash
g++ -O3 -std=c++17 -pthread \
    -I"$CONDA_PREFIX/include" \
    enamine_patch_features.cpp \
    -L"$CONDA_PREFIX/lib" \
    -Wl,-rpath,"$CONDA_PREFIX/lib" \
    -lbz2 \
    -o enamine_patch_features
```

## Example: patch M/H32

```bash
./enamine_patch_features \
  --old-root /work/umassmed/thymelab_umassmed/compressed_unmodified_64b_library/M/H32 \
  --new-root /PATH/TO/90b_library/M/H32 \
  --feature-root /work/umassmed/thymelab_umassmed/65b_library_extracted_features/features/M/H32 \
  --work-dir /scratch/$USER/enamine_90b_patch/M_H32 \
  --pipeline-script /work/umassmed/thymelab_umassmed/enamine_library_fragment_analysis/65b_fragment_processing/enamine_recipe_feature_pipeline_with_query.py \
  --python python \
  --threads 32 \
  --feature-workers 32 \
  --buckets 512 \
  --batch-size 3000000 \
  --resume
```

Use feature-generation flags/settings matching the original 64.9B run. For example, if the old run used a nondefault seed or wrote extra outputs, pass the corresponding options here.

## Checkpoint layout

Typical `--work-dir`:

```text
patch_config.txt
old_partitions/
  source_000000/
  source_000000.done
  ...
  PHASE.done
new_partitions/
  ...
novel_buckets/
  novel_bucket_000000.cxsmiles.bz2
  bucket_000000.done
  ...
  PHASE.done
patch_manifest.tsv
patch_manifest.meta
features.done
```

If Slurm kills the process, submit the same command again with `--resume`. Completed source files, comparison buckets, and feature batches are skipped.

A feature batch is considered complete only when its `batch_XXXXXX.stats.json` looks complete. Interrupted feature batches are rerun by the existing Python worker, which recreates its output files in write/truncate mode.

## Additive feature layout

If the old shard ends at:

```text
batch_000225/
```

and there are new compounds, the patch manifest begins at:

```text
batch_000226/
batch_000227/
...
```

The old batches are never changed.

## Buckets

`--buckets` is a RAM-vs-temporary-file tradeoff, not an input-batch count. Every source `.bz2` under the shard is still scanned. More buckets reduce the number of old SMILES loaded in RAM during any one membership comparison, but increase temporary file count.

Start around `512`. Increase to `1024` or `2048` if comparison buckets use too much RAM.

## Parallelism

- `--threads N`: source partitioning and bucket comparisons.
- `--feature-workers N`: simultaneous Python/RDKit feature tasks.

On a Slurm job requesting 64 CPUs, a reasonable initial maximum is:

```text
--threads 64 --feature-workers 64
```

The phases occur sequentially, so these do not imply 128 simultaneous CPUs. Storage I/O may become the limiting factor before CPU does.

## Consensus/global aggregation

This program appends compatible per-ligand feature batches and can pass `--write-summary-csvs` if that was part of the old feature run. It does **not** invent or rewrite a separate global consensus database because no global consensus file format was established in the existing pipeline. If you have a downstream global aggregation/consensus script, run it over the old + additive batches after the patch completes.

## Cleanup

By default, scratch checkpoints are retained for auditability/restart. After all feature batches are confirmed complete, `--cleanup-on-complete` removes the large partition and novel-bucket scratch directories while retaining the small config/manifest/checkpoint records.
