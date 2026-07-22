# Version 2 changes

## Resource request

The generated jobs now use:

```bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
```

They no longer use:

```bash
#SBATCH --exclusive
#SBATCH --mem=0
```

This matches the observed QOS limits:

```text
MaxJobsPU = 4
MaxTRESPU = cpu=64
```

Therefore, up to four jobs requesting 16 CPUs each can run concurrently.

## Internal extraction concurrency

Each job launches:

```text
LOCAL_WORKERS=8
```

Thus the intended maximum is:

```text
4 Slurm jobs x 8 local workers = 32 extraction processes
```

## Python safety changes

The V2 Python adapter:

- defaults to eight local workers;
- prints CPU and memory environment variables;
- refuses to run if local workers exceed the allocated CPU count;
- reuses manifests;
- skips completed batches based on valid stats JSON files.

## Resubmission

Cancel old pending jobs:

```bash
squeue -h -u "$USER" -o "%i %j" |
awk '$2 ~ /^node_[MS]_H[0-9]+$/ {print $1}' |
xargs -r scancel
```

Place the V2 Python file in the repository's `65b_fragment_processing`
directory, then run:

```bash
chmod +x submit_enamine_H_full_node_jobs_v2.sh
./submit_enamine_H_full_node_jobs_v2.sh
```

Existing valid completed batches are preserved and skipped.
