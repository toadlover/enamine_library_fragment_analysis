#!/usr/bin/env python3
"""
Node-parallel Slurm adapter for enamine_recipe_feature_pipeline_with_query.py.

Version 2 changes
-----------------
* Designed for explicit partial-node allocations rather than --exclusive jobs.
* Defaults to 8 local extraction workers.
* Reports Slurm CPU and memory allocation at startup.
* Refuses to start when --local-workers exceeds the allocated CPU count.
* Preserves restart-safe manifest reuse and completed-batch skipping.

The chemistry, manifest format, feature generation, binary files, and query
behavior continue to come directly from the original pipeline.
"""

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_PIPELINE_NAME = "enamine_recipe_feature_pipeline_with_query.py"
DEFAULT_LOCAL_WORKERS = 8


def load_pipeline(path):
    pipeline_path = Path(path).expanduser().resolve()
    if not pipeline_path.is_file():
        raise SystemExit("Original pipeline not found: {}".format(pipeline_path))

    spec = importlib.util.spec_from_file_location(
        "enamine_recipe_pipeline_original", str(pipeline_path)
    )
    if spec is None or spec.loader is None:
        raise SystemExit("Could not import pipeline: {}".format(pipeline_path))

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, pipeline_path


def add_common_worker_args(parser, pipeline):
    parser.add_argument(
        "--linker-max-heavy-atoms",
        type=int,
        default=pipeline.DEFAULT_LINKER_MAX_HEAVY_ATOMS,
    )
    parser.add_argument(
        "--minhash-perm",
        type=int,
        default=pipeline.DEFAULT_MINHASH_PERM,
    )
    parser.add_argument(
        "--lsh-rows-per-band",
        type=int,
        default=pipeline.DEFAULT_LSH_ROWS_PER_BAND,
    )
    parser.add_argument("--feature-hash-seed", type=int, default=0)
    parser.add_argument("--map-include-ligand-name", action="store_true")
    parser.add_argument("--map-include-smiles", action="store_true")
    parser.add_argument("--write-recipes", action="store_true")
    parser.add_argument("--write-readable-features", action="store_true")
    parser.add_argument("--write-readable-lsh", action="store_true")
    parser.add_argument("--write-summary-csvs", action="store_true")
    parser.add_argument("--write-feature-token-map", action="store_true")
    parser.add_argument("--write-lsh-memberships-bin", action="store_true")


def append_enabled_flags(command, args):
    flags = [
        "--map-include-ligand-name",
        "--map-include-smiles",
        "--write-recipes",
        "--write-readable-features",
        "--write-readable-lsh",
        "--write-summary-csvs",
        "--write-feature-token-map",
        "--write-lsh-memberships-bin",
    ]
    for flag in flags:
        attr = flag.lstrip("-").replace("-", "_")
        if getattr(args, attr, False):
            command.append(flag)


def batch_is_complete(output_root, task_id):
    batch_name = "batch_{:06d}".format(task_id)
    stats_path = (
        Path(output_root)
        / batch_name
        / "{}.stats.json".format(batch_name)
    )

    if not stats_path.is_file():
        return False

    try:
        with open(str(stats_path), "r") as handle:
            json.load(handle)
        return True
    except Exception:
        return False


def build_worker_command(args, task_id, adapter_path):
    command = [
        args.python_executable,
        str(adapter_path),
        "--pipeline-script",
        str(args.pipeline_script),
        "worker",
        "--manifest",
        str(args.manifest),
        "--task-id",
        str(task_id),
        "--output-root",
        str(args.output_root),
        "--linker-max-heavy-atoms",
        str(args.linker_max_heavy_atoms),
        "--minhash-perm",
        str(args.minhash_perm),
        "--lsh-rows-per-band",
        str(args.lsh_rows_per_band),
        "--feature-hash-seed",
        str(args.feature_hash_seed),
    ]
    append_enabled_flags(command, args)
    return command


def run_one_worker(args, task_id, adapter_path):
    command = build_worker_command(args, task_id, adapter_path)
    started = time.time()

    print(
        "[node-worker] starting task {}: {}".format(
            task_id, " ".join(command)
        ),
        flush=True,
    )

    completed = subprocess.run(command)
    elapsed = time.time() - started

    if completed.returncode != 0:
        raise RuntimeError(
            "Task {} failed with exit code {} after {:.1f} seconds".format(
                task_id, completed.returncode, elapsed
            )
        )

    print(
        "[node-worker] completed task {} in {:.1f} seconds".format(
            task_id, elapsed
        ),
        flush=True,
    )
    return task_id


def create_or_reuse_manifest(args, pipeline):
    output_root = Path(args.output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    manifest_path = output_root / "manifest.tsv"
    count_path = output_root / "input_file_counts.tsv"

    if (
        args.reuse_manifest
        and manifest_path.is_file()
        and count_path.is_file()
    ):
        task_ids = set()
        with open(str(manifest_path), "r") as handle:
            handle.readline()
            for line in handle:
                if line.strip():
                    task_ids.add(int(line.split("\t", 1)[0]))

        if not task_ids:
            raise SystemExit(
                "Existing manifest has no tasks: {}".format(manifest_path)
            )

        print("Reusing manifest: {}".format(manifest_path), file=sys.stderr)
        return str(manifest_path), len(task_ids)

    files = pipeline.discover_bz2_files(args.input_dir, args.pattern)
    if not files:
        raise SystemExit(
            "No files matching {!r} found under {}".format(
                args.pattern, args.input_dir
            )
        )

    manifest, counts, num_tasks, total_ligands = pipeline.make_manifest(
        files,
        str(output_root),
        args.batch_size,
        whole_files=args.whole_files,
    )

    print("Manifest: {}".format(manifest), file=sys.stderr)
    print("Input counts: {}".format(counts), file=sys.stderr)
    print("Total ligands: {}".format(total_ligands), file=sys.stderr)
    print("Tasks: {}".format(num_tasks), file=sys.stderr)
    return manifest, num_tasks


def slurm_allocated_cpus():
    for variable in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        value = os.environ.get(variable)
        if value:
            try:
                return int(value)
            except ValueError:
                pass
    return None


def report_and_validate_resources(local_workers):
    allocated_cpus = slurm_allocated_cpus()

    print("Slurm resource summary")
    print("  job ID: {}".format(os.environ.get("SLURM_JOB_ID", "not under Slurm")))
    print("  node: {}".format(os.environ.get("SLURMD_NODENAME", os.uname().nodename)))
    print("  SLURM_CPUS_PER_TASK: {}".format(
        os.environ.get("SLURM_CPUS_PER_TASK", "unset")
    ))
    print("  SLURM_CPUS_ON_NODE: {}".format(
        os.environ.get("SLURM_CPUS_ON_NODE", "unset")
    ))
    print("  SLURM_MEM_PER_NODE: {}".format(
        os.environ.get("SLURM_MEM_PER_NODE", "unset")
    ))
    print("  requested local workers: {}".format(local_workers))

    if allocated_cpus is not None and local_workers > allocated_cpus:
        raise SystemExit(
            "--local-workers={} exceeds allocated CPUs={}. "
            "Reduce local workers or request more CPUs.".format(
                local_workers, allocated_cpus
            )
        )


def node_main(args, pipeline, adapter_path):
    manifest_path, num_tasks = create_or_reuse_manifest(args, pipeline)
    args.manifest = manifest_path

    all_task_ids = list(range(1, num_tasks + 1))
    if args.skip_completed:
        task_ids = [
            task_id
            for task_id in all_task_ids
            if not batch_is_complete(args.output_root, task_id)
        ]
    else:
        task_ids = all_task_ids

    if not task_ids:
        print("All {} tasks are already complete.".format(num_tasks))
        return

    local_workers = args.local_workers
    report_and_validate_resources(local_workers)
    local_workers = max(1, min(local_workers, len(task_ids)))

    print("Node-parallel execution")
    print("  manifest: {}".format(manifest_path))
    print("  total manifest tasks: {}".format(num_tasks))
    print("  tasks remaining: {}".format(len(task_ids)))
    print("  concurrent local workers: {}".format(local_workers))
    print("  output root: {}".format(args.output_root))

    failures = []
    with ThreadPoolExecutor(max_workers=local_workers) as executor:
        futures = {
            executor.submit(run_one_worker, args, task_id, adapter_path): task_id
            for task_id in task_ids
        }

        for future in as_completed(futures):
            task_id = futures[future]
            try:
                future.result()
            except Exception as exc:
                failures.append((task_id, str(exc)))
                print(
                    "[node-worker] ERROR task {}: {}".format(task_id, exc),
                    file=sys.stderr,
                    flush=True,
                )
                if args.fail_fast:
                    for pending in futures:
                        pending.cancel()
                    break

    if failures:
        print("{} task(s) failed:".format(len(failures)), file=sys.stderr)
        for task_id, error in failures:
            print("  task {}: {}".format(task_id, error), file=sys.stderr)
        raise SystemExit(1)

    print("All requested node-local tasks completed successfully.")


def worker_main(args, pipeline):
    task_id = args.task_id
    if task_id is None:
        env_value = os.environ.get("SLURM_ARRAY_TASK_ID")
        if not env_value:
            raise SystemExit(
                "worker mode requires --task-id or SLURM_ARRAY_TASK_ID"
            )
        task_id = int(env_value)

    task = pipeline.read_manifest_task(args.manifest, task_id)
    pipeline.process_manifest_task(task, args)


def query_main(args, pipeline):
    pipeline.query_batch(args)


def build_parser(pipeline, default_pipeline_path):
    parser = argparse.ArgumentParser(
        description="Node-parallel Slurm adapter for Enamine feature extraction."
    )
    parser.add_argument(
        "--pipeline-script",
        default=str(default_pipeline_path),
    )

    subparsers = parser.add_subparsers(dest="mode", required=True)

    node = subparsers.add_parser("node")
    node.add_argument("input_dir")
    node.add_argument("--output-root", required=True)
    node.add_argument("--pattern", default="*.cxsmiles.bz2")
    node.add_argument(
        "--batch-size",
        type=int,
        default=pipeline.DEFAULT_BATCH_SIZE,
    )
    node.add_argument("--whole-files", action="store_true")
    node.add_argument(
        "--local-workers",
        type=int,
        default=DEFAULT_LOCAL_WORKERS,
    )
    node.add_argument(
        "--python-executable",
        default=sys.executable,
    )
    node.add_argument("--reuse-manifest", action="store_true")
    node.add_argument("--skip-completed", action="store_true")
    node.add_argument("--fail-fast", action="store_true")
    add_common_worker_args(node, pipeline)

    worker = subparsers.add_parser("worker")
    worker.add_argument("--manifest", required=True)
    worker.add_argument("--task-id", type=int, default=None)
    worker.add_argument("--output-root", required=True)
    add_common_worker_args(worker, pipeline)

    query = subparsers.add_parser("query")
    query.add_argument("--smiles", required=True)
    query.add_argument("--batch-dir", default=None)
    query.add_argument("--batch-prefix", default=None)
    query.add_argument("--feature-hashes-bin", default=None)
    query.add_argument("--lsh-packed-bin", default=None)
    query.add_argument("--ligand-map", default=None)
    query.add_argument("--output-prefix", default=None)
    query.add_argument("--top-n", type=int, default=100)
    query.add_argument("--write-all-feature-matches", action="store_true")
    query.add_argument("--match-feature-types", default="all")
    query.add_argument("--linker-max-heavy-atoms", type=int, default=None)
    query.add_argument("--minhash-perm", type=int, default=None)
    query.add_argument("--lsh-rows-per-band", type=int, default=None)
    query.add_argument("--feature-hash-seed", type=int, default=None)

    return parser


def main():
    adapter_path = Path(__file__).resolve()
    default_pipeline_path = adapter_path.with_name(DEFAULT_PIPELINE_NAME)

    bootstrap = argparse.ArgumentParser(add_help=False)
    bootstrap.add_argument(
        "--pipeline-script", default=str(default_pipeline_path)
    )
    bootstrap_args, _ = bootstrap.parse_known_args()

    pipeline, pipeline_path = load_pipeline(bootstrap_args.pipeline_script)
    parser = build_parser(pipeline, pipeline_path)
    args = parser.parse_args()
    args.pipeline_script = str(
        Path(args.pipeline_script).expanduser().resolve()
    )

    if args.mode == "node":
        node_main(args, pipeline, adapter_path)
    elif args.mode == "worker":
        worker_main(args, pipeline)
    elif args.mode == "query":
        query_main(args, pipeline)
    else:
        parser.error("Unknown mode: {}".format(args.mode))


if __name__ == "__main__":
    main()
