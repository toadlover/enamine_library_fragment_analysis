#!/usr/bin/env python3
"""
Slurm-native scheduler adapter for:
    enamine_recipe_feature_pipeline_with_query.py

Place this file in the same directory as the original pipeline, or pass
--pipeline-script /path/to/enamine_recipe_feature_pipeline_with_query.py.

The chemistry, feature generation, binary formats, manifest construction, and
query implementation are imported directly from the original pipeline. This
adapter replaces only the LSF-specific scheduler interface with Slurm:

    LSB_JOBINDEX       -> SLURM_ARRAY_TASK_ID
    bsub               -> sbatch
    -q                 -> --partition
    -W                 -> --time
    -M / rusage[mem=]  -> --mem
    %J_%I              -> %A_%a

Modes
-----
submit
    Discover input files, build the same manifest as the original pipeline,
    write a Slurm array-job script, and optionally submit it with sbatch.

worker
    Process one manifest task. The task ID is read from --task-id or
    SLURM_ARRAY_TASK_ID.

query
    Run the original query implementation without scheduler changes.
"""

import argparse
import importlib.util
import os
import shlex
import subprocess
import sys
from pathlib import Path


DEFAULT_PIPELINE_NAME = "enamine_recipe_feature_pipeline_with_query.py"


def load_pipeline(path):
    """Load the original pipeline as a Python module."""
    pipeline_path = Path(path).expanduser().resolve()
    if not pipeline_path.is_file():
        raise SystemExit(
            "Could not find the original pipeline at:\n"
            f"  {pipeline_path}\n"
            "Place both scripts together or provide --pipeline-script."
        )

    spec = importlib.util.spec_from_file_location("enamine_recipe_pipeline_original", str(pipeline_path))
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not load Python module from {pipeline_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, pipeline_path


def shell_join(parts):
    """Python 3.8-compatible shell command quoting."""
    return " ".join(shlex.quote(str(part)) for part in parts)


def add_common_worker_args(parser, pipeline):
    """Mirror worker options from the original pipeline."""
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


def append_enabled_worker_flags(command, args):
    for flag in [
        "--map-include-ligand-name",
        "--map-include-smiles",
        "--write-recipes",
        "--write-readable-features",
        "--write-readable-lsh",
        "--write-summary-csvs",
        "--write-feature-token-map",
        "--write-lsh-memberships-bin",
    ]:
        attribute = flag.lstrip("-").replace("-", "_")
        if getattr(args, attribute):
            command.append(flag)


def build_worker_command(args, adapter_path, manifest_path):
    command = [
        args.python_executable,
        str(adapter_path),
        "--pipeline-script",
        str(args.pipeline_script),
        "worker",
        "--manifest",
        str(manifest_path),
        "--task-id",
        "${SLURM_ARRAY_TASK_ID}",
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
    append_enabled_worker_flags(command, args)

    # Preserve expansion of the Slurm variable while quoting all other tokens.
    return " ".join(
        token if token == "${SLURM_ARRAY_TASK_ID}" else shlex.quote(str(token))
        for token in command
    )


def write_sbatch_script(args, manifest_path, num_tasks, adapter_path):
    output_root = Path(args.output_root).expanduser().resolve()
    logs_dir = output_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    array_spec = f"1-{num_tasks}"
    if args.array_throttle and args.array_throttle > 0:
        array_spec += f"%{args.array_throttle}"

    worker_command = build_worker_command(
        args=args,
        adapter_path=adapter_path,
        manifest_path=manifest_path,
    )

    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={args.job_name}",
        f"#SBATCH --array={array_spec}",
        f"#SBATCH --output={logs_dir}/%A_%a.out",
        f"#SBATCH --error={logs_dir}/%A_%a.err",
        f"#SBATCH --cpus-per-task={args.cpus_per_task}",
    ]

    if args.partition:
        lines.append(f"#SBATCH --partition={args.partition}")
    if args.walltime:
        lines.append(f"#SBATCH --time={args.walltime}")
    if args.memory:
        lines.append(f"#SBATCH --mem={args.memory}")
    if args.account:
        lines.append(f"#SBATCH --account={args.account}")
    if args.qos:
        lines.append(f"#SBATCH --qos={args.qos}")
    if args.constraint:
        lines.append(f"#SBATCH --constraint={args.constraint}")
    if args.mail_user:
        lines.append(f"#SBATCH --mail-user={args.mail_user}")
        lines.append(f"#SBATCH --mail-type={args.mail_type}")

    lines.extend([
        "",
        "set -euo pipefail",
        "",
        'echo "Job ID: ${SLURM_JOB_ID}"',
        'echo "Array task: ${SLURM_ARRAY_TASK_ID}"',
        'echo "Node: $(hostname)"',
        'echo "Started: $(date --iso-8601=seconds)"',
        "",
    ])

    if args.module_command:
        lines.append(args.module_command)
    if args.environment_command:
        lines.append(args.environment_command)

    lines.extend([
        "",
        worker_command,
        "",
        'echo "Finished: $(date --iso-8601=seconds)"',
        "",
    ])

    script_path = (
        Path(args.sbatch_script).expanduser().resolve()
        if args.sbatch_script
        else output_root / "submit_recipe_features.slurm"
    )
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("\n".join(lines))
    script_path.chmod(0o755)
    return script_path


def submit_main(args, pipeline, adapter_path):
    files = pipeline.discover_bz2_files(args.input_dir, args.pattern)
    if not files:
        raise SystemExit(
            f"No files matching {args.pattern!r} found under {args.input_dir}"
        )

    manifest_path, count_path, num_tasks, total_ligands = pipeline.make_manifest(
        files,
        args.output_root,
        args.batch_size,
        whole_files=args.whole_files,
    )

    sbatch_script = write_sbatch_script(
        args=args,
        manifest_path=manifest_path,
        num_tasks=num_tasks,
        adapter_path=adapter_path,
    )

    print(f"Manifest: {manifest_path}", file=sys.stderr)
    print(f"Input counts: {count_path}", file=sys.stderr)
    print(f"Total ligands: {total_ligands}", file=sys.stderr)
    print(f"Array tasks: {num_tasks}", file=sys.stderr)
    print(f"Slurm script: {sbatch_script}", file=sys.stderr)
    print(shell_join(["sbatch", str(sbatch_script)]))

    if args.submit:
        subprocess.run(["sbatch", str(sbatch_script)], check=True)


def worker_main(args, pipeline):
    task_id = args.task_id
    if task_id is None:
        env_value = os.environ.get("SLURM_ARRAY_TASK_ID")
        if not env_value:
            raise SystemExit(
                "worker mode needs --task-id or SLURM_ARRAY_TASK_ID"
            )
        task_id = int(env_value)

    task = pipeline.read_manifest_task(args.manifest, task_id)
    pipeline.process_manifest_task(task, args)


def query_main(args, pipeline):
    pipeline.query_batch(args)


def build_parser(pipeline, default_pipeline_path):
    parser = argparse.ArgumentParser(
        description="Slurm adapter for the Enamine recipe/feature/compact-LSH pipeline."
    )
    parser.add_argument(
        "--pipeline-script",
        default=str(default_pipeline_path),
        help="Path to the original enamine_recipe_feature_pipeline_with_query.py.",
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    submit = subparsers.add_parser(
        "submit",
        help="Create the manifest and a Slurm array script; optionally call sbatch.",
    )
    submit.add_argument("input_dir")
    submit.add_argument("--output-root", required=True)
    submit.add_argument("--pattern", default="*.bz2")
    submit.add_argument(
        "--batch-size",
        type=int,
        default=pipeline.DEFAULT_BATCH_SIZE,
    )
    submit.add_argument("--whole-files", action="store_true")

    submit.add_argument("--partition", default=None)
    submit.add_argument("--walltime", default=None, help="Slurm time, e.g. 12:00:00 or 2-00:00:00.")
    submit.add_argument("--memory", default=None, help="Slurm memory, e.g. 16G or 64000M.")
    submit.add_argument("--cpus-per-task", type=int, default=1)
    submit.add_argument("--array-throttle", type=int, default=1500)
    submit.add_argument("--job-name", default="recipe_feat")
    submit.add_argument("--account", default=None)
    submit.add_argument("--qos", default=None)
    submit.add_argument("--constraint", default=None)
    submit.add_argument("--mail-user", default=None)
    submit.add_argument(
        "--mail-type",
        default="FAIL",
        help="Used only with --mail-user; e.g. FAIL,END or ALL.",
    )
    submit.add_argument(
        "--python-executable",
        default=sys.executable,
        help="Python executable used inside each Slurm worker.",
    )
    submit.add_argument(
        "--module-command",
        default=None,
        help='Optional shell line placed before execution, e.g. "module load miniconda3".',
    )
    submit.add_argument(
        "--environment-command",
        default=None,
        help='Optional shell line, e.g. "source ~/miniconda3/etc/profile.d/conda.sh && conda activate rdkit".',
    )
    submit.add_argument(
        "--sbatch-script",
        default=None,
        help="Optional output path for the generated Slurm script.",
    )
    submit.add_argument("--submit", action="store_true")
    add_common_worker_args(submit, pipeline)

    worker = subparsers.add_parser(
        "worker",
        help="Process one manifest task, normally from a Slurm array.",
    )
    worker.add_argument("--manifest", required=True)
    worker.add_argument("--task-id", type=int, default=None)
    worker.add_argument("--output-root", required=True)
    add_common_worker_args(worker, pipeline)

    query = subparsers.add_parser(
        "query",
        help="Run the original compact-batch SMILES query mode.",
    )
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

    # Resolve --pipeline-script before constructing defaults from the imported module.
    bootstrap = argparse.ArgumentParser(add_help=False)
    bootstrap.add_argument("--pipeline-script", default=str(default_pipeline_path))
    bootstrap_args, _ = bootstrap.parse_known_args()

    pipeline, pipeline_path = load_pipeline(bootstrap_args.pipeline_script)
    parser = build_parser(pipeline, pipeline_path)
    args = parser.parse_args()

    # Normalize this value after the complete parse.
    args.pipeline_script = str(Path(args.pipeline_script).expanduser().resolve())

    if args.mode == "submit":
        submit_main(args, pipeline, adapter_path)
    elif args.mode == "worker":
        worker_main(args, pipeline)
    elif args.mode == "query":
        query_main(args, pipeline)
    else:
        parser.error(f"Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()
