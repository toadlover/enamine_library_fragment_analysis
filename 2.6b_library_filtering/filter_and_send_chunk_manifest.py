#!/usr/bin/env python3
"""Filter one Enamine 2.6B chunk into retained params plus a SMILES manifest.

Output for each chunk:
    chunk_<CHUNK>_shorthand_params.tar.gz
    chunk_<CHUNK>_ligands.csv
    chunk_<CHUNK>_checksums.sha256

Failure policy:
- Processing/input/compression failures: remove the entire workspace.
- Missing shorthand params for an otherwise retained ligand: skip that ligand,
  count it, and continue.
- Transfer or remote-verification failure: remove all extracted/intermediate
  data and retain only the final transfer-ready files in TRANSFER_RECOVERY_ROOT.
- Successful transfer: remove everything locally.

The input library retains its source superchunk/chunk layout. The remote output
uses one directory per five-digit chunk with no superchunk subdivision.
"""

import csv
import hashlib
import os
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from typing import Iterable, List, Sequence, Set, Tuple


FILTER_ROOT = Path("/pi/summer.thyme-umw/2.6b_library_filtering/ari_test")
LIBRARY_ROOT = Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion")

WORK_ROOT = Path(
    "/pi/summer.thyme-umw/2.6b_library_filtering/filtering_space_manifest"
)
TRANSFER_RECOVERY_ROOT = Path(
    "/pi/summer.thyme-umw/2.6b_library_filtering/transfer_failed_manifest"
)

REMOTE_HOST = "ari_ginsparg_umassmed@login.aicr.ai"
REMOTE_ROOT = Path(
    "/work/umassmed/thymelab_umassmed/2.6b_filtered_params_and_smiles"
)
SSH_KEY = Path.home() / ".ssh" / "id_ed25519_aicr"

TRANSFER_ATTEMPTS = 3
TRANSFER_RETRY_SECONDS = 30

LIGAND_COLUMN = 1
RETENTION_COLUMN = 2
SMILES_COLUMN = 11

# Set during main() so signal handlers can clean the active workspace.
ACTIVE_WORKSPACE = None
FINAL_OUTPUTS_READY = False


class TransferFailure(RuntimeError):
    """Raised only after final outputs exist but transfer/publication fails."""


def remove_path(path: Path) -> None:
    """Remove a file, symlink, or directory if it exists."""
    try:
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(str(path))
    except FileNotFoundError:
        pass


def signal_handler(signum, _frame) -> None:
    """Clean processing data on termination; preserve finalized outputs only."""
    global ACTIVE_WORKSPACE, FINAL_OUTPUTS_READY

    signal_name = signal.Signals(signum).name
    print("ERROR: received {}".format(signal_name), file=sys.stderr, flush=True)

    if ACTIVE_WORKSPACE is not None and ACTIVE_WORKSPACE.exists():
        if FINAL_OUTPUTS_READY:
            print(
                "Final outputs exist; attempting compact transfer-failure retention.",
                file=sys.stderr,
                flush=True,
            )
            try:
                preserve_final_outputs(ACTIVE_WORKSPACE)
            except Exception as exc:
                print(
                    "Could not preserve final outputs after signal: {}".format(exc),
                    file=sys.stderr,
                    flush=True,
                )
                remove_path(ACTIVE_WORKSPACE)
        else:
            remove_path(ACTIVE_WORKSPACE)

    os._exit(128 + signum)


def install_signal_handlers() -> None:
    for signal_number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(signal_number, signal_handler)


def run_checked(command: Sequence[str], description: str) -> None:
    """Run a command and raise a descriptive error if it fails."""
    print("RUN:", " ".join(shlex.quote(part) for part in command), flush=True)
    try:
        subprocess.run(list(command), check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "{} failed with exit code {}".format(description, exc.returncode)
        ) from exc


def log_disk_usage(path: Path, label: str) -> None:
    """Print current workspace size and filesystem free space."""
    try:
        usage = shutil.disk_usage(str(path))
        workspace_bytes = 0

        for root, _directories, filenames in os.walk(str(path)):
            for filename in filenames:
                file_path = Path(root) / filename
                try:
                    workspace_bytes += file_path.stat().st_size
                except FileNotFoundError:
                    pass

        gib = 1024.0 ** 3
        print(
            "DISK {}: workspace={:.2f} GiB, filesystem_free={:.2f} GiB".format(
                label,
                workspace_bytes / gib,
                usage.free / gib,
            ),
            flush=True,
        )
    except Exception as exc:
        print(
            "WARNING: could not measure disk usage at {}: {}".format(label, exc),
            flush=True,
        )


def safe_extract(archive: Path, destination: Path) -> None:
    """Extract directly from a source archive while rejecting path traversal."""
    destination_resolved = destination.resolve()

    with tarfile.open(str(archive), "r:gz") as tar_handle:
        for member in tar_handle.getmembers():
            member_path = (destination / member.name).resolve()
            try:
                member_path.relative_to(destination_resolved)
            except ValueError:
                raise RuntimeError(
                    "Unsafe path {!r} in archive {}".format(member.name, archive)
                )

        tar_handle.extractall(str(destination))


def create_tar_gz(output_path: Path, source_path: Path) -> None:
    """Create an archive through a temporary file and rename atomically."""
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    remove_path(temporary_path)

    with tarfile.open(str(temporary_path), "w:gz") as tar_handle:
        tar_handle.add(str(source_path), arcname=source_path.name)

    temporary_path.replace(output_path)


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def normalize_chunk(raw_chunk: str) -> str:
    if not raw_chunk.isdigit():
        raise ValueError("Chunk must be numeric; received {!r}".format(raw_chunk))

    chunk_number = int(raw_chunk)
    if chunk_number < 0 or chunk_number > 53085:
        raise ValueError("Chunk must be between 00000 and 53085")

    return "{:05d}".format(chunk_number)


def read_blacklist(path: Path) -> Set[str]:
    """Read ligand names directly from the source blacklist file."""
    if not path.is_file():
        raise FileNotFoundError("Blacklist file does not exist: {}".format(path))

    blacklist = set()

    with path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue

            ligand = row[0].strip()
            if ligand:
                blacklist.add(ligand)

    return blacklist


def find_params_directories(workspace: Path) -> List[Path]:
    directories = sorted(
        path
        for path in workspace.glob("condensed*/single_conf_params")
        if path.is_dir()
    )

    if not directories:
        raise RuntimeError(
            "No condensed*/single_conf_params directory was extracted in {}".format(
                workspace
            )
        )

    return directories


def find_param_file(
    params_directories: Iterable[Path],
    ligand: str,
) -> Path:
    filename = "{}_shorthand_params.txt".format(ligand)
    matches = [
        directory / filename
        for directory in params_directories
        if (directory / filename).is_file()
    ]

    if not matches:
        raise FileNotFoundError(filename)

    if len(matches) > 1:
        raise RuntimeError(
            "Multiple shorthand parameter files found for {}: {}".format(
                ligand,
                ", ".join(str(path) for path in matches),
            )
        )

    return matches[0]


def process_similarity_report(
    report_path: Path,
    blacklist: Set[str],
    params_directories: Sequence[Path],
    params_output: Path,
    manifest_path: Path,
) -> Tuple[int, int, int, int]:
    """Move retained params and write a two-column manifest.

    Ligands missing params are skipped rather than failing the chunk. Because
    the manifest describes the transferred parameter collection, skipped
    missing-param ligands are not written to the manifest.
    """
    retained_count = 0
    removed_count = 0
    blacklisted_count = 0
    missing_params_count = 0
    missing_preview = []

    with report_path.open("r", newline="") as report_handle, manifest_path.open(
        "w", newline=""
    ) as manifest_handle:
        reader = csv.reader(report_handle)
        writer = csv.writer(manifest_handle, lineterminator="\n")
        writer.writerow(["ligand_name", "smiles"])

        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue

            if row[0].startswith("global_"):
                continue

            required_column = max(
                LIGAND_COLUMN,
                RETENTION_COLUMN,
                SMILES_COLUMN,
            )
            if len(row) <= required_column:
                raise RuntimeError(
                    "Malformed similarity-report row {}: expected at least {} "
                    "columns, found {}".format(
                        line_number,
                        required_column + 1,
                        len(row),
                    )
                )

            ligand = row[LIGAND_COLUMN].strip()
            retention_status = row[RETENTION_COLUMN].strip().lower()
            smiles = row[SMILES_COLUMN].strip()

            if not ligand:
                raise RuntimeError(
                    "Empty ligand name in similarity-report row {}".format(line_number)
                )

            if ligand in blacklist:
                blacklisted_count += 1
                continue

            if retention_status == "removed":
                removed_count += 1
                continue

            if not smiles:
                raise RuntimeError(
                    "Empty SMILES for retained ligand {} on row {}".format(
                        ligand,
                        line_number,
                    )
                )

            try:
                source_param = find_param_file(params_directories, ligand)
            except FileNotFoundError:
                missing_params_count += 1
                if len(missing_preview) < 20:
                    missing_preview.append(ligand)
                continue

            destination_param = params_output / source_param.name
            shutil.move(str(source_param), str(destination_param))
            writer.writerow([ligand, smiles])
            retained_count += 1

    print(
        "FILTER SUMMARY: retained={}, similarity_removed={}, blacklisted={}, "
        "missing_params_skipped={}".format(
            retained_count,
            removed_count,
            blacklisted_count,
            missing_params_count,
        ),
        flush=True,
    )

    if missing_preview:
        print(
            "First missing-param ligands skipped: {}".format(
                ", ".join(missing_preview)
            ),
            flush=True,
        )

    if retained_count == 0:
        raise RuntimeError("No ligands with available params were retained")

    return (
        retained_count,
        removed_count,
        blacklisted_count,
        missing_params_count,
    )


def remove_extracted_condensed(workspace: Path) -> None:
    """Delete all expanded source data after retained params have been moved."""
    for condensed_path in workspace.glob("condensed*"):
        remove_path(condensed_path)


def final_output_paths(workspace: Path) -> List[Path]:
    """Return only fully finalized transfer-ready files in the workspace."""
    return sorted(
        path
        for path in workspace.iterdir()
        if path.is_file()
        and (
            path.name.endswith("_shorthand_params.tar.gz")
            or path.name.endswith("_ligands.csv")
            or path.name.endswith("_checksums.sha256")
        )
    )


def preserve_final_outputs(workspace: Path) -> Path:
    """Retain only finalized outputs after a transfer/publication failure."""
    outputs = final_output_paths(workspace)

    if len(outputs) != 3:
        raise RuntimeError(
            "Expected three finalized output files, found {}".format(len(outputs))
        )

    recovery_directory = (
        TRANSFER_RECOVERY_ROOT
        / workspace.parent.name
        / workspace.name
    )
    recovery_directory.parent.mkdir(parents=True, exist_ok=True)

    if recovery_directory.exists():
        remove_path(recovery_directory)

    recovery_directory.mkdir()

    for output in outputs:
        shutil.move(str(output), str(recovery_directory / output.name))

    # Delete every non-finalized artifact and the original workspace.
    remove_path(workspace)

    try:
        workspace.parent.rmdir()
    except OSError:
        pass

    print(
        "Transfer-ready files retained at {}".format(recovery_directory),
        file=sys.stderr,
        flush=True,
    )
    return recovery_directory


def transfer_and_verify(
    files: Sequence[Path],
    checksum_file: Path,
    remote_directory: Path,
    job_token: str,
) -> None:
    """Upload to remote staging, verify, and atomically publish."""
    remote_directory_text = str(remote_directory)
    staging_directory = remote_directory / ".incoming_{}".format(job_token)
    staging_text = str(staging_directory)

    ssh_base = [
        "ssh",
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=30",
        "-i",
        str(SSH_KEY),
        REMOTE_HOST,
    ]

    try:
        run_checked(
            ssh_base
            + [
                "mkdir -p {final} && rm -rf {staging} && mkdir -p {staging}".format(
                    final=shlex.quote(remote_directory_text),
                    staging=shlex.quote(staging_text),
                )
            ],
            "remote staging-directory creation",
        )

        transfer_sources = [str(path) for path in files] + [str(checksum_file)]
        rsync_command = [
            "rsync",
            "-av",
            "--partial",
            "--checksum",
            "--timeout=600",
            "-e",
            "ssh -o BatchMode=yes -o ConnectTimeout=30 -i {}".format(
                shlex.quote(str(SSH_KEY))
            ),
        ] + transfer_sources + [
            "{}:{}/".format(REMOTE_HOST, staging_text)
        ]

        last_error = None
        for attempt in range(1, TRANSFER_ATTEMPTS + 1):
            try:
                run_checked(
                    rsync_command,
                    "rsync transfer attempt {}".format(attempt),
                )
                last_error = None
                break
            except RuntimeError as exc:
                last_error = exc

                if attempt == TRANSFER_ATTEMPTS:
                    break

                sleep_seconds = TRANSFER_RETRY_SECONDS * attempt
                print(
                    "Transfer attempt {} failed; retrying in {} seconds".format(
                        attempt,
                        sleep_seconds,
                    ),
                    flush=True,
                )
                time.sleep(sleep_seconds)

        if last_error is not None:
            raise last_error

        publish_names = [path.name for path in files] + [checksum_file.name]

        verification_command = (
            "set -e; "
            "cd {staging}; "
            "sha256sum -c {checksum}; "
            "rm -f {old_files}; "
            "mv -f {files} {final}/; "
            "cd /; "
            "rmdir {staging}"
        ).format(
            staging=shlex.quote(staging_text),
            checksum=shlex.quote(checksum_file.name),
            old_files=" ".join(
                shlex.quote(str(remote_directory / name))
                for name in publish_names
            ),
            files=" ".join(shlex.quote(name) for name in publish_names),
            final=shlex.quote(remote_directory_text),
        )

        run_checked(
            ssh_base + [verification_command],
            "remote checksum verification and publication",
        )

    except Exception as exc:
        raise TransferFailure(str(exc)) from exc


def main() -> int:
    global ACTIVE_WORKSPACE, FINAL_OUTPUTS_READY

    if len(sys.argv) != 2:
        print("Usage: {} CHUNK_ID".format(sys.argv[0]), file=sys.stderr)
        return 2

    try:
        working_chunk = normalize_chunk(sys.argv[1])
    except ValueError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    working_superchunk = str(int(working_chunk[:3]))
    lsf_job_id = os.environ.get("LSB_JOBID", "manual")
    lsf_job_index = os.environ.get("LSB_JOBINDEX", "0")
    job_token = "{}_{}".format(lsf_job_id, lsf_job_index)

    workspace = WORK_ROOT / working_chunk / "job_{}".format(job_token)
    ACTIVE_WORKSPACE = workspace
    FINAL_OUTPUTS_READY = False

    print("Source superchunk: {}".format(working_superchunk), flush=True)
    print("Chunk: {}".format(working_chunk), flush=True)
    print("Workspace: {}".format(workspace), flush=True)

    if workspace.exists():
        print(
            "Removing stale workspace from an earlier failed processing run: "
            "{}".format(workspace),
            flush=True,
        )
        remove_path(workspace)

    workspace.mkdir(parents=True)

    try:
        filtered_archive = (
            FILTER_ROOT
            / working_superchunk
            / working_chunk
            / "chunk_{}_similarity_report.csv.tar.gz".format(working_chunk)
        )
        unfiltered_location = LIBRARY_ROOT / working_superchunk / working_chunk
        params_archives = sorted(
            unfiltered_location.glob("condensed_params_and_db_*.tar.gz")
        )
        blacklist_source = unfiltered_location / "blacklist_file.csv"

        if not filtered_archive.is_file():
            raise FileNotFoundError(
                "Similarity-report archive does not exist: {}".format(
                    filtered_archive
                )
            )
        if not params_archives:
            raise FileNotFoundError(
                "No condensed_params_and_db_*.tar.gz files found in {}".format(
                    unfiltered_location
                )
            )

        blacklist = read_blacklist(blacklist_source)

        print(
            "Extracting directly from source archives; no compressed inputs "
            "are copied into scratch",
            flush=True,
        )
        safe_extract(filtered_archive, workspace)
        for archive_number, archive in enumerate(params_archives, start=1):
            print(
                "Extracting params archive {}/{}: {}".format(
                    archive_number,
                    len(params_archives),
                    archive,
                ),
                flush=True,
            )
            safe_extract(archive, workspace)

        log_disk_usage(workspace, "after extraction")

        report_path = (
            workspace
            / "chunk_{}_similarity_report.csv".format(working_chunk)
        )
        if not report_path.is_file():
            raise FileNotFoundError(
                "Extracted similarity report not found: {}".format(report_path)
            )

        params_output = workspace / "shorthand_params"
        params_output.mkdir()

        manifest_path = (
            workspace
            / "chunk_{}_ligands.csv".format(working_chunk)
        )

        params_directories = find_params_directories(workspace)

        print(
            "Filtering ligands and writing name/SMILES manifest",
            flush=True,
        )
        process_similarity_report(
            report_path=report_path,
            blacklist=blacklist,
            params_directories=params_directories,
            params_output=params_output,
            manifest_path=manifest_path,
        )

        # Delete expanded source data before compression. The retained files were
        # moved within the same filesystem into params_output.
        remove_extracted_condensed(workspace)
        remove_path(report_path)
        log_disk_usage(workspace, "after deleting expanded source data")

        print("Creating retained-params archive", flush=True)
        params_archive = (
            workspace
            / "chunk_{}_shorthand_params.tar.gz".format(working_chunk)
        )
        create_tar_gz(params_archive, params_output)

        # Once the archive exists, the uncompressed retained params are no
        # longer needed.
        remove_path(params_output)

        checksum_file = (
            workspace
            / "chunk_{}_checksums.sha256".format(working_chunk)
        )
        output_files = [params_archive, manifest_path]

        with checksum_file.open("w") as handle:
            for output_file in output_files:
                handle.write(
                    "{}  {}\n".format(
                        sha256sum(output_file),
                        output_file.name,
                    )
                )

        FINAL_OUTPUTS_READY = True
        log_disk_usage(workspace, "final transfer-ready outputs")

        remote_directory = REMOTE_ROOT / working_chunk

        print("Transferring and verifying outputs", flush=True)
        transfer_and_verify(
            files=output_files,
            checksum_file=checksum_file,
            remote_directory=remote_directory,
            job_token=job_token,
        )

    except TransferFailure as exc:
        print("TRANSFER ERROR: {}".format(exc), file=sys.stderr, flush=True)
        try:
            preserve_final_outputs(workspace)
        except Exception as preserve_exc:
            print(
                "ERROR: could not retain compact final outputs: {}".format(
                    preserve_exc
                ),
                file=sys.stderr,
                flush=True,
            )
            remove_path(workspace)
        return 1

    except Exception as exc:
        print("PROCESSING ERROR: {}".format(exc), file=sys.stderr, flush=True)
        print(
            "Removing all local data for failed chunk {}".format(working_chunk),
            file=sys.stderr,
            flush=True,
        )
        remove_path(workspace)
        try:
            workspace.parent.rmdir()
        except OSError:
            pass
        return 1

    print(
        "Transfer verified; deleting all local data",
        flush=True,
    )
    remove_path(workspace)
    try:
        workspace.parent.rmdir()
    except OSError:
        pass

    print(
        "Chunk {} completed successfully".format(working_chunk),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    install_signal_handlers()
    sys.exit(main())
