#!/usr/bin/env python3
"""Filter one Enamine 2.6B chunk, package retained ligands, and transfer it safely.

The script is designed for LSF array execution. Each invocation uses an isolated
job-specific workspace. Local data are removed only after the remote archives
have passed SHA-256 verification.
"""

import csv
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from typing import Iterable, List, Sequence

from rdkit import Chem
from rdkit.Chem import AllChem


FILTER_ROOT = Path("/pi/summer.thyme-umw/2.6b_library_filtering/ari_test")
LIBRARY_ROOT = Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion")
WORK_ROOT = Path("/pi/summer.thyme-umw/2.6b_library_filtering/filtering_space")

REMOTE_HOST = "ari_ginsparg_umassmed@login.aicr.ai"
REMOTE_ROOT = Path("/work/umassmed/thymelab_umassmed/2.6b_conformer_library")
SSH_KEY = Path.home() / ".ssh" / "id_ed25519_aicr"
TRANSFER_ATTEMPTS = 3
TRANSFER_RETRY_SECONDS = 30


def run_checked(command: Sequence[str], description: str) -> None:
    """Run a command and raise a descriptive error if it fails."""
    print("RUN:", " ".join(shlex.quote(part) for part in command), flush=True)
    try:
        subprocess.run(list(command), check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "{} failed with exit code {}".format(description, exc.returncode)
        ) from exc


def safe_extract(archive: Path, destination: Path) -> None:
    """Extract a tar archive while rejecting path traversal entries."""
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
    """Create an archive through a temporary file and rename it atomically."""
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if temporary_path.exists():
        temporary_path.unlink()

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
        raise ValueError("Chunk must be numeric, received {!r}".format(raw_chunk))

    chunk_number = int(raw_chunk)
    if chunk_number < 0 or chunk_number > 53085:
        raise ValueError("Chunk must be between 00000 and 53085")

    return "{:05d}".format(chunk_number)


def copy_required_file(source: Path, destination: Path) -> Path:
    if not source.is_file():
        raise FileNotFoundError("Required input does not exist: {}".format(source))
    target = destination / source.name
    shutil.copy2(str(source), str(target))
    return target


def read_blacklist(path: Path) -> set:
    blacklist = set()
    with path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if row and row[0].strip():
                blacklist.add(row[0].strip())
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


def find_param_file(params_directories: Iterable[Path], ligand: str) -> Path:
    filename = "{}_shorthand_params.txt".format(ligand)
    matches = [directory / filename for directory in params_directories]
    matches = [path for path in matches if path.is_file()]

    if not matches:
        raise FileNotFoundError("Missing shorthand parameter file for {}".format(ligand))
    if len(matches) > 1:
        raise RuntimeError(
            "Multiple shorthand parameter files found for {}: {}".format(
                ligand, ", ".join(str(path) for path in matches)
            )
        )
    return matches[0]


def make_sdf(smiles: str, ligand: str, output_path: Path) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES for {}: {}".format(ligand, smiles), flush=True)
        return False

    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", ligand)

    try:
        embed_status = AllChem.EmbedMolecule(
            mol,
            randomSeed=42,
            useRandomCoords=True,
        )
        if embed_status != 0:
            print("Could not generate a 3D conformer for {}".format(ligand), flush=True)
            return False

        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)

        writer = Chem.SDWriter(str(output_path))
        if writer is None:
            raise RuntimeError("Could not open SDF writer for {}".format(output_path))
        try:
            writer.write(mol)
        finally:
            writer.close()
    except Exception as exc:
        print("SDF generation failed for {}: {}".format(ligand, exc), flush=True)
        if output_path.exists():
            output_path.unlink()
        return False

    return True


def process_similarity_report(
    report_path: Path,
    blacklist: set,
    params_directories: Sequence[Path],
    shorthand_output: Path,
    sdf_output: Path,
) -> None:
    retained_count = 0
    blacklisted_count = 0
    removed_count = 0
    sdf_failure_count = 0
    missing_param_ligands = []

    with report_path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if row[0].startswith("global_"):
                continue
            if len(row) <= 11:
                raise RuntimeError(
                    "Malformed similarity-report row {}: expected at least 12 columns".format(
                        line_number
                    )
                )

            ligand = row[1].strip()
            retained = row[2].strip()
            smiles = row[11].strip()

            if ligand in blacklist:
                blacklisted_count += 1
                continue
            if retained == "removed":
                removed_count += 1
                continue

            try:
                source_param = find_param_file(params_directories, ligand)
            except FileNotFoundError:
                missing_param_ligands.append(ligand)
                continue

            destination_param = shorthand_output / source_param.name
            shutil.move(str(source_param), str(destination_param))

            sdf_path = sdf_output / "{}.sdf".format(ligand)
            if not make_sdf(smiles, ligand, sdf_path):
                sdf_failure_count += 1
                if destination_param.exists():
                    destination_param.unlink()
                continue

            retained_count += 1

    print(
        "Processing summary: retained={}, removed={}, blacklisted={}, sdf_failures={}, missing_params={}".format(
            retained_count,
            removed_count,
            blacklisted_count,
            sdf_failure_count,
            len(missing_param_ligands),
        ),
        flush=True,
    )

    if missing_param_ligands:
        preview = ", ".join(missing_param_ligands[:20])
        raise RuntimeError(
            "{} retained ligands were missing shorthand parameter files. First entries: {}".format(
                len(missing_param_ligands), preview
            )
        )
    if retained_count == 0:
        raise RuntimeError("No retained ligands were successfully produced")


def transfer_and_verify(
    files: Sequence[Path],
    checksum_file: Path,
    remote_directory: Path,
    job_token: str,
) -> None:
    """Upload into a staging directory, verify checksums, then publish atomically."""
    remote_directory_text = str(remote_directory)
    staging_directory = remote_directory / (".incoming_{}".format(job_token))
    staging_text = str(staging_directory)

    ssh_base = ["ssh", "-o", "BatchMode=yes", "-i", str(SSH_KEY), REMOTE_HOST]

    run_checked(
        ssh_base
        + [
            "mkdir -p {} && rm -rf {} && mkdir -p {}".format(
                shlex.quote(remote_directory_text),
                shlex.quote(staging_text),
                shlex.quote(staging_text),
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
        "-e",
        "ssh -o BatchMode=yes -i {}".format(shlex.quote(str(SSH_KEY))),
    ] + transfer_sources + ["{}:{}/".format(REMOTE_HOST, staging_text)]

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
                    attempt, sleep_seconds
                ),
                flush=True,
            )
            time.sleep(sleep_seconds)

    if last_error is not None:
        raise last_error

    archive_names = [path.name for path in files]
    publish_names = archive_names + [checksum_file.name]
    move_arguments = " ".join(shlex.quote(name) for name in publish_names)

    verification_command = (
        "set -e; "
        "cd {staging}; "
        "sha256sum -c {checksum}; "
        "mv -f {files} {final}/; "
        "cd /; "
        "rmdir {staging}"
    ).format(
        staging=shlex.quote(staging_text),
        checksum=shlex.quote(checksum_file.name),
        files=move_arguments,
        final=shlex.quote(remote_directory_text),
    )

    run_checked(
        ssh_base + [verification_command],
        "remote checksum verification and publication",
    )


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: {} CHUNK_ID".format(sys.argv[0]), file=sys.stderr)
        return 2

    working_chunk = normalize_chunk(sys.argv[1])
    working_superchunk = str(int(working_chunk[:3]))

    lsf_job_id = os.environ.get("LSB_JOBID", "manual")
    lsf_job_index = os.environ.get("LSB_JOBINDEX", "0")
    job_token = "{}_{}".format(lsf_job_id, lsf_job_index)

    workspace = (
        WORK_ROOT
        / working_superchunk
        / working_chunk
        / "job_{}".format(job_token)
    )

    print("Superchunk: {}".format(working_superchunk), flush=True)
    print("Chunk: {}".format(working_chunk), flush=True)
    print("Workspace: {}".format(workspace), flush=True)

    if workspace.exists():
        raise RuntimeError(
            "Workspace already exists; refusing to overwrite possible recovery data: {}".format(
                workspace
            )
        )

    workspace.mkdir(parents=True)

    try:
        filtered_archive_source = (
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

        if not params_archives:
            raise FileNotFoundError(
                "No condensed_params_and_db_*.tar.gz files found in {}".format(
                    unfiltered_location
                )
            )

        print("Copying required inputs", flush=True)
        filtered_archive = copy_required_file(filtered_archive_source, workspace)
        blacklist_file = copy_required_file(blacklist_source, workspace)
        copied_params_archives = [
            copy_required_file(path, workspace) for path in params_archives
        ]

        print("Extracting inputs", flush=True)
        safe_extract(filtered_archive, workspace)
        for archive in copied_params_archives:
            safe_extract(archive, workspace)

        for archive in [filtered_archive] + copied_params_archives:
            archive.unlink()

        report_path = workspace / "chunk_{}_similarity_report.csv".format(
            working_chunk
        )
        if not report_path.is_file():
            raise FileNotFoundError(
                "Extracted similarity report not found: {}".format(report_path)
            )

        shorthand_output = workspace / "shorthand_params"
        sdf_output = workspace / "sdfs"
        shorthand_output.mkdir()
        sdf_output.mkdir()

        blacklist = read_blacklist(blacklist_file)
        params_directories = find_params_directories(workspace)

        print("Filtering ligands and generating SDF files", flush=True)
        process_similarity_report(
            report_path,
            blacklist,
            params_directories,
            shorthand_output,
            sdf_output,
        )

        # Remove bulky extracted inputs before compression to reduce disk pressure.
        for condensed_path in workspace.glob("condensed*"):
            if condensed_path.is_dir():
                shutil.rmtree(str(condensed_path))
            else:
                condensed_path.unlink()
        report_path.unlink()
        blacklist_file.unlink()

        print("Creating output archives", flush=True)
        sdf_archive = workspace / "sdfs.tar.gz"
        params_archive = workspace / "shorthand_params.tar.gz"
        create_tar_gz(sdf_archive, sdf_output)
        create_tar_gz(params_archive, shorthand_output)

        shutil.rmtree(str(sdf_output))
        shutil.rmtree(str(shorthand_output))

        checksum_file = workspace / "checksums.sha256"
        with checksum_file.open("w") as handle:
            for archive in (sdf_archive, params_archive):
                handle.write("{}  {}\n".format(sha256sum(archive), archive.name))

        remote_directory = REMOTE_ROOT / working_superchunk / working_chunk
        print("Transferring and verifying outputs", flush=True)
        transfer_and_verify(
            [sdf_archive, params_archive],
            checksum_file,
            remote_directory,
            job_token,
        )

    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr, flush=True)
        print(
            "Workspace retained for inspection or restart: {}".format(workspace),
            file=sys.stderr,
            flush=True,
        )
        return 1

    print("Transfer verified; deleting local workspace", flush=True)
    shutil.rmtree(str(workspace))

    # Remove now-empty chunk/superchunk directories without touching nonempty ones.
    for parent in (workspace.parent, workspace.parent.parent):
        try:
            parent.rmdir()
        except OSError:
            pass

    print("Chunk {} completed successfully".format(working_chunk), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
