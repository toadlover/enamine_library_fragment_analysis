#!/usr/bin/env python3
"""
filter_enamine_chunk_similarity.py

Filter one Enamine REAL 2.6B chunk across its ten compressed SDF subchunks.

Designed for Python 3.8+ and RDKit.

Key behavior
------------
- A requested chunk number determines its superchunk:
      superchunk = chunk // 100
- The script locates split_new_named_0.sdf.tar.gz through
  split_new_named_9.sdf.tar.gz in the chunk directory.
- Missing archives are fatal by default.
- Archives are decompressed into a temporary directory under the output
  directory and cleaned automatically.
- Molecules are processed in subchunk and record order.
- Stereochemistry is ignored for both exact duplicate detection and Morgan
  fingerprints, so enantiomers and other stereochemical variants with the
  same connectivity are collapsed.
- Broader similarity filtering uses a greedy representative algorithm:
    * compare each new molecule against all representatives retained so far;
    * remove it if its maximum Tanimoto similarity is >= the cutoff;
    * otherwise retain it as a new representative.
- The algorithm is order-dependent, but avoids constructing an infeasibly
  large all-pairs distance matrix.
- Provenance is written into both the retained SDF and the CSV report.

Example
-------
python filter_enamine_chunk_similarity.py \
    26900 \
    /path/to/output/chunk_26900 \
    --similarity-cutoff 0.75 \
    --max-ligands 50000

The fixed library layout is assumed to be:
    /pi/summer.thyme-umw/enamine-REAL-2.6billion/<superchunk>/<chunk:05d>/
"""

import argparse
import csv
import json
import os
import shutil
import sys
import tarfile
import tempfile
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator


ARCHIVE_TEMPLATE = "split_new_named_{}.sdf.tar.gz"
EXPECTED_SUBCHUNKS = tuple(range(10))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter one Enamine chunk for stereo-insensitive exact duplicates "
            "and Morgan/Tanimoto-similar ligands."
        )
    )
    parser.add_argument(
        "chunk",
        type=int,
        help="Global chunk number, for example 0, 157, or 26900.",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory in which outputs and temporary decompression are created.",
    )
    parser.add_argument(
        "--library-root",
        type=Path,
        default=Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion"),
        help=(
            "Root of the Enamine REAL 2.6B library. "
            "Default: /pi/summer.thyme-umw/enamine-REAL-2.6billion"
        ),
    )
    parser.add_argument(
        "--similarity-cutoff",
        type=float,
        default=0.75,
        help="Remove a molecule when maximum Tanimoto similarity is >= this value.",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=2,
        help="Morgan fingerprint radius. Default: 2 (ECFP4-like).",
    )
    parser.add_argument(
        "--fp-size",
        type=int,
        default=2048,
        help="Morgan fingerprint bit length. Default: 2048.",
    )
    parser.add_argument(
        "--max-ligands",
        type=int,
        default=50000,
        help=(
            "Safety ceiling for successfully parsed ligands. The run aborts if "
            "this is exceeded. Set to 0 for no ceiling. Default: 50000."
        ),
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Print a progress line after this many input records. Default: 1000.",
    )
    parser.add_argument(
        "--allow-missing-subchunks",
        action="store_true",
        help="Continue when one or more of the ten expected archives are missing.",
    )
    parser.add_argument(
        "--keep-salts",
        action="store_true",
        help=(
            "Compare and write complete multicomponent records. By default only "
            "the largest covalent fragment is retained."
        ),
    )
    parser.add_argument(
        "--name-property",
        default="_Name",
        help="SDF property to use as the ligand name. Default: _Name.",
    )
    parser.add_argument(
        "--quiet-rdkit",
        action="store_true",
        help="Suppress RDKit warning and error log messages.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print("[{}] {}".format(timestamp, message), flush=True)


def superchunk_for_chunk(chunk: int) -> int:
    if chunk < 0:
        raise ValueError("Chunk number cannot be negative.")
    return chunk // 100


def resolve_chunk_directory(library_root: Path, chunk: int) -> Path:
    """
    Resolve the fixed Enamine hierarchy:

        <library_root>/<superchunk>/<chunk:05d>/

    where:
        superchunk = chunk // 100
    """
    root = library_root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(
            "Library root does not exist or is not a directory: {}".format(root)
        )

    superchunk = superchunk_for_chunk(chunk)
    superchunk_dir = root / str(superchunk)
    chunk_dir = superchunk_dir / "{:05d}".format(chunk)

    if not superchunk_dir.is_dir():
        raise FileNotFoundError(
            "Expected superchunk directory does not exist: {}".format(
                superchunk_dir
            )
        )

    if not chunk_dir.is_dir():
        raise FileNotFoundError(
            "Expected chunk directory does not exist: {}".format(chunk_dir)
        )

    return chunk_dir.resolve()

def collect_archives(
    chunk_dir: Path,
    allow_missing: bool,
) -> List[Tuple[int, Path]]:
    archives = []
    missing = []

    for subchunk in EXPECTED_SUBCHUNKS:
        archive = chunk_dir / ARCHIVE_TEMPLATE.format(subchunk)
        if archive.is_file():
            archives.append((subchunk, archive.resolve()))
        else:
            missing.append(archive)

    if missing and not allow_missing:
        raise FileNotFoundError(
            "Missing {} expected subchunk archive(s):\n  {}".format(
                len(missing), "\n  ".join(str(path) for path in missing)
            )
        )
    if not archives:
        raise FileNotFoundError(
            "No expected subchunk archives were found in {}".format(chunk_dir)
        )

    if missing:
        log(
            "WARNING: continuing without missing subchunks: {}".format(
                ", ".join(path.name for path in missing)
            )
        )

    return archives


def safe_sdf_members(archive: tarfile.TarFile) -> List[tarfile.TarInfo]:
    members = []
    for member in archive.getmembers():
        if not member.isfile():
            continue
        name = Path(member.name)
        lower = name.name.lower()
        if lower.endswith(".sdf") or lower.endswith(".sd"):
            members.append(member)
    return members


def decompress_archive_sdfs(
    archive_path: Path,
    destination: Path,
    subchunk: int,
) -> List[Path]:
    """
    Decompress only SDF members without using tar.extract(), preventing path
    traversal and ignoring unrelated files.
    """
    extracted = []
    with tarfile.open(str(archive_path), mode="r:gz") as archive:
        members = safe_sdf_members(archive)
        if not members:
            raise RuntimeError(
                "Archive contains no .sdf or .sd member: {}".format(archive_path)
            )

        for member_number, member in enumerate(members):
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError(
                    "Could not read {} from {}".format(member.name, archive_path)
                )
            output_name = "subchunk_{:01d}_member_{:03d}_{}".format(
                subchunk, member_number, Path(member.name).name
            )
            output_path = destination / output_name
            with source, output_path.open("wb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
            extracted.append(output_path)

    return extracted


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        return Chem.Mol(mol)
    return max(
        fragments,
        key=lambda fragment: (
            fragment.GetNumHeavyAtoms(),
            fragment.GetNumAtoms(),
            Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=False),
        ),
    )


def stereo_insensitive_key(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    Chem.RemoveStereochemistry(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def molecule_name(
    mol: Chem.Mol,
    global_index: int,
    name_property: str,
) -> str:
    if mol.HasProp(name_property):
        value = mol.GetProp(name_property).strip()
        if value:
            return value
    return "chunk_molecule_{:07d}".format(global_index + 1)


def report_fieldnames() -> List[str]:
    return [
        "global_input_index",
        "name",
        "status",
        "reason",
        "superchunk",
        "chunk",
        "subchunk",
        "archive_name",
        "archive_path",
        "sdf_member_name",
        "record_index_in_member",
        "stereo_insensitive_smiles",
        "representative_global_index",
        "representative_name",
        "representative_superchunk",
        "representative_chunk",
        "representative_subchunk",
        "representative_archive_name",
        "representative_record_index_in_member",
        "similarity_to_representative",
    ]


def make_provenance(
    global_index: int,
    name: str,
    superchunk: int,
    chunk: int,
    subchunk: int,
    archive_path: Path,
    sdf_member_path: Path,
    member_record_index: int,
) -> Dict[str, object]:
    return {
        "global_input_index": global_index,
        "name": name,
        "superchunk": superchunk,
        "chunk": chunk,
        "subchunk": subchunk,
        "archive_name": archive_path.name,
        "archive_path": str(archive_path),
        "sdf_member_name": sdf_member_path.name,
        "record_index_in_member": member_record_index,
    }


def representative_fields(rep: Dict[str, object]) -> Dict[str, object]:
    return {
        "representative_global_index": rep["global_input_index"],
        "representative_name": rep["name"],
        "representative_superchunk": rep["superchunk"],
        "representative_chunk": rep["chunk"],
        "representative_subchunk": rep["subchunk"],
        "representative_archive_name": rep["archive_name"],
        "representative_record_index_in_member": rep["record_index_in_member"],
    }


def write_report_row(
    writer: csv.DictWriter,
    provenance: Dict[str, object],
    status: str,
    reason: str,
    smiles_key: str,
    representative: Optional[Dict[str, object]],
    similarity: Optional[float],
) -> None:
    row = dict(provenance)
    row.update(
        {
            "status": status,
            "reason": reason,
            "stereo_insensitive_smiles": smiles_key,
            "representative_global_index": "",
            "representative_name": "",
            "representative_superchunk": "",
            "representative_chunk": "",
            "representative_subchunk": "",
            "representative_archive_name": "",
            "representative_record_index_in_member": "",
            "similarity_to_representative": (
                "" if similarity is None else "{:.6f}".format(similarity)
            ),
        }
    )
    if representative is not None:
        row.update(representative_fields(representative))
    writer.writerow(row)


def annotate_retained_molecule(
    mol: Chem.Mol,
    provenance: Dict[str, object],
    similarity_cutoff: float,
    representative_number: int,
) -> None:
    mol.SetProp("_Name", str(provenance["name"]))
    mol.SetIntProp("_SourceSuperchunk", int(provenance["superchunk"]))
    mol.SetIntProp("_SourceChunk", int(provenance["chunk"]))
    mol.SetIntProp("_SourceSubchunk", int(provenance["subchunk"]))
    mol.SetProp("_SourceArchive", str(provenance["archive_name"]))
    mol.SetProp("_SourceArchivePath", str(provenance["archive_path"]))
    mol.SetProp("_SourceSDFMember", str(provenance["sdf_member_name"]))
    mol.SetIntProp(
        "_SourceRecordIndexInMember",
        int(provenance["record_index_in_member"]),
    )
    mol.SetIntProp("_ChunkGlobalInputIndex", int(provenance["global_input_index"]))
    mol.SetIntProp("_SimilarityRepresentativeNumber", representative_number)
    mol.SetDoubleProp("_SimilarityCutoff", similarity_cutoff)
    mol.SetProp("_SimilarityMethod", "greedy_morgan_tanimoto_stereo_insensitive")


def validate_args(args: argparse.Namespace) -> None:
    if args.chunk < 0:
        raise ValueError("Chunk must be nonnegative.")
    if not 0.0 < args.similarity_cutoff <= 1.0:
        raise ValueError("--similarity-cutoff must be in (0, 1].")
    if args.radius < 1:
        raise ValueError("--radius must be at least 1.")
    if args.fp_size < 128:
        raise ValueError("--fp-size must be at least 128.")
    if args.max_ligands < 0:
        raise ValueError("--max-ligands cannot be negative.")
    if args.progress_every < 0:
        raise ValueError("--progress-every cannot be negative.")


def main() -> int:
    args = parse_args()

    try:
        validate_args(args)
    except ValueError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    if args.quiet_rdkit:
        RDLogger.DisableLog("rdApp.*")

    start_time = time.time()
    superchunk = superchunk_for_chunk(args.chunk)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = args.output_dir.resolve()

    output_sdf = output_dir / "chunk_{:05d}_filtered.sdf".format(args.chunk)
    output_csv = output_dir / "chunk_{:05d}_similarity_report.csv".format(args.chunk)
    output_summary = output_dir / "chunk_{:05d}_summary.json".format(args.chunk)

    partial_sdf = output_dir / (output_sdf.name + ".partial")
    partial_csv = output_dir / (output_csv.name + ".partial")
    partial_summary = output_dir / (output_summary.name + ".partial")

    for partial in (partial_sdf, partial_csv, partial_summary):
        if partial.exists():
            partial.unlink()

    try:
        chunk_dir = resolve_chunk_directory(args.library_root, args.chunk)
        archives = collect_archives(chunk_dir, args.allow_missing_subchunks)
    except (FileNotFoundError, RuntimeError, OSError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2

    log("Chunk: {} | superchunk: {}".format(args.chunk, superchunk))
    log("Chunk directory: {}".format(chunk_dir))
    log(
        "Found {}/10 expected subchunk archives.".format(len(archives))
    )
    log(
        "Similarity cutoff: {:.3f} | Morgan radius: {} | fingerprint bits: {}".format(
            args.similarity_cutoff, args.radius, args.fp_size
        )
    )

    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=args.radius,
        fpSize=args.fp_size,
        includeChirality=False,
    )

    retained_fingerprints: List[object] = []
    retained_provenance: List[Dict[str, object]] = []
    # Map every exact stereo-insensitive connectivity key encountered to the
    # final retained representative, not merely to the first raw occurrence.
    exact_key_to_representative: Dict[str, int] = {}

    total_records = 0
    parsed_records = 0
    unreadable_records = 0
    fragment_failures = 0
    exact_duplicates_removed = 0
    similarity_removed = 0
    retained_count = 0
    comparisons_performed = 0

    try:
        with tempfile.TemporaryDirectory(
            prefix=".chunk_{:05d}_work_".format(args.chunk),
            dir=str(output_dir),
        ) as temp_name:
            temp_dir = Path(temp_name)
            log("Temporary extraction directory: {}".format(temp_dir))

            sdf_writer = Chem.SDWriter(str(partial_sdf))
            if sdf_writer is None:
                raise RuntimeError(
                    "Could not open temporary SDF output: {}".format(partial_sdf)
                )

            try:
                with partial_csv.open("w", newline="", encoding="utf-8") as report_handle:
                    report_writer = csv.DictWriter(
                        report_handle, fieldnames=report_fieldnames()
                    )
                    report_writer.writeheader()

                    for archive_position, (subchunk, archive_path) in enumerate(
                        archives, start=1
                    ):
                        log(
                            "Subchunk {} ({}/{}): decompressing {}".format(
                                subchunk, archive_position, len(archives), archive_path.name
                            )
                        )
                        extracted_sdfs = decompress_archive_sdfs(
                            archive_path, temp_dir, subchunk
                        )
                        log(
                            "Subchunk {}: processing {} extracted SDF member(s).".format(
                                subchunk, len(extracted_sdfs)
                            )
                        )

                        for sdf_member_path in extracted_sdfs:
                            supplier = Chem.SDMolSupplier(
                                str(sdf_member_path),
                                sanitize=True,
                                removeHs=True,
                                strictParsing=False,
                            )

                            for member_record_index, mol in enumerate(supplier):
                                global_index = total_records
                                total_records += 1

                                fallback_name = "unreadable_{:07d}".format(global_index + 1)
                                base_provenance = make_provenance(
                                    global_index=global_index,
                                    name=fallback_name,
                                    superchunk=superchunk,
                                    chunk=args.chunk,
                                    subchunk=subchunk,
                                    archive_path=archive_path,
                                    sdf_member_path=sdf_member_path,
                                    member_record_index=member_record_index,
                                )

                                if mol is None:
                                    unreadable_records += 1
                                    write_report_row(
                                        report_writer,
                                        base_provenance,
                                        status="removed",
                                        reason="SDF parsing or sanitization failed",
                                        smiles_key="",
                                        representative=None,
                                        similarity=None,
                                    )
                                    continue

                                parsed_records += 1
                                if (
                                    args.max_ligands > 0
                                    and parsed_records > args.max_ligands
                                ):
                                    raise RuntimeError(
                                        "Parsed ligand count exceeded --max-ligands {}. "
                                        "No final outputs were installed. Increase the "
                                        "ceiling or use 0 to disable it.".format(
                                            args.max_ligands
                                        )
                                    )

                                name = molecule_name(
                                    mol, global_index, args.name_property
                                )
                                provenance = make_provenance(
                                    global_index=global_index,
                                    name=name,
                                    superchunk=superchunk,
                                    chunk=args.chunk,
                                    subchunk=subchunk,
                                    archive_path=archive_path,
                                    sdf_member_path=sdf_member_path,
                                    member_record_index=member_record_index,
                                )

                                working = Chem.Mol(mol)
                                if not args.keep_salts:
                                    try:
                                        working = largest_fragment(working)
                                    except Exception as exc:
                                        fragment_failures += 1
                                        write_report_row(
                                            report_writer,
                                            provenance,
                                            status="removed",
                                            reason="fragment processing failed: {}".format(
                                                exc
                                            ),
                                            smiles_key="",
                                            representative=None,
                                            similarity=None,
                                        )
                                        continue

                                try:
                                    smiles_key = stereo_insensitive_key(working)
                                    fingerprint = fingerprint_generator.GetFingerprint(
                                        working
                                    )
                                except Exception as exc:
                                    fragment_failures += 1
                                    write_report_row(
                                        report_writer,
                                        provenance,
                                        status="removed",
                                        reason="fingerprint/key generation failed: {}".format(
                                            exc
                                        ),
                                        smiles_key="",
                                        representative=None,
                                        similarity=None,
                                    )
                                    continue

                                exact_rep_index = exact_key_to_representative.get(
                                    smiles_key
                                )
                                if exact_rep_index is not None:
                                    representative = retained_provenance[
                                        exact_rep_index
                                    ]
                                    exact_duplicates_removed += 1
                                    write_report_row(
                                        report_writer,
                                        provenance,
                                        status="removed",
                                        reason=(
                                            "exact duplicate after removing stereochemistry"
                                        ),
                                        smiles_key=smiles_key,
                                        representative=representative,
                                        similarity=1.0,
                                    )
                                else:
                                    representative_position = None
                                    best_similarity = -1.0

                                    if retained_fingerprints:
                                        similarities = (
                                            DataStructs.BulkTanimotoSimilarity(
                                                fingerprint,
                                                retained_fingerprints,
                                            )
                                        )
                                        comparisons_performed += len(
                                            retained_fingerprints
                                        )
                                        best_similarity = max(similarities)
                                        if best_similarity >= args.similarity_cutoff:
                                            representative_position = similarities.index(
                                                best_similarity
                                            )

                                    if representative_position is not None:
                                        representative = retained_provenance[
                                            representative_position
                                        ]
                                        similarity_removed += 1
                                        # Future exact duplicates of this molecule should
                                        # point directly to the retained representative.
                                        exact_key_to_representative[
                                            smiles_key
                                        ] = representative_position
                                        write_report_row(
                                            report_writer,
                                            provenance,
                                            status="removed",
                                            reason=(
                                                "Morgan/Tanimoto similarity >= {:.3f}"
                                            ).format(args.similarity_cutoff),
                                            smiles_key=smiles_key,
                                            representative=representative,
                                            similarity=best_similarity,
                                        )
                                    else:
                                        representative_position = len(
                                            retained_fingerprints
                                        )
                                        retained_fingerprints.append(fingerprint)
                                        retained_provenance.append(provenance)
                                        exact_key_to_representative[
                                            smiles_key
                                        ] = representative_position
                                        retained_count += 1

                                        annotate_retained_molecule(
                                            working,
                                            provenance,
                                            args.similarity_cutoff,
                                            retained_count,
                                        )
                                        sdf_writer.write(working)
                                        write_report_row(
                                            report_writer,
                                            provenance,
                                            status="retained",
                                            reason="greedy similarity representative",
                                            smiles_key=smiles_key,
                                            representative=provenance,
                                            similarity=1.0,
                                        )

                                if (
                                    args.progress_every > 0
                                    and total_records % args.progress_every == 0
                                ):
                                    elapsed = max(time.time() - start_time, 1e-9)
                                    log(
                                        "Processed {:,} records | parsed {:,} | "
                                        "retained {:,} | exact removed {:,} | "
                                        "similarity removed {:,} | {:.1f} records/s".format(
                                            total_records,
                                            parsed_records,
                                            retained_count,
                                            exact_duplicates_removed,
                                            similarity_removed,
                                            total_records / elapsed,
                                        )
                                    )

                            # Delete each decompressed member as soon as processed.
                            try:
                                sdf_member_path.unlink()
                            except OSError:
                                pass

            finally:
                sdf_writer.close()

        elapsed_seconds = time.time() - start_time
        summary = {
            "library_root": str(args.library_root.resolve()),
            "chunk_directory": str(chunk_dir),
            "superchunk": superchunk,
            "chunk": args.chunk,
            "subchunks_processed": [subchunk for subchunk, _ in archives],
            "archives_processed": [str(path) for _, path in archives],
            "similarity_method": "greedy_morgan_tanimoto_stereo_insensitive",
            "similarity_cutoff": args.similarity_cutoff,
            "morgan_radius": args.radius,
            "fingerprint_bits": args.fp_size,
            "keep_salts": args.keep_salts,
            "total_sdf_records": total_records,
            "successfully_parsed_ligands": parsed_records,
            "unreadable_records": unreadable_records,
            "fragment_or_fingerprint_failures": fragment_failures,
            "exact_stereo_insensitive_duplicates_removed": exact_duplicates_removed,
            "similarity_removed": similarity_removed,
            "retained_ligands": retained_count,
            "fingerprint_comparisons_performed": comparisons_performed,
            "elapsed_seconds": elapsed_seconds,
            "records_per_second": (
                total_records / elapsed_seconds if elapsed_seconds > 0 else None
            ),
            "output_sdf": str(output_sdf),
            "output_report": str(output_csv),
        }
        with partial_summary.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        # Install outputs atomically only after a successful run.
        os.replace(str(partial_sdf), str(output_sdf))
        os.replace(str(partial_csv), str(output_csv))
        os.replace(str(partial_summary), str(output_summary))

    except Exception as exc:
        for partial in (partial_sdf, partial_csv, partial_summary):
            try:
                if partial.exists():
                    partial.unlink()
            except OSError:
                pass
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1

    log("Completed chunk {}.".format(args.chunk))
    log(
        "Records: {:,} | parsed: {:,} | retained: {:,} | "
        "exact removed: {:,} | similarity removed: {:,} | unreadable: {:,}".format(
            total_records,
            parsed_records,
            retained_count,
            exact_duplicates_removed,
            similarity_removed,
            unreadable_records,
        )
    )
    log("Fingerprint comparisons: {:,}".format(comparisons_performed))
    log("Filtered SDF: {}".format(output_sdf))
    log("CSV report: {}".format(output_csv))
    log("JSON summary: {}".format(output_summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
