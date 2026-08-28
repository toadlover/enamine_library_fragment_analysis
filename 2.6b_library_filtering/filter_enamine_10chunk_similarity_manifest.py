#!/usr/bin/env python3
"""
Filter a span of Enamine REAL 2.6B chunks by stereo-insensitive exact identity
and Morgan/Tanimoto similarity, producing a manifest of ligands to keep/remove.

The default span is 10 consecutive chunks (up to ~500,000 ligands when each
chunk contains ~50,000 ligands). Molecules are ranked before similarity
filtering so that the preferred representative is:

    1. fewest heavy atoms
    2. lowest molecular weight
    3. earliest original input order (deterministic tie-break)

This is NOT exhaustive connected-component clustering. It is a greedy
representative filter: molecules are considered smallest-first and a molecule
is removed when it is >= the similarity cutoff to any already-retained
representative. Thus every removed molecule points to a retained molecule that
is at least as preferred by the ranking above.

The manifest preserves source superchunk/chunk/subchunk/archive/member/record
provenance for every ligand.
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
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator


ARCHIVE_TEMPLATE = "split_new_named_{}.sdf.tar.gz"
EXPECTED_SUBCHUNKS = tuple(range(10))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a consecutive span of Enamine chunks for stereo-insensitive "
            "exact duplicates and Morgan/Tanimoto-similar ligands, writing a "
            "keep/remove manifest with full provenance."
        )
    )
    parser.add_argument(
        "start_chunk",
        type=int,
        help="First global chunk number in the span, e.g. 26900.",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory in which manifest and summary outputs are created.",
    )
    parser.add_argument(
        "--num-chunks",
        type=int,
        default=10,
        help="Number of consecutive chunks to compare together. Default: 10.",
    )
    parser.add_argument(
        "--library-root",
        type=Path,
        default=Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion"),
        help="Root of the Enamine REAL 2.6B library.",
    )
    parser.add_argument(
        "--similarity-cutoff",
        type=float,
        default=0.75,
        help="Remove a molecule when Tanimoto similarity to a retained representative is >= this value. Default: 0.75.",
    )
    parser.add_argument("--radius", type=int, default=2, help="Morgan radius. Default: 2 (ECFP4-like).")
    parser.add_argument("--fp-size", type=int, default=2048, help="Morgan fingerprint bit length. Default: 2048.")
    parser.add_argument(
        "--max-ligands",
        type=int,
        default=500000,
        help="Safety ceiling for successfully parsed ligands across the entire span. 0 disables. Default: 500000.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=5000,
        help="Print progress every N input records. 0 disables. Default: 5000.",
    )
    parser.add_argument(
        "--similarity-progress-every",
        type=int,
        default=1000,
        help="Print similarity-filter progress every N ranked molecules. 0 disables. Default: 1000.",
    )
    parser.add_argument(
        "--allow-missing-subchunks",
        action="store_true",
        help="Continue when one or more expected subchunk archives are missing.",
    )
    parser.add_argument(
        "--allow-missing-chunks",
        action="store_true",
        help="Continue when a requested chunk directory is missing.",
    )
    parser.add_argument(
        "--keep-salts",
        action="store_true",
        help="Compare complete multicomponent records. By default only the largest covalent fragment is used.",
    )
    parser.add_argument("--name-property", default="_Name", help="SDF property used as ligand name. Default: _Name.")
    parser.add_argument("--quiet-rdkit", action="store_true", help="Suppress RDKit warning/error logs.")
    return parser.parse_args()


def log(message: str) -> None:
    print("[{}] {}".format(time.strftime("%Y-%m-%d %H:%M:%S"), message), flush=True)


def validate_args(args: argparse.Namespace) -> None:
    if args.start_chunk < 0:
        raise ValueError("start_chunk must be nonnegative.")
    if args.num_chunks < 1:
        raise ValueError("--num-chunks must be at least 1.")
    if not 0.0 < args.similarity_cutoff <= 1.0:
        raise ValueError("--similarity-cutoff must be in (0, 1].")
    if args.radius < 1:
        raise ValueError("--radius must be at least 1.")
    if args.fp_size < 128:
        raise ValueError("--fp-size must be at least 128.")
    if args.max_ligands < 0:
        raise ValueError("--max-ligands cannot be negative.")
    if args.progress_every < 0 or args.similarity_progress_every < 0:
        raise ValueError("progress intervals cannot be negative.")


def superchunk_for_chunk(chunk: int) -> int:
    return chunk // 100


def resolve_chunk_directory(library_root: Path, chunk: int) -> Path:
    root = library_root.resolve()
    if not root.is_dir():
        raise FileNotFoundError("Library root does not exist: {}".format(root))
    superchunk = superchunk_for_chunk(chunk)
    superchunk_dir = root / str(superchunk)
    chunk_dir = superchunk_dir / "{:05d}".format(chunk)
    if not chunk_dir.is_dir():
        raise FileNotFoundError("Expected chunk directory does not exist: {}".format(chunk_dir))
    return chunk_dir.resolve()


def collect_archives(chunk_dir: Path, allow_missing: bool) -> List[Tuple[int, Path]]:
    archives: List[Tuple[int, Path]] = []
    missing: List[Path] = []
    for subchunk in EXPECTED_SUBCHUNKS:
        archive = chunk_dir / ARCHIVE_TEMPLATE.format(subchunk)
        if archive.is_file():
            archives.append((subchunk, archive.resolve()))
        else:
            missing.append(archive)
    if missing and not allow_missing:
        raise FileNotFoundError(
            "Missing {} expected subchunk archive(s) in {}:\n  {}".format(
                len(missing), chunk_dir, "\n  ".join(str(p) for p in missing)
            )
        )
    if not archives:
        raise FileNotFoundError("No expected subchunk archives found in {}".format(chunk_dir))
    if missing:
        log("WARNING: chunk {} missing subchunks: {}".format(chunk_dir.name, ", ".join(p.name for p in missing)))
    return archives


def safe_sdf_members(archive: tarfile.TarFile) -> List[tarfile.TarInfo]:
    result = []
    for member in archive.getmembers():
        if not member.isfile():
            continue
        lower = Path(member.name).name.lower()
        if lower.endswith(".sdf") or lower.endswith(".sd"):
            result.append(member)
    return result


def decompress_archive_sdfs(archive_path: Path, destination: Path, chunk: int, subchunk: int) -> List[Path]:
    extracted: List[Path] = []
    with tarfile.open(str(archive_path), mode="r:gz") as archive:
        members = safe_sdf_members(archive)
        if not members:
            raise RuntimeError("Archive contains no SDF member: {}".format(archive_path))
        for member_number, member in enumerate(members):
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError("Could not read {} from {}".format(member.name, archive_path))
            out = destination / "chunk_{:05d}_subchunk_{}_member_{:03d}_{}".format(
                chunk, subchunk, member_number, Path(member.name).name
            )
            with source, out.open("wb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
            extracted.append(out)
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


def molecule_name(mol: Chem.Mol, global_index: int, name_property: str) -> str:
    if mol.HasProp(name_property):
        value = mol.GetProp(name_property).strip()
        if value:
            return value
    return "span_molecule_{:07d}".format(global_index + 1)


def report_fieldnames() -> List[str]:
    return [
        "global_input_index", "name", "status", "reason",
        "superchunk", "chunk", "subchunk", "archive_name", "archive_path",
        "sdf_member_name", "record_index_in_member",
        "heavy_atoms", "molecular_weight", "stereo_insensitive_smiles",
        "representative_global_index", "representative_name",
        "representative_superchunk", "representative_chunk", "representative_subchunk",
        "representative_archive_name", "representative_record_index_in_member",
        "representative_heavy_atoms", "representative_molecular_weight",
        "similarity_to_representative",
    ]


def make_base_record(global_index: int, name: str, chunk: int, subchunk: int,
                     archive_path: Path, sdf_member_path: Path, member_record_index: int) -> Dict[str, object]:
    return {
        "global_input_index": global_index,
        "name": name,
        "superchunk": superchunk_for_chunk(chunk),
        "chunk": chunk,
        "subchunk": subchunk,
        "archive_name": archive_path.name,
        "archive_path": str(archive_path),
        "sdf_member_name": sdf_member_path.name,
        "record_index_in_member": member_record_index,
        "heavy_atoms": "",
        "molecular_weight": "",
        "stereo_insensitive_smiles": "",
        "fingerprint": None,
        "status": "",
        "reason": "",
        "representative_index": None,
        "similarity": None,
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
        "representative_heavy_atoms": rep["heavy_atoms"],
        "representative_molecular_weight": "{:.6f}".format(float(rep["molecular_weight"])),
    }


def write_manifest(path: Path, records: List[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=report_fieldnames())
        writer.writeheader()
        for rec in sorted(records, key=lambda r: int(r["global_input_index"])):
            row = {k: rec.get(k, "") for k in report_fieldnames()}
            row["molecular_weight"] = "" if rec["molecular_weight"] == "" else "{:.6f}".format(float(rec["molecular_weight"]))
            row["similarity_to_representative"] = "" if rec["similarity"] is None else "{:.6f}".format(float(rec["similarity"]))
            for key in [
                "representative_global_index", "representative_name", "representative_superchunk",
                "representative_chunk", "representative_subchunk", "representative_archive_name",
                "representative_record_index_in_member", "representative_heavy_atoms",
                "representative_molecular_weight",
            ]:
                row[key] = ""
            rep_index = rec.get("representative_index")
            if rep_index is not None:
                row.update(representative_fields(records[int(rep_index)]))
            row.pop("fingerprint", None)
            row.pop("representative_index", None)
            row.pop("similarity", None)
            writer.writerow({k: row.get(k, "") for k in report_fieldnames()})


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
    chunks_requested = list(range(args.start_chunk, args.start_chunk + args.num_chunks))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = args.output_dir.resolve()
    span_label = "chunks_{:05d}_{:05d}".format(chunks_requested[0], chunks_requested[-1])
    output_csv = output_dir / (span_label + "_similarity_manifest.csv")
    output_summary = output_dir / (span_label + "_summary.json")
    partial_csv = output_dir / (output_csv.name + ".partial")
    partial_summary = output_dir / (output_summary.name + ".partial")
    for p in (partial_csv, partial_summary):
        if p.exists():
            p.unlink()

    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=args.radius, fpSize=args.fp_size, includeChirality=False
    )

    records: List[Dict[str, object]] = []
    parsed_records = unreadable_records = processing_failures = 0
    total_sdf_records = 0
    chunks_processed: List[int] = []
    archives_processed: List[str] = []

    temp_dir = Path(tempfile.mkdtemp(prefix=".{}_work_".format(span_label), dir=str(output_dir)))
    log("Chunk span: {}".format(", ".join(str(c) for c in chunks_requested)))
    log("Temporary extraction directory: {}".format(temp_dir))
    log("Similarity cutoff: {:.3f}; representative ranking: heavy atoms, then molecular weight".format(args.similarity_cutoff))

    try:
        for chunk_position, chunk in enumerate(chunks_requested, start=1):
            try:
                chunk_dir = resolve_chunk_directory(args.library_root, chunk)
                archives = collect_archives(chunk_dir, args.allow_missing_subchunks)
            except FileNotFoundError as exc:
                if args.allow_missing_chunks:
                    log("WARNING: skipping missing chunk {}: {}".format(chunk, exc))
                    continue
                raise

            chunks_processed.append(chunk)
            log("Chunk {} ({}/{}): found {} subchunk archive(s).".format(chunk, chunk_position, len(chunks_requested), len(archives)))

            for subchunk, archive_path in archives:
                archives_processed.append(str(archive_path))
                extracted_sdfs = decompress_archive_sdfs(archive_path, temp_dir, chunk, subchunk)
                for sdf_member_path in extracted_sdfs:
                    supplier = Chem.SDMolSupplier(
                        str(sdf_member_path), sanitize=True, removeHs=True, strictParsing=False
                    )
                    for member_record_index, mol in enumerate(supplier):
                        global_index = total_sdf_records
                        total_sdf_records += 1
                        rec = make_base_record(
                            global_index, "unreadable_{:07d}".format(global_index + 1),
                            chunk, subchunk, archive_path, sdf_member_path, member_record_index
                        )
                        records.append(rec)

                        if mol is None:
                            unreadable_records += 1
                            rec["status"] = "removed"
                            rec["reason"] = "SDF parsing or sanitization failed"
                            continue

                        parsed_records += 1
                        if args.max_ligands > 0 and parsed_records > args.max_ligands:
                            raise RuntimeError(
                                "Parsed ligand count exceeded --max-ligands {:,}; increase the ceiling or use 0 to disable.".format(args.max_ligands)
                            )
                        rec["name"] = molecule_name(mol, global_index, args.name_property)

                        try:
                            working = Chem.Mol(mol)
                            if not args.keep_salts:
                                working = largest_fragment(working)
                            rec["heavy_atoms"] = int(working.GetNumHeavyAtoms())
                            rec["molecular_weight"] = float(Descriptors.MolWt(working))
                            rec["stereo_insensitive_smiles"] = stereo_insensitive_key(working)
                            rec["fingerprint"] = fingerprint_generator.GetFingerprint(working)
                        except Exception as exc:
                            processing_failures += 1
                            rec["status"] = "removed"
                            rec["reason"] = "molecule processing failed: {}".format(exc)

                        if args.progress_every > 0 and total_sdf_records % args.progress_every == 0:
                            elapsed = max(time.time() - start_time, 1e-9)
                            log("Loaded {:,} records | parsed {:,} | {:.1f} records/s".format(total_sdf_records, parsed_records, total_sdf_records / elapsed))

                    try:
                        sdf_member_path.unlink()
                    except OSError:
                        pass

        valid_indices = [
            i for i, rec in enumerate(records)
            if rec["fingerprint"] is not None and rec["status"] == ""
        ]
        valid_indices.sort(key=lambda i: (
            int(records[i]["heavy_atoms"]),
            float(records[i]["molecular_weight"]),
            int(records[i]["global_input_index"]),
        ))

        log("Loaded {:,} usable ligands. Beginning smallest-first greedy similarity filtering.".format(len(valid_indices)))

        retained_record_indices: List[int] = []
        retained_fingerprints: List[object] = []
        exact_key_to_rep_record: Dict[str, int] = {}
        exact_duplicates_removed = similarity_removed = comparisons_performed = 0

        for rank_position, record_index in enumerate(valid_indices, start=1):
            rec = records[record_index]
            key = str(rec["stereo_insensitive_smiles"])

            exact_rep_record = exact_key_to_rep_record.get(key)
            if exact_rep_record is not None:
                rec["status"] = "remove"
                rec["reason"] = "exact duplicate after removing stereochemistry"
                rec["representative_index"] = exact_rep_record
                rec["similarity"] = 1.0
                exact_duplicates_removed += 1
            else:
                best_similarity = -1.0
                representative_position: Optional[int] = None
                if retained_fingerprints:
                    similarities = DataStructs.BulkTanimotoSimilarity(
                        rec["fingerprint"], retained_fingerprints
                    )
                    comparisons_performed += len(retained_fingerprints)
                    best_similarity = max(similarities)
                    if best_similarity >= args.similarity_cutoff:
                        representative_position = similarities.index(best_similarity)

                if representative_position is None:
                    rec["status"] = "keep"
                    rec["reason"] = "smallest-first greedy similarity representative"
                    rec["representative_index"] = record_index
                    rec["similarity"] = 1.0
                    retained_record_indices.append(record_index)
                    retained_fingerprints.append(rec["fingerprint"])
                    exact_key_to_rep_record[key] = record_index
                else:
                    rep_record_index = retained_record_indices[representative_position]
                    rec["status"] = "remove"
                    rec["reason"] = "Morgan/Tanimoto similarity >= {:.3f}".format(args.similarity_cutoff)
                    rec["representative_index"] = rep_record_index
                    rec["similarity"] = best_similarity
                    exact_key_to_rep_record[key] = rep_record_index
                    similarity_removed += 1

            if args.similarity_progress_every > 0 and rank_position % args.similarity_progress_every == 0:
                log(
                    "Similarity pass {:,}/{:,} | keep {:,} | exact remove {:,} | similar remove {:,} | comparisons {:,}".format(
                        rank_position, len(valid_indices), len(retained_record_indices), exact_duplicates_removed,
                        similarity_removed, comparisons_performed
                    )
                )

        write_manifest(partial_csv, records)

        elapsed_seconds = time.time() - start_time
        summary = {
            "library_root": str(args.library_root.resolve()),
            "chunks_requested": chunks_requested,
            "chunks_processed": chunks_processed,
            "archives_processed": archives_processed,
            "similarity_method": "smallest_first_greedy_morgan_tanimoto_stereo_insensitive",
            "clustering_semantics": "not exhaustive connected-component clustering; each removed ligand maps to a retained representative at or above cutoff",
            "representative_priority": ["fewest_heavy_atoms", "lowest_molecular_weight", "earliest_input_order"],
            "similarity_cutoff": args.similarity_cutoff,
            "morgan_radius": args.radius,
            "fingerprint_bits": args.fp_size,
            "keep_salts": args.keep_salts,
            "total_sdf_records": total_sdf_records,
            "successfully_parsed_ligands": parsed_records,
            "unreadable_records": unreadable_records,
            "processing_failures": processing_failures,
            "exact_stereo_insensitive_duplicates_removed": exact_duplicates_removed,
            "similarity_removed": similarity_removed,
            "retained_ligands": len(retained_record_indices),
            "fingerprint_comparisons_performed": comparisons_performed,
            "theoretical_all_pairs_for_usable_ligands": len(valid_indices) * (len(valid_indices) - 1) // 2,
            "elapsed_seconds": elapsed_seconds,
            "output_manifest": str(output_csv),
        }
        with partial_summary.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        os.replace(str(partial_csv), str(output_csv))
        os.replace(str(partial_summary), str(output_summary))

    except Exception as exc:
        for p in (partial_csv, partial_summary):
            try:
                if p.exists():
                    p.unlink()
            except OSError:
                pass
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
    finally:
        try:
            shutil.rmtree(str(temp_dir), ignore_errors=True)
        except Exception:
            pass

    log("Completed span {}-{}.".format(chunks_requested[0], chunks_requested[-1]))
    log("Manifest: {}".format(output_csv))
    log("Summary: {}".format(output_summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
