#!/usr/bin/env python3
"""
filter_similar_sdf.py

Remove structurally redundant molecules from an SDF file.

Workflow
--------
1. Read and sanitize molecules.
2. Optionally retain only the largest covalent fragment.
3. Collapse exact duplicates after removing stereochemistry. This guarantees
   that enantiomers/diastereomers differing only in assigned stereochemistry
   are treated as duplicates.
4. Generate stereo-insensitive Morgan fingerprints.
5. Cluster remaining molecules using the Butina algorithm and a configurable
   Tanimoto similarity cutoff.
6. Keep one representative (the Butina cluster centroid) from each cluster.
7. Write retained molecules to SDF and all decisions to CSV.

Python 3.8+ compatible.

Example
-------
python filter_similar_sdf.py \
    input.sdf \
    filtered.sdf \
    --report similarity_report.csv \
    --similarity-cutoff 0.75
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remove enantiomers, stereo variants, exact duplicates, and "
            "fingerprint-similar molecules from an SDF."
        )
    )
    parser.add_argument("input_sdf", type=Path, help="Input SDF file.")
    parser.add_argument("output_sdf", type=Path, help="Output SDF of retained molecules.")
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("similarity_report.csv"),
        help="CSV report describing retained and removed molecules.",
    )
    parser.add_argument(
        "--similarity-cutoff",
        type=float,
        default=0.75,
        help=(
            "Tanimoto similarity at or above which molecules are clustered "
            "as similar. Default: 0.75."
        ),
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
        "--keep-salts",
        action="store_true",
        help=(
            "Do not reduce multicomponent records to their largest fragment. "
            "By default, the largest covalent fragment is used for comparisons "
            "and written to the output."
        ),
    )
    parser.add_argument(
        "--name-property",
        default="_Name",
        help="SDF property used as the ligand name. Default: _Name.",
    )
    return parser.parse_args()


def molecule_name(mol: Chem.Mol, index: int, property_name: str) -> str:
    if mol.HasProp(property_name):
        value = mol.GetProp(property_name).strip()
        if value:
            return value
    return "molecule_{:06d}".format(index + 1)


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    """Return the fragment with the greatest heavy-atom count."""
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        return Chem.Mol(mol)
    return max(
        fragments,
        key=lambda fragment: (
            fragment.GetNumHeavyAtoms(),
            fragment.GetNumAtoms(),
            Chem.MolToSmiles(fragment, isomericSmiles=False),
        ),
    )


def stereo_insensitive_key(mol: Chem.Mol) -> str:
    """
    Generate an exact connectivity key with stereochemistry removed.

    Isotopes, formal charges, bond orders, and atom identities remain relevant.
    """
    copy = Chem.Mol(mol)
    Chem.RemoveStereochemistry(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def load_molecules(
    input_path: Path,
    keep_salts: bool,
    name_property: str,
) -> Tuple[List[Chem.Mol], List[Dict[str, object]]]:
    supplier = Chem.SDMolSupplier(
        str(input_path),
        sanitize=True,
        removeHs=True,
        strictParsing=False,
    )

    molecules: List[Chem.Mol] = []
    records: List[Dict[str, object]] = []

    for sdf_index, mol in enumerate(supplier):
        if mol is None:
            records.append(
                {
                    "input_index": sdf_index,
                    "name": "unreadable_record_{:06d}".format(sdf_index + 1),
                    "status": "removed",
                    "reason": "SDF parsing or sanitization failed",
                    "representative_input_index": "",
                    "representative_name": "",
                    "similarity_to_representative": "",
                    "stereo_insensitive_smiles": "",
                    "cluster_id": "",
                }
            )
            continue

        working = Chem.Mol(mol)
        if not keep_salts:
            try:
                working = largest_fragment(working)
            except Exception as exc:
                records.append(
                    {
                        "input_index": sdf_index,
                        "name": molecule_name(mol, sdf_index, name_property),
                        "status": "removed",
                        "reason": "fragment processing failed: {}".format(exc),
                        "representative_input_index": "",
                        "representative_name": "",
                        "similarity_to_representative": "",
                        "stereo_insensitive_smiles": "",
                        "cluster_id": "",
                    }
                )
                continue

        name = molecule_name(mol, sdf_index, name_property)
        working.SetProp("_Name", name)
        working.SetIntProp("_OriginalSDFIndex", sdf_index)

        molecules.append(working)
        records.append(
            {
                "input_index": sdf_index,
                "name": name,
                "status": "pending",
                "reason": "",
                "representative_input_index": "",
                "representative_name": "",
                "similarity_to_representative": "",
                "stereo_insensitive_smiles": stereo_insensitive_key(working),
                "cluster_id": "",
            }
        )

    return molecules, records


def collapse_exact_stereo_insensitive_duplicates(
    molecules: Sequence[Chem.Mol],
    records_by_input_index: Dict[int, Dict[str, object]],
) -> List[Chem.Mol]:
    unique: List[Chem.Mol] = []
    first_by_key: Dict[str, Chem.Mol] = {}

    for mol in molecules:
        input_index = mol.GetIntProp("_OriginalSDFIndex")
        key = stereo_insensitive_key(mol)

        if key not in first_by_key:
            first_by_key[key] = mol
            unique.append(mol)
            continue

        representative = first_by_key[key]
        representative_index = representative.GetIntProp("_OriginalSDFIndex")
        representative_name = representative.GetProp("_Name")

        record = records_by_input_index[input_index]
        record["status"] = "removed"
        record["reason"] = "exact duplicate after removing stereochemistry"
        record["representative_input_index"] = representative_index
        record["representative_name"] = representative_name
        record["similarity_to_representative"] = "1.000000"

    return unique


def butina_clusters(
    fingerprints: Sequence[object],
    similarity_cutoff: float,
) -> Tuple[Tuple[int, ...], ...]:
    """
    Cluster fingerprints using a condensed lower-triangle distance matrix.

    Butina expects a distance threshold, so:
        distance_threshold = 1 - Tanimoto similarity cutoff
    """
    distances: List[float] = []
    for i in range(1, len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[i], fingerprints[:i]
        )
        distances.extend(1.0 - similarity for similarity in similarities)

    return Butina.ClusterData(
        distances,
        len(fingerprints),
        1.0 - similarity_cutoff,
        isDistData=True,
        reordering=True,
    )


def filter_by_similarity(
    molecules: Sequence[Chem.Mol],
    records_by_input_index: Dict[int, Dict[str, object]],
    similarity_cutoff: float,
    radius: int,
    fp_size: int,
) -> List[Chem.Mol]:
    if not molecules:
        return []

    # Chirality is intentionally excluded. Enantiomers therefore generate
    # identical fingerprints, while close analogs can exceed the cutoff.
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        includeChirality=False,
    )
    fingerprints = [generator.GetFingerprint(mol) for mol in molecules]
    clusters = butina_clusters(fingerprints, similarity_cutoff)

    retained: List[Chem.Mol] = []

    for cluster_number, cluster in enumerate(clusters, start=1):
        representative_position = cluster[0]
        representative = molecules[representative_position]
        representative_fp = fingerprints[representative_position]
        representative_index = representative.GetIntProp("_OriginalSDFIndex")
        representative_name = representative.GetProp("_Name")

        representative.SetIntProp("_SimilarityCluster", cluster_number)
        representative.SetIntProp("_SimilarityClusterSize", len(cluster))
        representative.SetDoubleProp("_SimilarityCutoff", similarity_cutoff)
        retained.append(representative)

        representative_record = records_by_input_index[representative_index]
        representative_record["status"] = "retained"
        representative_record["reason"] = "cluster representative"
        representative_record["representative_input_index"] = representative_index
        representative_record["representative_name"] = representative_name
        representative_record["similarity_to_representative"] = "1.000000"
        representative_record["cluster_id"] = cluster_number

        for member_position in cluster[1:]:
            member = molecules[member_position]
            member_index = member.GetIntProp("_OriginalSDFIndex")
            similarity = DataStructs.TanimotoSimilarity(
                fingerprints[member_position], representative_fp
            )

            record = records_by_input_index[member_index]
            record["status"] = "removed"
            record["reason"] = (
                "Morgan/Tanimoto similarity >= {:.3f}".format(similarity_cutoff)
            )
            record["representative_input_index"] = representative_index
            record["representative_name"] = representative_name
            record["similarity_to_representative"] = "{:.6f}".format(similarity)
            record["cluster_id"] = cluster_number

    # Preserve input order in the filtered SDF.
    retained.sort(key=lambda mol: mol.GetIntProp("_OriginalSDFIndex"))
    return retained


def write_sdf(path: Path, molecules: Sequence[Chem.Mol]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    if writer is None:
        raise RuntimeError("Could not open output SDF: {}".format(path))

    try:
        for mol in molecules:
            writer.write(mol)
    finally:
        writer.close()


def write_report(path: Path, records: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "input_index",
        "name",
        "status",
        "reason",
        "representative_input_index",
        "representative_name",
        "similarity_to_representative",
        "stereo_insensitive_smiles",
        "cluster_id",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in sorted(records, key=lambda item: int(item["input_index"])):
            writer.writerow(record)


def main() -> int:
    args = parse_args()

    if not args.input_sdf.is_file():
        print("ERROR: input file does not exist: {}".format(args.input_sdf), file=sys.stderr)
        return 2
    if not 0.0 < args.similarity_cutoff <= 1.0:
        print("ERROR: --similarity-cutoff must be in (0, 1].", file=sys.stderr)
        return 2
    if args.radius < 1:
        print("ERROR: --radius must be at least 1.", file=sys.stderr)
        return 2
    if args.fp_size < 128:
        print("ERROR: --fp-size must be at least 128.", file=sys.stderr)
        return 2

    molecules, records = load_molecules(
        args.input_sdf,
        keep_salts=args.keep_salts,
        name_property=args.name_property,
    )
    records_by_input_index = {
        int(record["input_index"]): record for record in records
    }

    unique_molecules = collapse_exact_stereo_insensitive_duplicates(
        molecules, records_by_input_index
    )
    retained = filter_by_similarity(
        unique_molecules,
        records_by_input_index,
        similarity_cutoff=args.similarity_cutoff,
        radius=args.radius,
        fp_size=args.fp_size,
    )

    # Any still-pending record indicates an internal logic error.
    for record in records:
        if record["status"] == "pending":
            record["status"] = "removed"
            record["reason"] = "internal filtering error"

    write_sdf(args.output_sdf, retained)
    write_report(args.report, records)

    parsed_count = len(molecules)
    failed_count = sum(
        1 for record in records
        if str(record["reason"]).startswith("SDF parsing")
    )
    exact_duplicate_count = sum(
        1 for record in records
        if record["reason"] == "exact duplicate after removing stereochemistry"
    )
    similarity_removed_count = sum(
        1 for record in records
        if str(record["reason"]).startswith("Morgan/Tanimoto")
    )

    print("Input records:              {}".format(len(records)))
    print("Successfully parsed:        {}".format(parsed_count))
    print("Unreadable records:         {}".format(failed_count))
    print("Stereo/exact duplicates:    {}".format(exact_duplicate_count))
    print("Similarity-based removals:  {}".format(similarity_removed_count))
    print("Retained molecules:         {}".format(len(retained)))
    print("Filtered SDF:               {}".format(args.output_sdf))
    print("Decision report:            {}".format(args.report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
