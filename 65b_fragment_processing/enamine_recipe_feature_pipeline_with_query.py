#!/usr/bin/env python3
"""
Large-scale Enamine recipe + feature + compact LSH extraction pipeline,
with an added query mode for searching compact batch binary files by SMILES.

Modes
-----
1) submit
   Recursively discover .bz2 files, create a manifest, and print/submit an LSF
   array job.

2) worker
   Process one manifest task and write compact production outputs by default:

      batch_XXXXXX.ligand_map.tsv.bz2
      batch_XXXXXX.feature_hashes.bin
      batch_XXXXXX.lsh_packed.bin
      batch_XXXXXX.binary_schema.json
      batch_XXXXXX.stats.json

3) query
   Given a SMILES string and one batch directory/prefix, decompose the query with
   the same fragment/linker recipe logic, hash its features, and scan the compact
   batch binaries to answer:

      * Which query fragment/linker/motif features exist in this batch?
      * How many batch ligands contain each query feature?
      * Which batch ligands share the most query features?
      * Which batch ligands share LSH buckets with the query?

The query mode does not need readable feature TSVs. It works directly from
feature_hashes.bin and lsh_packed.bin. It can join top hits back to the
ligand_map.tsv.bz2 for source-file/row references.

Important notes
---------------
* local_ligand_id is a batch-local metadata key, not a model feature.
* Hashes/LSH buckets are comparable across runs only if feature schema,
  RDKit behavior, linker cutoff, hash seed, and LSH settings are the same.
* Binary feature/LSH files are compact storage. Convert them to human-readable
  form or tensors downstream as needed.
"""

import argparse
import bz2
import csv
import hashlib
import json
import os
import shlex
import struct
import subprocess
import sys
from collections import defaultdict, deque
from contextlib import ExitStack
from heapq import heappush, heappushpop, nlargest
from itertools import combinations
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdBase

RDLogger.DisableLog("rdApp.warning")

ROT_BOND_SMARTS = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")
DEFAULT_LINKER_MAX_HEAVY_ATOMS = 4
DEFAULT_BATCH_SIZE = 3_000_000
MAX_EXAMPLES = 3

# Production defaults. 64 x 8 gives 8 LSH bands per ligand.
DEFAULT_MINHASH_PERM = 64
DEFAULT_LSH_ROWS_PER_BAND = 8
FEATURE_SCHEMA_VERSION = "recipe_features_v1"
SCRIPT_VERSION = "compact_lsh_query_v1"

FEATURE_RECORD_HEADER = struct.Struct("<IH")  # uint32 ligand_id, uint16 n_features
LSH_MEMBER_RECORD = struct.Struct("<HQI")     # uint16 band, uint64 bucket_hash, uint32 ligand_id


# ----------------------------- small utilities -----------------------------


def stable_hash_u64(text, seed=0):
    """Stable 64-bit unsigned hash, little-endian, used for feature and LSH hashes."""
    data = (str(seed) + "\x1f" + text).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(data, digest_size=8).digest(), "little")


def qualify_ligand(source_file, ligand_name):
    return f"{source_file}::{ligand_name}"


def sorted_count_items(counter):
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))


def count_keys(counter):
    return ";".join(k for k, v in sorted_count_items(counter))


def count_values(counter):
    return ";".join(str(v) for k, v in sorted_count_items(counter))


def add_count_entry(summary_dict, key, ligand_name):
    if key not in summary_dict:
        summary_dict[key] = {
            "count": 1,
            "examples": [ligand_name],
            "connected_fragments": defaultdict(int),
            "connected_linkers": defaultdict(int),
            "connected_terminal_fragments": defaultdict(int),
            "connection_patterns": defaultdict(int),
        }
    else:
        summary_dict[key]["count"] += 1
        if len(summary_dict[key]["examples"]) < MAX_EXAMPLES:
            summary_dict[key]["examples"].append(ligand_name)


def add_neighbors(summary_dict, key, fragments=None, linkers=None, terminals=None):
    if key not in summary_dict:
        return
    for value in fragments or []:
        summary_dict[key]["connected_fragments"][value] += 1
    for value in linkers or []:
        summary_dict[key]["connected_linkers"][value] += 1
    for value in terminals or []:
        summary_dict[key]["connected_terminal_fragments"][value] += 1


def add_pattern(summary_dict, key, pattern):
    if key in summary_dict:
        summary_dict[key]["connection_patterns"][pattern] += 1


def parse_smiles_line(line, line_number):
    line = line.strip()
    if not line:
        return None, None
    fields = line.split("\t")
    if len(fields) < 2:
        fields = line.split()
    if len(fields) < 2:
        return None, None
    smiles = fields[0]
    ligand_name = fields[1]
    if smiles.lower() in {"smiles", "cxsmiles"}:
        return None, None
    if not ligand_name:
        ligand_name = f"line_{line_number}"
    return smiles, ligand_name


# ------------------------- molecule recipe extraction ------------------------


def heavy_atom_count_for_atom_set(mol, atom_indices):
    return sum(1 for idx in atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)


def mol_fragment_smiles(mol, atom_indices):
    return Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=sorted(atom_indices),
        canonical=True,
        isomericSmiles=False,
    )


def get_cut_components(mol):
    rot_bonds = list(mol.GetSubstructMatches(ROT_BOND_SMARTS))
    em = Chem.EditableMol(mol)
    for atom1, atom2 in rot_bonds:
        em.RemoveBond(atom1, atom2)
    cut_mol = em.GetMol()

    component_atom_tuples = Chem.GetMolFrags(cut_mol, asMols=False, sanitizeFrags=False)
    component_atoms = [set(x) for x in component_atom_tuples]

    atom_to_component = {}
    for component_id, atom_set in enumerate(component_atoms):
        for atom_idx in atom_set:
            atom_to_component[atom_idx] = component_id

    component_graph = defaultdict(set)
    cut_edges = []
    for atom1, atom2 in rot_bonds:
        c1 = atom_to_component[atom1]
        c2 = atom_to_component[atom2]
        if c1 == c2:
            continue
        component_graph[c1].add(c2)
        component_graph[c2].add(c1)
        cut_edges.append({"atom1": atom1, "atom2": atom2, "component1": c1, "component2": c2})

    for component_id in range(len(component_atoms)):
        component_graph[component_id] = set(component_graph[component_id])

    return rot_bonds, component_atoms, component_graph, cut_edges


def classify_and_recipe_molecule(mol, linker_max_heavy_atoms):
    rot_bonds, component_atoms, component_graph, cut_edges = get_cut_components(mol)

    component_info = {}
    for component_id, atom_set in enumerate(component_atoms):
        heavy_atoms = heavy_atom_count_for_atom_set(mol, atom_set)
        degree = len(component_graph[component_id])
        smiles = mol_fragment_smiles(mol, atom_set)
        is_small = heavy_atoms <= linker_max_heavy_atoms
        component_info[component_id] = {
            "component_id": component_id,
            "atoms": sorted(atom_set),
            "heavy_atoms": heavy_atoms,
            "degree": degree,
            "smiles": smiles,
            "is_small": is_small,
            "role": None,
        }

    candidate_linker_components = {
        cid for cid, info in component_info.items()
        if info["is_small"] and info["degree"] >= 2
    }

    linker_groups = []
    visited = set()
    for start in sorted(candidate_linker_components):
        if start in visited:
            continue
        group = set()
        queue = deque([start])
        while queue:
            current = queue.popleft()
            if current in visited or current not in candidate_linker_components:
                continue
            visited.add(current)
            group.add(current)
            for nbr in component_graph[current]:
                if nbr in candidate_linker_components and nbr not in visited:
                    queue.append(nbr)
        if group:
            linker_groups.append(group)

    true_linker_groups = []
    component_to_linker_group = {}
    for group in linker_groups:
        outside_neighbors = set()
        for cid in group:
            for nbr in component_graph[cid]:
                if nbr not in group:
                    outside_neighbors.add(nbr)
        if len(outside_neighbors) >= 2:
            lid = len(true_linker_groups)
            true_linker_groups.append(group)
            for cid in group:
                component_to_linker_group[cid] = lid

    for cid, info in component_info.items():
        if cid in component_to_linker_group:
            info["role"] = "linker_component"
        elif len(component_atoms) == 1:
            info["role"] = "singleton_molecule"
        elif info["is_small"] and info["degree"] <= 1:
            info["role"] = "terminal_group"
        else:
            info["role"] = "fragment"

    linker_objects = []
    for lid, group in enumerate(true_linker_groups):
        linker_atoms = set()
        connected_components = set()
        for cid in group:
            linker_atoms.update(component_info[cid]["atoms"])
            for nbr in component_graph[cid]:
                if nbr not in group:
                    connected_components.add(nbr)
        linker_smiles = mol_fragment_smiles(mol, linker_atoms)
        linker_objects.append({
            "linker_id": f"L{lid}",
            "component_ids": sorted(group),
            "smiles": linker_smiles,
            "heavy_atoms": heavy_atom_count_for_atom_set(mol, linker_atoms),
            "connected_component_ids": sorted(connected_components),
        })

    non_linker_components = []
    component_id_to_recipe_id = {}
    for cid, info in sorted(component_info.items()):
        if info["role"] == "linker_component":
            continue
        fid = f"F{len(non_linker_components)}"
        component_id_to_recipe_id[cid] = fid
        non_linker_components.append({
            "fragment_id": fid,
            "component_id": cid,
            "smiles": info["smiles"],
            "role": info["role"],
            "heavy_atoms": info["heavy_atoms"],
            "degree": info["degree"],
        })

    for linker in linker_objects:
        linker["connected_fragment_ids"] = [
            component_id_to_recipe_id[c]
            for c in linker["connected_component_ids"]
            if c in component_id_to_recipe_id
        ]
        linker["connected_fragment_smiles"] = [
            component_info[c]["smiles"]
            for c in linker["connected_component_ids"]
            if c in component_id_to_recipe_id
        ]

    edges = []
    seen_edges = set()
    for edge in cut_edges:
        c1 = edge["component1"]
        c2 = edge["component2"]
        l1 = component_to_linker_group.get(c1)
        l2 = component_to_linker_group.get(c2)

        if l1 is None and l2 is None:
            node1 = component_id_to_recipe_id.get(c1)
            node2 = component_id_to_recipe_id.get(c2)
            etype = "direct"
        elif l1 is not None and l2 is None:
            node1 = f"L{l1}"
            node2 = component_id_to_recipe_id.get(c2)
            etype = "linker_attachment"
        elif l1 is None and l2 is not None:
            node1 = component_id_to_recipe_id.get(c1)
            node2 = f"L{l2}"
            etype = "linker_attachment"
        else:
            continue

        if node1 is None or node2 is None:
            continue
        ekey = tuple(sorted([node1, node2]) + [etype])
        if ekey in seen_edges:
            continue
        seen_edges.add(ekey)
        edges.append({"source": node1, "target": node2, "type": etype, "cut_atoms": [edge["atom1"], edge["atom2"]]})

    recipe = {
        "fragments": non_linker_components,
        "linkers": linker_objects,
        "edges": edges,
        "num_rotatable_bonds_cut": len(rot_bonds),
        "linker_max_heavy_atoms": linker_max_heavy_atoms,
    }
    return component_info, linker_objects, recipe


# -------------------------- features and summaries --------------------------


def connection_pattern(fragment_a, linker_smiles, fragment_b):
    left, right = sorted([fragment_a, fragment_b])
    return f"{left}:{linker_smiles}:{right}"


def recipe_maps(recipe):
    fragments_by_id = {f["fragment_id"]: f for f in recipe["fragments"]}
    linkers_by_id = {l["linker_id"]: l for l in recipe["linkers"]}
    return fragments_by_id, linkers_by_id


def make_feature_tokens(recipe):
    """Return recipe-derived feature tokens used for feature hashes and LSH."""
    fragments_by_id, linkers_by_id = recipe_maps(recipe)
    tokens = set()

    for frag in recipe["fragments"]:
        role = frag["role"]
        smi = frag["smiles"]
        tokens.add(f"COMP:{role}:{smi}")
        if role == "fragment":
            tokens.add(f"FRAG:{smi}")
        elif role == "terminal_group":
            tokens.add(f"TERM:{smi}")
        elif role == "singleton_molecule":
            tokens.add(f"SINGLETON:{smi}")

    for linker in recipe["linkers"]:
        tokens.add(f"LINKER:{linker['smiles']}")

    for edge in recipe["edges"]:
        s = edge["source"]
        t = edge["target"]
        if edge["type"] == "direct" and s in fragments_by_id and t in fragments_by_id:
            a = fragments_by_id[s]["smiles"]
            b = fragments_by_id[t]["smiles"]
            left, right = sorted([a, b])
            tokens.add(f"DIRECT:{left}:[direct]:{right}")
            tokens.add(f"PAIR:{left}|{right}")

    for linker in recipe["linkers"]:
        connected_ids = [fid for fid in linker.get("connected_fragment_ids", []) if fid in fragments_by_id]
        for fa_id, fb_id in combinations(sorted(connected_ids), 2):
            a = fragments_by_id[fa_id]["smiles"]
            b = fragments_by_id[fb_id]["smiles"]
            left, right = sorted([a, b])
            tokens.add(f"MOTIF:{left}:{linker['smiles']}:{right}")
            tokens.add(f"PAIR:{left}|{right}")
            tokens.add(f"FRAG_LINKER:{left}:{linker['smiles']}")
            tokens.add(f"FRAG_LINKER:{right}:{linker['smiles']}")

    topo = "TOPO:nfrag={}:nlink={}:ndirect={}".format(
        len(recipe["fragments"]),
        len(recipe["linkers"]),
        sum(1 for e in recipe["edges"] if e["type"] == "direct"),
    )
    tokens.add(topo)
    return sorted(tokens)


def feature_type_from_token(token):
    return token.split(":", 1)[0]


def feature_hashes(tokens, seed=0):
    return sorted({stable_hash_u64(tok, seed=seed) for tok in tokens})


def minhash_signature(tokens, num_perm):
    if not tokens:
        return [0] * num_perm
    sig = []
    for seed in range(num_perm):
        sig.append(min(stable_hash_u64(tok, seed=seed) for tok in tokens))
    return sig


def lsh_bucket_hashes(tokens, num_perm=DEFAULT_MINHASH_PERM, rows_per_band=DEFAULT_LSH_ROWS_PER_BAND):
    sig = minhash_signature(tokens, num_perm)
    if num_perm % rows_per_band != 0:
        raise ValueError("num_perm must be divisible by rows_per_band")
    hashes = []
    band_count = num_perm // rows_per_band
    for band in range(band_count):
        start = band * rows_per_band
        band_values = sig[start:start + rows_per_band]
        band_text = ",".join(str(x) for x in band_values)
        hashes.append(stable_hash_u64(band_text, seed=band))
    return hashes


def update_summaries(qualified_ligand, recipe, fragment_summary, linker_summary, connection_summary):
    fragments_by_id, linkers_by_id = recipe_maps(recipe)

    fragment_key_by_id = {}
    for frag in recipe["fragments"]:
        key = (frag["smiles"], frag["role"], frag["heavy_atoms"])
        fragment_key_by_id[frag["fragment_id"]] = key
        add_count_entry(fragment_summary, key, qualified_ligand)

    linker_key_by_id = {}
    for linker in recipe["linkers"]:
        key = (linker["smiles"], linker["heavy_atoms"])
        linker_key_by_id[linker["linker_id"]] = key
        add_count_entry(linker_summary, key, qualified_ligand)

    for edge in recipe["edges"]:
        s = edge["source"]
        t = edge["target"]
        if edge["type"] == "direct" and s in fragments_by_id and t in fragments_by_id:
            sf = fragments_by_id[s]
            tf = fragments_by_id[t]
            skey = fragment_key_by_id[s]
            tkey = fragment_key_by_id[t]

            add_neighbors(fragment_summary, skey, terminals=[tf["smiles"]] if tf["role"] == "terminal_group" else None,
                          fragments=[tf["smiles"]] if tf["role"] != "terminal_group" else None)
            add_neighbors(fragment_summary, tkey, terminals=[sf["smiles"]] if sf["role"] == "terminal_group" else None,
                          fragments=[sf["smiles"]] if sf["role"] != "terminal_group" else None)

            pat = connection_pattern(sf["smiles"], "[direct]", tf["smiles"])
            add_pattern(fragment_summary, skey, pat)
            add_pattern(fragment_summary, tkey, pat)
            conn_key = ("direct", "[direct]", ";".join(sorted([sf["smiles"], tf["smiles"]])))
            add_count_entry(connection_summary, conn_key, qualified_ligand)

        elif edge["type"] == "linker_attachment":
            if s in linkers_by_id and t in fragments_by_id:
                lid, fid = s, t
            elif t in linkers_by_id and s in fragments_by_id:
                lid, fid = t, s
            else:
                continue
            linker = linkers_by_id[lid]
            frag = fragments_by_id[fid]
            add_neighbors(fragment_summary, fragment_key_by_id[fid], linkers=[linker["smiles"]])
            add_neighbors(linker_summary, linker_key_by_id[lid],
                          terminals=[frag["smiles"]] if frag["role"] == "terminal_group" else None,
                          fragments=[frag["smiles"]] if frag["role"] != "terminal_group" else None)

    for linker in recipe["linkers"]:
        lid = linker["linker_id"]
        lkey = linker_key_by_id[lid]
        connected_ids = [fid for fid in linker.get("connected_fragment_ids", []) if fid in fragments_by_id]
        connected_smiles = sorted([fragments_by_id[fid]["smiles"] for fid in connected_ids])
        add_count_entry(connection_summary, ("linker", linker["smiles"], ";".join(connected_smiles)), qualified_ligand)
        for fa_id, fb_id in combinations(sorted(connected_ids), 2):
            fa = fragments_by_id[fa_id]
            fb = fragments_by_id[fb_id]
            pat = connection_pattern(fa["smiles"], linker["smiles"], fb["smiles"])
            add_pattern(linker_summary, lkey, pat)
            add_pattern(fragment_summary, fragment_key_by_id[fa_id], pat)
            add_pattern(fragment_summary, fragment_key_by_id[fb_id], pat)
            add_neighbors(fragment_summary, fragment_key_by_id[fa_id],
                          terminals=[fb["smiles"]] if fb["role"] == "terminal_group" else None,
                          fragments=[fb["smiles"]] if fb["role"] != "terminal_group" else None)
            add_neighbors(fragment_summary, fragment_key_by_id[fb_id],
                          terminals=[fa["smiles"]] if fa["role"] == "terminal_group" else None,
                          fragments=[fa["smiles"]] if fa["role"] != "terminal_group" else None)


# --------------------------- compact binary writing --------------------------


def write_feature_hash_record(handle, local_ligand_id, hashes):
    if len(hashes) > 65535:
        raise ValueError(f"Too many feature hashes for ligand {local_ligand_id}: {len(hashes)}")
    handle.write(FEATURE_RECORD_HEADER.pack(local_ligand_id, len(hashes)))
    if hashes:
        handle.write(struct.pack("<" + "Q" * len(hashes), *hashes))


def write_lsh_packed_record(handle, local_ligand_id, bucket_hashes):
    handle.write(struct.pack("<I" + "Q" * len(bucket_hashes), local_ligand_id, *bucket_hashes))


def write_lsh_membership_records(handle, local_ligand_id, bucket_hashes):
    for band, bucket_hash in enumerate(bucket_hashes):
        handle.write(LSH_MEMBER_RECORD.pack(band, bucket_hash, local_ligand_id))


def write_binary_schema(prefix, args, band_count, paths):
    schema = {
        "script_version": SCRIPT_VERSION,
        "feature_schema_version": FEATURE_SCHEMA_VERSION,
        "rdkit_version": rdBase.rdkitVersion,
        "endianness": "little",
        "local_ligand_id": "uint32; unique only within this batch directory",
        "global_ligand_reference": "Use batch/task id plus local_ligand_id, or join through ligand_map.tsv.bz2.",
        "feature_hashes_bin": {
            "path": paths.get("feature_hashes_bin"),
            "record": [
                {"name": "local_ligand_id", "type": "uint32"},
                {"name": "num_feature_hashes", "type": "uint16"},
                {"name": "feature_hashes", "type": "uint64[num_feature_hashes]"},
            ],
            "feature_hash": "blake2b-64 over recipe feature token with feature_hash_seed",
            "feature_hash_seed": args.feature_hash_seed,
        },
        "lsh_packed_bin": {
            "path": paths.get("lsh_packed_bin"),
            "record": [
                {"name": "local_ligand_id", "type": "uint32"},
                {"name": "bucket_hashes", "type": f"uint64[{band_count}]"},
            ],
            "minhash_perm": args.minhash_perm,
            "lsh_rows_per_band": args.lsh_rows_per_band,
            "band_count": band_count,
            "bucket_hash": "blake2b-64 of comma-joined minhash values in each band, seeded by band index",
        },
        "lsh_memberships_bin": {
            "path": paths.get("lsh_memberships_bin"),
            "written": bool(paths.get("lsh_memberships_bin")),
            "record": [
                {"name": "band", "type": "uint16"},
                {"name": "bucket_hash", "type": "uint64"},
                {"name": "local_ligand_id", "type": "uint32"},
            ],
            "note": "Optional unfolded representation for direct bucket grouping. Packed LSH can be unfolded later, so this is off by default.",
        },
        "ligand_map": {
            "path": paths.get("ligand_map"),
            "columns": paths.get("ligand_map_columns"),
            "source_data_row": "0-based index among parsed data rows in the source bz2 file; header/non-data rows are ignored.",
        },
    }
    path = str(prefix) + ".binary_schema.json"
    with open(path, "w") as handle:
        json.dump(schema, handle, indent=2, sort_keys=True)
    return path


# ----------------------------- file processing ------------------------------


def iter_ligands_in_range(input_path, start_row, end_row):
    """Yield data rows with 0-based data-row index in [start_row, end_row)."""
    data_index = 0
    with bz2.open(input_path, "rt", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            smiles, ligand_name = parse_smiles_line(line, line_number)
            if smiles is None:
                continue
            if data_index >= end_row:
                break
            if data_index >= start_row:
                yield data_index, smiles, ligand_name
            data_index += 1


def count_data_rows(input_path):
    count = 0
    with bz2.open(input_path, "rt", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            smiles, ligand_name = parse_smiles_line(line, line_number)
            if smiles is not None:
                count += 1
    return count


def ligand_map_header(args):
    columns = ["local_ligand_id", "source_file", "source_data_row"]
    if args.map_include_ligand_name:
        columns.append("ligand_name")
    if args.map_include_smiles:
        columns.append("smiles")
    return columns


def ligand_map_row(local_ligand_id, source_file, row_index, smiles, ligand_name, args):
    row = [str(local_ligand_id), source_file, str(row_index)]
    if args.map_include_ligand_name:
        row.append(ligand_name)
    if args.map_include_smiles:
        row.append(smiles)
    return row


def process_manifest_task(task, args):
    if args.minhash_perm % args.lsh_rows_per_band != 0:
        raise SystemExit("--minhash-perm must be divisible by --lsh-rows-per-band")
    band_count = args.minhash_perm // args.lsh_rows_per_band
    if band_count > 65535:
        raise SystemExit("Too many LSH bands for uint16 band id")

    out_dir = Path(args.output_root) / f"batch_{task['task_id']:06d}"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / f"batch_{task['task_id']:06d}"

    ligand_map_path = str(prefix) + ".ligand_map.tsv.bz2"
    feature_hashes_path = str(prefix) + ".feature_hashes.bin"
    lsh_packed_path = str(prefix) + ".lsh_packed.bin"
    lsh_memberships_path = str(prefix) + ".lsh_memberships.bin" if args.write_lsh_memberships_bin else None
    recipe_path = str(prefix) + ".recipes.tsv.bz2" if args.write_recipes else None
    readable_feature_path = str(prefix) + ".features.tsv.bz2" if args.write_readable_features else None
    readable_lsh_path = str(prefix) + ".lsh_buckets.tsv.bz2" if args.write_readable_lsh else None
    feature_token_map_path = str(prefix) + ".feature_token_map.tsv.bz2" if args.write_feature_token_map else None

    fragment_summary = {} if args.write_summary_csvs else None
    linker_summary = {} if args.write_summary_csvs else None
    connection_summary = {} if args.write_summary_csvs else None
    unique_feature_tokens = {} if args.write_feature_token_map else None

    total_seen = valid_mols = invalid_mols = failed_mols = 0
    local_ligand_id = 0
    feature_count_sum = 0
    feature_count_max = 0

    with ExitStack() as stack:
        ligand_map_handle = stack.enter_context(bz2.open(ligand_map_path, "wt"))
        ligand_map_columns = ligand_map_header(args)
        ligand_map_handle.write("\t".join(ligand_map_columns) + "\n")

        feature_hashes_handle = stack.enter_context(open(feature_hashes_path, "wb"))
        lsh_packed_handle = stack.enter_context(open(lsh_packed_path, "wb"))
        lsh_memberships_handle = stack.enter_context(open(lsh_memberships_path, "wb")) if lsh_memberships_path else None

        recipe_handle = stack.enter_context(bz2.open(recipe_path, "wt")) if recipe_path else None
        if recipe_handle:
            recipe_handle.write("batch_id\tlocal_ligand_id\tsource_file\tsource_data_row\tsmiles\tligand_name\trecipe_json\n")

        readable_feature_handle = stack.enter_context(bz2.open(readable_feature_path, "wt")) if readable_feature_path else None
        if readable_feature_handle:
            readable_feature_handle.write("batch_id\tlocal_ligand_id\tfeature_type\tfeature_key\n")

        readable_lsh_handle = stack.enter_context(bz2.open(readable_lsh_path, "wt")) if readable_lsh_path else None
        if readable_lsh_handle:
            readable_lsh_handle.write("batch_id\tlocal_ligand_id\tlsh_scheme\tband\tbucket_hash_u64\n")

        for item in task["items"]:
            input_path = item["input_path"]
            source_file = os.path.basename(input_path)
            for row_index, smiles, ligand_name in iter_ligands_in_range(input_path, item["start_row"], item["end_row"]):
                total_seen += 1
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    invalid_mols += 1
                    continue
                try:
                    component_info, linker_objects, recipe = classify_and_recipe_molecule(mol, args.linker_max_heavy_atoms)

                    if args.write_summary_csvs:
                        qualified = qualify_ligand(source_file, ligand_name)
                        update_summaries(qualified, recipe, fragment_summary, linker_summary, connection_summary)

                    tokens = make_feature_tokens(recipe)
                    hashes = feature_hashes(tokens, seed=args.feature_hash_seed)
                    buckets = lsh_bucket_hashes(tokens, args.minhash_perm, args.lsh_rows_per_band)

                    ligand_map_handle.write("\t".join(ligand_map_row(local_ligand_id, source_file, row_index, smiles, ligand_name, args)) + "\n")
                    write_feature_hash_record(feature_hashes_handle, local_ligand_id, hashes)
                    write_lsh_packed_record(lsh_packed_handle, local_ligand_id, buckets)
                    if lsh_memberships_handle:
                        write_lsh_membership_records(lsh_memberships_handle, local_ligand_id, buckets)

                    if recipe_handle:
                        recipe_json = json.dumps(recipe, sort_keys=True, separators=(",", ":"))
                        recipe_handle.write(f"{task['task_id']}\t{local_ligand_id}\t{source_file}\t{row_index}\t{smiles}\t{ligand_name}\t{recipe_json}\n")

                    if readable_feature_handle:
                        for token in tokens:
                            readable_feature_handle.write(f"{task['task_id']}\t{local_ligand_id}\t{feature_type_from_token(token)}\t{token}\n")

                    if readable_lsh_handle:
                        scheme = f"minhash_{args.minhash_perm}x{args.lsh_rows_per_band}"
                        for band, bucket_hash in enumerate(buckets):
                            readable_lsh_handle.write(f"{task['task_id']}\t{local_ligand_id}\t{scheme}\t{band}\t{bucket_hash}\n")

                    if unique_feature_tokens is not None:
                        for token in tokens:
                            unique_feature_tokens.setdefault(stable_hash_u64(token, seed=args.feature_hash_seed), token)

                    n_features = len(hashes)
                    feature_count_sum += n_features
                    feature_count_max = max(feature_count_max, n_features)
                    valid_mols += 1
                    local_ligand_id += 1

                except Exception as exc:
                    failed_mols += 1
                    sys.stderr.write(f"Failed {source_file} row {row_index} ligand {ligand_name}: {exc}\n")

    if unique_feature_tokens is not None:
        with bz2.open(feature_token_map_path, "wt") as handle:
            handle.write("feature_hash_u64\tfeature_type\tfeature_token\n")
            for h, token in sorted(unique_feature_tokens.items()):
                handle.write(f"{h}\t{feature_type_from_token(token)}\t{token}\n")

    if args.write_summary_csvs:
        write_summary_csvs(str(prefix), fragment_summary, linker_summary, connection_summary)

    paths = {
        "ligand_map": ligand_map_path,
        "ligand_map_columns": ligand_map_columns,
        "feature_hashes_bin": feature_hashes_path,
        "lsh_packed_bin": lsh_packed_path,
        "lsh_memberships_bin": lsh_memberships_path,
    }
    schema_path = write_binary_schema(prefix, args, band_count, paths)

    stats_path = str(prefix) + ".stats.json"
    stats = {
        "script_version": SCRIPT_VERSION,
        "feature_schema_version": FEATURE_SCHEMA_VERSION,
        "rdkit_version": rdBase.rdkitVersion,
        "task_id": task["task_id"],
        "items": task["items"],
        "total_rows_seen_in_range": total_seen,
        "valid_molecules": valid_mols,
        "invalid_molecules": invalid_mols,
        "failed_molecules": failed_mols,
        "local_ligand_ids_assigned": local_ligand_id,
        "avg_feature_hashes_per_valid_ligand": (feature_count_sum / valid_mols) if valid_mols else 0,
        "max_feature_hashes_per_valid_ligand": feature_count_max,
        "minhash_perm": args.minhash_perm,
        "lsh_rows_per_band": args.lsh_rows_per_band,
        "lsh_band_count": band_count,
        "lsh_record_size_bytes": 4 + 8 * band_count,
        "expected_lsh_packed_size_bytes": valid_mols * (4 + 8 * band_count),
        "ligand_map_path": ligand_map_path,
        "feature_hashes_bin_path": feature_hashes_path,
        "lsh_packed_bin_path": lsh_packed_path,
        "lsh_memberships_bin_path": lsh_memberships_path,
        "schema_path": schema_path,
        "recipe_path": recipe_path,
        "readable_feature_path": readable_feature_path,
        "readable_lsh_path": readable_lsh_path,
        "feature_token_map_path": feature_token_map_path,
        "summary_csvs_written": bool(args.write_summary_csvs),
    }
    try:
        stats["actual_lsh_packed_size_bytes"] = os.path.getsize(lsh_packed_path)
        stats["actual_feature_hashes_size_bytes"] = os.path.getsize(feature_hashes_path)
    except OSError:
        pass
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2, sort_keys=True)

    sys.stderr.write(
        f"Done task {task['task_id']}: valid={valid_mols} invalid={invalid_mols} failed={failed_mols} "
        f"bands={band_count} avg_features={stats['avg_feature_hashes_per_valid_ligand']:.2f}\n"
    )


def write_summary_csvs(prefix, fragment_summary, linker_summary, connection_summary):
    with open(prefix + ".fragments.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "fragment_smiles", "role", "heavy_atoms", "occurrences",
            "connected_fragment_smiles", "connected_fragment_counts",
            "connected_linker_smiles", "connected_linker_counts",
            "connected_terminal_fragment_smiles", "connected_terminal_fragment_counts",
            "connection_pattern_entries", "connection_pattern_counts", "example_ligands",
        ])
        for (smiles, role, heavy), info in sorted(fragment_summary.items(), key=lambda x: (-x[1]["count"], x[0])):
            writer.writerow([
                smiles, role, heavy, info["count"],
                count_keys(info["connected_fragments"]), count_values(info["connected_fragments"]),
                count_keys(info["connected_linkers"]), count_values(info["connected_linkers"]),
                count_keys(info["connected_terminal_fragments"]), count_values(info["connected_terminal_fragments"]),
                count_keys(info["connection_patterns"]), count_values(info["connection_patterns"]),
                ";".join(info["examples"]),
            ])

    with open(prefix + ".linkers.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "linker_smiles", "heavy_atoms", "occurrences",
            "connected_fragment_smiles", "connected_fragment_counts",
            "connected_linker_smiles", "connected_linker_counts",
            "connected_terminal_fragment_smiles", "connected_terminal_fragment_counts",
            "connection_pattern_entries", "connection_pattern_counts", "example_ligands",
        ])
        for (smiles, heavy), info in sorted(linker_summary.items(), key=lambda x: (-x[1]["count"], x[0])):
            writer.writerow([
                smiles, heavy, info["count"],
                count_keys(info["connected_fragments"]), count_values(info["connected_fragments"]),
                count_keys(info["connected_linkers"]), count_values(info["connected_linkers"]),
                count_keys(info["connected_terminal_fragments"]), count_values(info["connected_terminal_fragments"]),
                count_keys(info["connection_patterns"]), count_values(info["connection_patterns"]),
                ";".join(info["examples"]),
            ])

    with open(prefix + ".linker_connections.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["connection_type", "linker_smiles", "connected_fragment_smiles", "occurrences", "example_ligands"])
        for (ctype, linker, connected), info in sorted(connection_summary.items(), key=lambda x: (-x[1]["count"], x[0])):
            writer.writerow([ctype, linker, connected, info["count"], ";".join(info["examples"])])


# ------------------------------- query mode ---------------------------------


def infer_batch_prefix(args):
    if args.batch_prefix:
        return Path(args.batch_prefix)
    if not args.batch_dir:
        raise SystemExit("query mode requires --batch-dir or --batch-prefix")
    batch_dir = Path(args.batch_dir)
    return batch_dir / batch_dir.name


def load_batch_schema(prefix):
    schema_path = Path(str(prefix) + ".binary_schema.json")
    if not schema_path.exists():
        return None
    with open(schema_path) as handle:
        return json.load(handle)


def schema_value(schema, path, default=None):
    if schema is None:
        return default
    current = schema
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def iter_feature_hash_records(path):
    with open(path, "rb") as handle:
        while True:
            header = handle.read(FEATURE_RECORD_HEADER.size)
            if not header:
                break
            if len(header) != FEATURE_RECORD_HEADER.size:
                raise IOError(f"Truncated feature hash header in {path}")
            local_ligand_id, n_features = FEATURE_RECORD_HEADER.unpack(header)
            raw = handle.read(8 * n_features)
            if len(raw) != 8 * n_features:
                raise IOError(f"Truncated feature hash record for ligand {local_ligand_id} in {path}")
            hashes = struct.unpack("<" + "Q" * n_features, raw) if n_features else ()
            yield local_ligand_id, hashes


def iter_lsh_packed_records(path, band_count):
    record = struct.Struct("<I" + "Q" * band_count)
    with open(path, "rb") as handle:
        while True:
            raw = handle.read(record.size)
            if not raw:
                break
            if len(raw) != record.size:
                raise IOError(f"Truncated LSH packed record in {path}")
            values = record.unpack(raw)
            yield values[0], values[1:]


def read_ligand_map_for_ids(path, wanted_ids):
    wanted = set(wanted_ids)
    result = {}
    if not wanted or not path or not Path(path).exists():
        return result
    with bz2.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                lid = int(row["local_ligand_id"])
            except Exception:
                continue
            if lid in wanted:
                result[lid] = row
                if len(result) == len(wanted):
                    break
    return result


def parse_feature_type_filter(value):
    if not value or value.lower() in {"all", "*"}:
        return None
    return {x.strip() for x in value.split(",") if x.strip()}


def decompose_query_smiles(smiles, linker_max_heavy_atoms):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"Could not parse query SMILES: {smiles}")
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    _, _, recipe = classify_and_recipe_molecule(mol, linker_max_heavy_atoms)
    tokens = make_feature_tokens(recipe)
    return canonical_smiles, recipe, tokens


def filter_tokens_by_type(tokens, allowed_types):
    if allowed_types is None:
        return list(tokens)
    return [tok for tok in tokens if feature_type_from_token(tok) in allowed_types]


def heap_add_top(heap, max_size, item):
    # item is ordered by score tuple first. Keep the largest items.
    if max_size is None or max_size <= 0:
        return
    if len(heap) < max_size:
        heappush(heap, item)
    elif item > heap[0]:
        heappushpop(heap, item)


def write_query_recipe(prefix, canonical_smiles, recipe, tokens, hashes, buckets):
    with open(str(prefix) + ".query_recipe.json", "w") as handle:
        json.dump({
            "canonical_query_smiles": canonical_smiles,
            "recipe": recipe,
            "feature_tokens": tokens,
            "feature_hashes_u64": hashes,
            "lsh_bucket_hashes_u64_by_band": buckets,
        }, handle, indent=2, sort_keys=True)


def query_batch(args):
    prefix = infer_batch_prefix(args)
    schema = load_batch_schema(prefix)

    feature_hashes_path = Path(args.feature_hashes_bin or schema_value(schema, ["feature_hashes_bin", "path"], str(prefix) + ".feature_hashes.bin"))
    lsh_packed_path = Path(args.lsh_packed_bin or schema_value(schema, ["lsh_packed_bin", "path"], str(prefix) + ".lsh_packed.bin"))
    ligand_map_path = Path(args.ligand_map or schema_value(schema, ["ligand_map", "path"], str(prefix) + ".ligand_map.tsv.bz2"))

    if not feature_hashes_path.exists():
        raise SystemExit(f"Missing feature hash file: {feature_hashes_path}")
    if not lsh_packed_path.exists():
        sys.stderr.write(f"Warning: missing LSH packed file; LSH matching disabled: {lsh_packed_path}\n")
        lsh_packed_path = None

    linker_max_heavy_atoms = args.linker_max_heavy_atoms
    if linker_max_heavy_atoms is None:
        linker_max_heavy_atoms = schema_value(schema, ["feature_schema", "linker_max_heavy_atoms"], None)
    if linker_max_heavy_atoms is None:
        linker_max_heavy_atoms = DEFAULT_LINKER_MAX_HEAVY_ATOMS

    feature_hash_seed = args.feature_hash_seed
    if feature_hash_seed is None:
        feature_hash_seed = schema_value(schema, ["feature_hashes_bin", "feature_hash_seed"], 0)

    minhash_perm = args.minhash_perm
    if minhash_perm is None:
        minhash_perm = schema_value(schema, ["lsh_packed_bin", "minhash_perm"], DEFAULT_MINHASH_PERM)

    lsh_rows_per_band = args.lsh_rows_per_band
    if lsh_rows_per_band is None:
        lsh_rows_per_band = schema_value(schema, ["lsh_packed_bin", "lsh_rows_per_band"], DEFAULT_LSH_ROWS_PER_BAND)

    band_count = minhash_perm // lsh_rows_per_band
    schema_band_count = schema_value(schema, ["lsh_packed_bin", "band_count"], None)
    if schema_band_count is not None:
        band_count = int(schema_band_count)

    canonical_smiles, recipe, all_tokens = decompose_query_smiles(args.smiles, linker_max_heavy_atoms)
    allowed_types = parse_feature_type_filter(args.match_feature_types)
    query_tokens = filter_tokens_by_type(all_tokens, allowed_types)
    query_hash_by_token = {tok: stable_hash_u64(tok, seed=feature_hash_seed) for tok in query_tokens}
    query_token_by_hash = defaultdict(list)
    for tok, h in query_hash_by_token.items():
        query_token_by_hash[h].append(tok)
    query_hash_set = set(query_token_by_hash)
    query_buckets = lsh_bucket_hashes(query_tokens, minhash_perm, lsh_rows_per_band)

    out_prefix = Path(args.output_prefix) if args.output_prefix else Path(str(prefix) + ".query")
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    write_query_recipe(out_prefix, canonical_smiles, recipe, query_tokens, sorted(query_hash_set), query_buckets)

    # Scan feature hashes: count per query feature and rank ligands by shared query hashes.
    feature_presence_counts = defaultdict(int)
    top_feature_heap = []
    all_feature_matches = [] if args.write_all_feature_matches else None
    total_records = 0
    ligands_with_any_shared_feature = 0

    for local_ligand_id, hashes in iter_feature_hash_records(feature_hashes_path):
        total_records += 1
        shared = query_hash_set.intersection(hashes)
        if not shared:
            continue
        ligands_with_any_shared_feature += 1
        for h in shared:
            feature_presence_counts[h] += 1
        shared_count = len(shared)
        score_tuple = (shared_count, local_ligand_id, tuple(sorted(shared)))
        heap_add_top(top_feature_heap, args.top_n, score_tuple)
        if all_feature_matches is not None:
            all_feature_matches.append(score_tuple)

    top_feature_matches = sorted(nlargest(args.top_n, top_feature_heap), reverse=True)

    # Scan LSH packed: count shared band-specific buckets and rank ligands.
    top_lsh_heap = []
    lsh_counts_for_top_feature_ligands = {}
    top_feature_ligand_ids = {x[1] for x in top_feature_matches}
    lsh_ligands_with_any_shared_bucket = 0

    if lsh_packed_path is not None:
        for local_ligand_id, buckets in iter_lsh_packed_records(lsh_packed_path, band_count):
            shared_bands = [band for band, bucket in enumerate(buckets) if band < len(query_buckets) and bucket == query_buckets[band]]
            if not shared_bands:
                continue
            lsh_ligands_with_any_shared_bucket += 1
            shared_band_count = len(shared_bands)
            heap_add_top(top_lsh_heap, args.top_n, (shared_band_count, local_ligand_id, tuple(shared_bands)))
            if local_ligand_id in top_feature_ligand_ids:
                lsh_counts_for_top_feature_ligands[local_ligand_id] = shared_band_count

    top_lsh_matches = sorted(nlargest(args.top_n, top_lsh_heap), reverse=True)

    # Resolve selected ligand IDs to source metadata.
    selected_ids = {x[1] for x in top_feature_matches} | {x[1] for x in top_lsh_matches}
    ligand_map = read_ligand_map_for_ids(ligand_map_path, selected_ids)

    # Write query components.
    components_path = str(out_prefix) + ".query_components.tsv"
    with open(components_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["kind", "local_id", "role", "smiles", "heavy_atoms", "degree", "connected_fragment_smiles"])
        for frag in recipe["fragments"]:
            writer.writerow(["fragment_or_terminal", frag["fragment_id"], frag["role"], frag["smiles"], frag["heavy_atoms"], frag["degree"], ""])
        for linker in recipe["linkers"]:
            writer.writerow(["linker", linker["linker_id"], "linker", linker["smiles"], linker["heavy_atoms"], "", ";".join(linker.get("connected_fragment_smiles", []))])

    # Write query feature presence counts.
    presence_path = str(out_prefix) + ".query_feature_presence.tsv"
    with open(presence_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature_type", "feature_token", "feature_hash_u64", "batch_ligands_containing_feature", "present_in_batch"])
        for token, h in sorted(query_hash_by_token.items(), key=lambda x: (feature_type_from_token(x[0]), x[0])):
            count = feature_presence_counts.get(h, 0)
            writer.writerow([feature_type_from_token(token), token, h, count, int(count > 0)])

    def ligand_metadata_columns():
        # Stable set; missing ligand map fields are left blank.
        return ["source_file", "source_data_row", "ligand_name", "smiles"]

    def ligand_metadata_values(local_ligand_id):
        row = ligand_map.get(local_ligand_id, {})
        return [row.get("source_file", ""), row.get("source_data_row", ""), row.get("ligand_name", ""), row.get("smiles", "")]

    # Write top feature matches.
    feature_matches_path = str(out_prefix) + ".top_feature_matches.tsv"
    with open(feature_matches_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "local_ligand_id", "shared_query_feature_count", "shared_lsh_band_count", "shared_feature_hashes_u64", "shared_feature_tokens"] + ligand_metadata_columns())
        for rank, (shared_count, local_ligand_id, shared_hashes) in enumerate(top_feature_matches, start=1):
            shared_tokens = []
            for h in shared_hashes:
                shared_tokens.extend(query_token_by_hash[h])
            writer.writerow([
                rank,
                local_ligand_id,
                shared_count,
                lsh_counts_for_top_feature_ligands.get(local_ligand_id, 0),
                ";".join(str(x) for x in shared_hashes),
                ";".join(sorted(shared_tokens)),
            ] + ligand_metadata_values(local_ligand_id))

    # Optional all feature matches. This can be large.
    all_feature_matches_path = None
    if all_feature_matches is not None:
        all_feature_matches_path = str(out_prefix) + ".all_feature_matches.tsv"
        all_feature_matches_sorted = sorted(all_feature_matches, reverse=True)
        with open(all_feature_matches_path, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["local_ligand_id", "shared_query_feature_count", "shared_feature_hashes_u64", "shared_feature_tokens"])
            for shared_count, local_ligand_id, shared_hashes in all_feature_matches_sorted:
                shared_tokens = []
                for h in shared_hashes:
                    shared_tokens.extend(query_token_by_hash[h])
                writer.writerow([local_ligand_id, shared_count, ";".join(str(x) for x in shared_hashes), ";".join(sorted(shared_tokens))])

    # Write top LSH matches.
    lsh_matches_path = str(out_prefix) + ".top_lsh_matches.tsv"
    with open(lsh_matches_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "local_ligand_id", "shared_lsh_band_count", "shared_lsh_bands"] + ligand_metadata_columns())
        for rank, (shared_band_count, local_ligand_id, shared_bands) in enumerate(top_lsh_matches, start=1):
            writer.writerow([rank, local_ligand_id, shared_band_count, ";".join(str(x) for x in shared_bands)] + ligand_metadata_values(local_ligand_id))

    # Write summary JSON.
    summary_path = str(out_prefix) + ".query_summary.json"
    summary = {
        "query_smiles_input": args.smiles,
        "canonical_query_smiles": canonical_smiles,
        "batch_prefix": str(prefix),
        "feature_hashes_bin": str(feature_hashes_path),
        "lsh_packed_bin": str(lsh_packed_path) if lsh_packed_path else None,
        "ligand_map": str(ligand_map_path),
        "linker_max_heavy_atoms": linker_max_heavy_atoms,
        "feature_hash_seed": feature_hash_seed,
        "minhash_perm": minhash_perm,
        "lsh_rows_per_band": lsh_rows_per_band,
        "lsh_band_count": band_count,
        "match_feature_types": sorted(allowed_types) if allowed_types else "all",
        "query_feature_count": len(query_tokens),
        "query_feature_hash_count": len(query_hash_set),
        "batch_feature_records_scanned": total_records,
        "batch_ligands_with_any_shared_query_feature": ligands_with_any_shared_feature,
        "batch_ligands_with_any_shared_lsh_bucket": lsh_ligands_with_any_shared_bucket,
        "outputs": {
            "query_recipe_json": str(out_prefix) + ".query_recipe.json",
            "query_components_tsv": components_path,
            "query_feature_presence_tsv": presence_path,
            "top_feature_matches_tsv": feature_matches_path,
            "top_lsh_matches_tsv": lsh_matches_path,
            "all_feature_matches_tsv": all_feature_matches_path,
        },
    }
    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    sys.stderr.write("Query complete.\n")
    sys.stderr.write(f"Feature records scanned: {total_records}\n")
    sys.stderr.write(f"Ligands with any shared query feature: {ligands_with_any_shared_feature}\n")
    sys.stderr.write(f"Ligands with any shared LSH bucket: {lsh_ligands_with_any_shared_bucket}\n")
    sys.stderr.write(f"Query feature presence: {presence_path}\n")
    sys.stderr.write(f"Top feature matches: {feature_matches_path}\n")
    sys.stderr.write(f"Top LSH matches: {lsh_matches_path}\n")


# --------------------------- manifest and LSF mode --------------------------


def discover_bz2_files(input_dir, pattern):
    root = Path(input_dir)
    return sorted(str(p) for p in root.rglob(pattern) if p.is_file())


def make_manifest(files, output_root, batch_size, whole_files=False):
    manifest_path = Path(output_root) / "manifest.tsv"
    count_path = Path(output_root) / "input_file_counts.tsv"
    Path(output_root).mkdir(parents=True, exist_ok=True)

    file_counts = []
    for i, path in enumerate(files, start=1):
        n = count_data_rows(path)
        file_counts.append((path, n))
        sys.stderr.write(f"Counted {i}/{len(files)}: {path} rows={n}\n")

    with open(count_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["input_path", "data_rows"])
        for path, n in file_counts:
            writer.writerow([path, n])

    tasks = []
    current = []
    current_n = 0

    if whole_files:
        for path, n in file_counts:
            if current and current_n + n > batch_size:
                tasks.append(current)
                current = []
                current_n = 0
            current.append((path, 0, n, n))
            current_n += n
        if current:
            tasks.append(current)
    else:
        for path, n in file_counts:
            start = 0
            while start < n:
                if current_n >= batch_size:
                    tasks.append(current)
                    current = []
                    current_n = 0
                take = min(batch_size - current_n, n - start)
                current.append((path, start, start + take, take))
                current_n += take
                start += take
        if current:
            tasks.append(current)

    with open(manifest_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["task_id", "input_path", "start_row", "end_row", "num_ligands"])
        for task_id, task_items in enumerate(tasks, start=1):
            for path, start, end, n in task_items:
                writer.writerow([task_id, path, start, end, n])

    return str(manifest_path), str(count_path), len(tasks), sum(n for path, n in file_counts)


def read_manifest_task(manifest_path, task_id):
    items = []
    with open(manifest_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["task_id"]) == task_id:
                items.append({
                    "input_path": row["input_path"],
                    "start_row": int(row["start_row"]),
                    "end_row": int(row["end_row"]),
                    "num_ligands": int(row["num_ligands"]),
                })
    if not items:
        raise ValueError(f"No manifest entries found for task_id {task_id}")
    return {"task_id": task_id, "items": items}


def add_worker_flag(worker_cmd, args, flag_name, attr_name=None):
    attr = attr_name or flag_name.lstrip("-").replace("-", "_")
    if getattr(args, attr):
        worker_cmd.append(flag_name)


def build_bsub_command(args, manifest_path, num_tasks):
    logs = Path(args.output_root) / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    script_path = os.path.abspath(sys.argv[0])
    array_spec = f"recipe_feat[1-{num_tasks}]"
    if args.array_throttle:
        array_spec += f"%{args.array_throttle}"

    worker_cmd = [
        sys.executable,
        script_path,
        "worker",
        "--manifest", manifest_path,
        "--task-id", "$LSB_JOBINDEX",
        "--output-root", args.output_root,
        "--linker-max-heavy-atoms", str(args.linker_max_heavy_atoms),
        "--minhash-perm", str(args.minhash_perm),
        "--lsh-rows-per-band", str(args.lsh_rows_per_band),
        "--feature-hash-seed", str(args.feature_hash_seed),
    ]
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
        add_worker_flag(worker_cmd, args, flag)

    quoted_worker_cmd = " ".join(shlex.quote(x) if x != "$LSB_JOBINDEX" else x for x in worker_cmd)

    bsub_cmd = [
        "bsub",
        "-q", args.queue,
        "-J", array_spec,
        "-oo", str(logs / "%J_%I.out"),
        "-eo", str(logs / "%J_%I.err"),
    ]
    if args.memory_mb:
        bsub_cmd += ["-M", str(args.memory_mb), "-R", f"rusage[mem={args.memory_mb}]"]
    if args.walltime:
        bsub_cmd += ["-W", args.walltime]
    bsub_cmd.append(quoted_worker_cmd)
    return bsub_cmd


# ---------------------------------- main ------------------------------------


def submit_main(args):
    files = discover_bz2_files(args.input_dir, args.pattern)
    if not files:
        raise SystemExit(f"No files matching {args.pattern} found under {args.input_dir}")

    manifest_path, count_path, num_tasks, total_ligands = make_manifest(
        files,
        args.output_root,
        args.batch_size,
        whole_files=args.whole_files,
    )
    sys.stderr.write(f"Manifest: {manifest_path}\n")
    sys.stderr.write(f"Input counts: {count_path}\n")
    sys.stderr.write(f"Total ligands: {total_ligands}\n")
    sys.stderr.write(f"Tasks: {num_tasks}\n")

    cmd = build_bsub_command(args, manifest_path, num_tasks)
    print(" ".join(shlex.quote(x) for x in cmd))
    if args.submit:
        subprocess.run(cmd, check=True)


def worker_main(args):
    task_id = args.task_id
    if task_id is None:
        env_id = os.environ.get("LSB_JOBINDEX")
        if not env_id:
            raise SystemExit("worker mode needs --task-id or LSB_JOBINDEX")
        task_id = int(env_id)
    task = read_manifest_task(args.manifest, task_id)
    process_manifest_task(task, args)


def query_main(args):
    query_batch(args)


def add_common_worker_args(parser):
    parser.add_argument("--linker-max-heavy-atoms", type=int, default=DEFAULT_LINKER_MAX_HEAVY_ATOMS)
    parser.add_argument("--minhash-perm", type=int, default=DEFAULT_MINHASH_PERM,
                        help="Number of MinHash permutations. Default: %(default)s")
    parser.add_argument("--lsh-rows-per-band", type=int, default=DEFAULT_LSH_ROWS_PER_BAND,
                        help="Rows per LSH band. Band count = minhash_perm / rows_per_band. Default: %(default)s")
    parser.add_argument("--feature-hash-seed", type=int, default=0)
    parser.add_argument("--map-include-ligand-name", action="store_true",
                        help="Add ligand_name to ligand_map.tsv.bz2. Off by default to reduce output size.")
    parser.add_argument("--map-include-smiles", action="store_true",
                        help="Add original SMILES to ligand_map.tsv.bz2. Off by default to reduce output size.")
    parser.add_argument("--write-recipes", action="store_true",
                        help="Write full readable recipe JSON TSV. Useful for debugging; off by default.")
    parser.add_argument("--write-readable-features", action="store_true",
                        help="Write readable feature-key TSV. Useful for debugging; off by default.")
    parser.add_argument("--write-readable-lsh", action="store_true",
                        help="Write readable LSH bucket TSV. Bulky; off by default.")
    parser.add_argument("--write-summary-csvs", action="store_true",
                        help="Write fragments/linkers/linker_connections summary CSVs. Off by default.")
    parser.add_argument("--write-feature-token-map", action="store_true",
                        help="Write per-batch map from feature_hash_u64 to feature token. Useful for interpretation/debugging; off by default.")
    parser.add_argument("--write-lsh-memberships-bin", action="store_true",
                        help="Also write unfolded binary records band,bucket_hash,ligand_id. Packed LSH can be unfolded later, so off by default.")


def main():
    parser = argparse.ArgumentParser(description="Recipe/feature/compact-LSH extraction and query pipeline for Enamine bz2 files.")
    sub = parser.add_subparsers(dest="mode", required=True)

    p_submit = sub.add_parser("submit", help="Discover inputs, create manifest, and print/submit an LSF array job.")
    p_submit.add_argument("input_dir", help="Directory to recursively search for input .bz2 files.")
    p_submit.add_argument("--output-root", required=True, help="Directory for manifest and batch outputs.")
    p_submit.add_argument("--pattern", default="*.bz2", help="Recursive filename pattern. Default: %(default)s")
    p_submit.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Target ligands per worker task. Default: %(default)s")
    p_submit.add_argument("--whole-files", action="store_true", help="Do not split bz2 files by row range; group whole files into batches instead.")
    p_submit.add_argument("--queue", default="long", help="LSF queue. Default: %(default)s")
    p_submit.add_argument("--walltime", default=None, help="Optional LSF -W walltime, e.g. 12:00 or 24:00.")
    p_submit.add_argument("--memory-mb", type=int, default=None, help="Optional LSF memory request in MB.")
    p_submit.add_argument("--array-throttle", type=int, default=1500, help="LSF array throttle. Default: %(default)s")
    p_submit.add_argument("--submit", action="store_true", help="Actually run bsub. Without this, only print the command.")
    add_common_worker_args(p_submit)
    p_submit.set_defaults(func=submit_main)

    p_worker = sub.add_parser("worker", help="Process one manifest task. Usually called by LSF array.")
    p_worker.add_argument("--manifest", required=True)
    p_worker.add_argument("--task-id", type=int, default=None)
    p_worker.add_argument("--output-root", required=True)
    add_common_worker_args(p_worker)
    p_worker.set_defaults(func=worker_main)

    p_query = sub.add_parser("query", help="Query one compact batch by decomposing a SMILES and scanning binary feature/LSH files.")
    p_query.add_argument("--smiles", required=True, help="Query SMILES to decompose into fragments/linkers/features.")
    p_query.add_argument("--batch-dir", default=None, help="Batch directory, e.g. /out/batch_000046. Prefix inferred as /out/batch_000046/batch_000046.")
    p_query.add_argument("--batch-prefix", default=None, help="Explicit batch prefix, e.g. /out/batch_000046/batch_000046.")
    p_query.add_argument("--feature-hashes-bin", default=None, help="Optional explicit path to feature_hashes.bin.")
    p_query.add_argument("--lsh-packed-bin", default=None, help="Optional explicit path to lsh_packed.bin.")
    p_query.add_argument("--ligand-map", default=None, help="Optional explicit path to ligand_map.tsv.bz2.")
    p_query.add_argument("--output-prefix", default=None, help="Output prefix for query result files. Default: <batch_prefix>.query")
    p_query.add_argument("--top-n", type=int, default=100, help="Number of top feature/LSH matches to write. Default: %(default)s")
    p_query.add_argument("--write-all-feature-matches", action="store_true", help="Write every ligand sharing at least one query feature. Can be large.")
    p_query.add_argument("--match-feature-types", default="all",
                         help="Comma-separated feature types used for matching, e.g. FRAG,LINKER,TERM,MOTIF. Default: all")
    p_query.add_argument("--linker-max-heavy-atoms", type=int, default=None,
                         help="Query linker cutoff. Default from batch schema if available, else 4.")
    p_query.add_argument("--minhash-perm", type=int, default=None,
                         help="Query MinHash permutations. Default from batch schema if available, else 64.")
    p_query.add_argument("--lsh-rows-per-band", type=int, default=None,
                         help="Query rows per LSH band. Default from batch schema if available, else 8.")
    p_query.add_argument("--feature-hash-seed", type=int, default=None,
                         help="Feature hash seed. Default from batch schema if available, else 0.")
    p_query.set_defaults(func=query_main)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
