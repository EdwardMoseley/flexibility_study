#!/usr/bin/env python3
import argparse
import csv
import os
import sys
from pathlib import Path

import pandas as pd

from build_ias_manifest import write_manifest


def resolve_pdb_path(pdb_id, base_dir):
    candidates = [
        Path(pdb_id),
        Path(f"{pdb_id}.pdb"),
        Path(f"{pdb_id}.clean.pdb"),
        Path(f"{pdb_id}.PDB"),
        Path(base_dir) / pdb_id,
        Path(base_dir) / f"{pdb_id}.pdb",
        Path(base_dir) / f"{pdb_id}.clean.pdb",
        Path(base_dir) / f"{pdb_id}.PDB",
        Path(base_dir) / "pdbs_for_analysis" / f"{pdb_id}.pdb",
        Path(base_dir) / "PDBs" / f"{pdb_id}.pdb",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def build_all_manifests(csv_path, output_dir, base_dir=None, limit=None):
    base_dir = Path(base_dir or os.getcwd())
    df = pd.read_csv(csv_path)
    pdb_ids = [str(x).strip() for x in df["pdb"].dropna().astype(str).unique() if str(x).strip()]

    if limit is not None:
        pdb_ids = pdb_ids[:limit]

    os.makedirs(output_dir, exist_ok=True)
    written = []

    for pdb_id in pdb_ids:
        pdb_path = resolve_pdb_path(pdb_id, base_dir)
        if not pdb_path:
            print(f"Skipping {pdb_id}: no local PDB file found")
            continue

        manifest_path = Path(output_dir) / f"{pdb_id}.manifest.csv"
        print(f"Processing {pdb_id} -> {pdb_path}")
        rows = write_manifest(pdb_path, str(manifest_path))
        if rows is None:
            print(f"No manifest generated for {pdb_id}; continuing")
            continue
        written.append((pdb_id, manifest_path, len(rows)))

    return written


def merge_manifest_files(manifest_paths, output_path):
    manifest_paths = [Path(path) for path in manifest_paths if Path(path).exists()]
    fieldnames = ["pdb", "mode", "protein_flex_residue", "ligand_flex_residue"]

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for manifest_path in manifest_paths:
            with open(manifest_path, newline="") as source:
                reader = csv.DictReader(source)
                for row in reader:
                    writer.writerow({key: row.get(key, "") for key in fieldnames})

    if not manifest_paths:
        print(f"No manifest files were found to merge; wrote empty manifest to {output_path}")
    else:
        print(f"Merged {len(manifest_paths)} manifest files into {output_path}")

    return output_path


def main():
    parser = argparse.ArgumentParser(description="Generate IAS manifests from a SKEMPI2 CSV")
    parser.add_argument("csv_path")
    parser.add_argument("output_dir")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--base-dir", default=None)
    parser.add_argument("--aggregate-manifest", default=None)
    args = parser.parse_args()

    written = build_all_manifests(
        args.csv_path,
        args.output_dir,
        base_dir=args.base_dir or Path(args.csv_path).parent,
        limit=args.limit,
    )

    if args.aggregate_manifest:
        manifest_paths = [str(path) for _, path, _ in written]
        merged_path = merge_manifest_files(manifest_paths, args.aggregate_manifest)
        print(f"Wrote merged manifest to {merged_path}")


if __name__ == "__main__":
    main()
