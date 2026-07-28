#!/usr/bin/env python3
import csv
import os
import sys

from convex_hull import run_convex_hull, build_experiment_rows


def write_manifest(pdb_path, manifest_path):
    try:
        _, interchain = run_convex_hull(pdb_path)
    except BaseException as exc:
        print(f"Skipping {pdb_path}: convex-hull failed ({type(exc).__name__}: {exc})")
        return None

    rows = build_experiment_rows(interchain)

    manifest_dir = os.path.dirname(manifest_path)
    if manifest_dir:
        os.makedirs(manifest_dir, exist_ok=True)

    with open(manifest_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb", "mode", "protein_flex_residue", "ligand_flex_residue"])
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "pdb": pdb_path,
                **row,
            })

    print(f"Wrote {len(rows)} rows to {manifest_path}")
    return rows


def main():
    if len(sys.argv) < 3:
        print("Usage: python build_ias_manifest.py <pdb_path> <manifest_path>")
        sys.exit(1)

    pdb_path = sys.argv[1]
    manifest_path = sys.argv[2]
    write_manifest(pdb_path, manifest_path)


if __name__ == "__main__":
    main()
