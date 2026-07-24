#!/usr/bin/env python3
"""Build manifests from a folder of cleaned PDBs.

Writes a CSV manifest (header + one row per experiment).

Usage:
    python3 build_manifests.py --pdb-dir PATH --out manifests.csv [--dry-run] [--limit N]
"""
import argparse
import glob
import csv
import json
import os
import sys


def find_pdbs(pdb_dir):
    pattern = os.path.join(pdb_dir, "*.clean.pdb")
    return sorted(glob.glob(pattern))


def build_baseline_row(source_pdb, manifest_id):
    return [
        manifest_id,
        os.path.abspath(source_pdb),
        "baseline",
        "",
        "",
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir", required=True, help="Directory with .clean.pdb files")
    parser.add_argument("--out", default="manifests.jsonl", help="Output JSONL manifest path")
    parser.add_argument("--dry-run", action="store_true", help="Do not call convex_hull; write baseline rows only")
    parser.add_argument("--limit", type=int, default=0, help="Limit number of PDBs processed (0 = all)")
    args = parser.parse_args()

    pdbs = find_pdbs(args.pdb_dir)
    if args.limit and args.limit > 0:
        pdbs = pdbs[: args.limit]

    if not pdbs:
        print(f"No .clean.pdb files found in {args.pdb_dir}")
        sys.exit(0)

    failures = []
    manifest_id = 0
    out_dir = os.path.dirname(os.path.abspath(args.out)) or "."
    os.makedirs(out_dir, exist_ok=True)

    header = [
        "manifest_id",
        "source_pdb",
        "mode",
        "protein_flex_residue",
        "ligand_flex_residue",
    ]

    with open(args.out, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)

        for pdb in pdbs:
            try:
                if args.dry_run:
                    row = build_baseline_row(pdb, manifest_id)
                    writer.writerow(row)
                    manifest_id += 1
                    continue

                # non-dry run: import and call convex_hull.run_convex_hull
                here = os.path.dirname(os.path.abspath(__file__))
                if here not in sys.path:
                    sys.path.insert(0, here)

                import convex_hull

                contacts, interchain = convex_hull.run_convex_hull(pdb, outfolder="pdb_hulls")
                rows = convex_hull.build_experiment_rows(interchain)
                for r in rows:
                    r_out = [
                        manifest_id,
                        os.path.abspath(pdb),
                        r.get("mode", "baseline"),
                        r.get("protein_flex_residue", ""),
                        r.get("ligand_flex_residue", ""),
                    ]
                    writer.writerow(r_out)
                    manifest_id += 1

            except Exception as e:
                failures.append({"pdb": pdb, "error": str(e)})
                print(f"Failed processing {pdb}: {e}")

    if failures:
        fail_path = os.path.join(out_dir, "manifests_failures.jsonl")
        with open(fail_path, "w") as ff:
            for f in failures:
                ff.write(json.dumps(f) + "\n")
        print(f"Wrote {len(failures)} failure(s) to {fail_path}")

    print(f"Wrote CSV manifest to {args.out} ({manifest_id} rows)")


if __name__ == "__main__":
    main()
