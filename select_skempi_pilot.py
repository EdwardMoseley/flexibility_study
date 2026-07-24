#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def pick_pdb(rows, chain, min_mutations, exclude_pdbs, explicit_pdb=None):
    grouped = defaultdict(set)

    for row in rows:
        pdb = (row.get("pdb") or "").strip()
        mut = (row.get("Mutation.s._cleaned") or "").strip()

        if not pdb or not mut or len(mut) < 3:
            continue
        if mut[1].upper() != chain.upper():
            continue
        if pdb in exclude_pdbs:
            continue

        grouped[pdb].add(mut)

    if explicit_pdb:
        muts = sorted(grouped.get(explicit_pdb, set()))
        if len(muts) < min_mutations:
            raise ValueError(
                f"Requested pdb {explicit_pdb} has {len(muts)} matching mutations; need at least {min_mutations}."
            )
        return explicit_pdb, muts

    for pdb in sorted(grouped):
        muts = sorted(grouped[pdb])
        if len(muts) >= min_mutations:
            return pdb, muts

    raise ValueError("No PDB met selection criteria.")


def main():
    parser = argparse.ArgumentParser(description="Select one SKEMPI PDB and emit a pilot subset CSV")
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    parser.add_argument("--chain", default="B", help="Mutation chain letter from Mutation.s._cleaned (default: B)")
    parser.add_argument("--min-mutations", type=int, default=2)
    parser.add_argument("--mutations-per-pdb", type=int, default=2)
    parser.add_argument("--pdb", default=None, help="Force a specific PDB ID")
    parser.add_argument("--exclude-pdb", action="append", default=[])
    args = parser.parse_args()

    input_csv = Path(args.input_csv)
    output_csv = Path(args.output_csv)

    with open(input_csv, newline="") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = list(rows[0].keys()) if rows else []

    if not rows or not fieldnames:
        raise ValueError(f"No rows found in {input_csv}")

    selected_pdb, selected_mutations = pick_pdb(
        rows,
        chain=args.chain,
        min_mutations=args.min_mutations,
        exclude_pdbs=set(args.exclude_pdb),
        explicit_pdb=args.pdb,
    )

    keep_mutations = set(selected_mutations[: args.mutations_per_pdb])
    subset = [
        row for row in rows
        if (row.get("pdb") or "").strip() == selected_pdb
        and (row.get("Mutation.s._cleaned") or "").strip() in keep_mutations
    ]

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(subset)

    print(f"selected_pdb={selected_pdb}")
    print(f"selected_mutations={sorted(keep_mutations)}")
    print(f"rows_written={len(subset)}")
    print(f"output_csv={output_csv}")


if __name__ == "__main__":
    main()
