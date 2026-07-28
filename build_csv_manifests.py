#!/usr/bin/env python3
import argparse
import ast
import csv
import os
import subprocess
import sys
from pathlib import Path
import re
import shutil


def resolve_pdb_path(pdb_id, base_dir, cleaned_pdb_dir=None):
    cleaned_base = Path(cleaned_pdb_dir) if cleaned_pdb_dir else None
    candidates = [
        (cleaned_base / f"{pdb_id}.clean.pdb") if cleaned_base else None,
        (cleaned_base / f"{pdb_id}.pdb") if cleaned_base else None,
        Path(pdb_id),
        Path(f"{pdb_id}.pdb"),
        Path(f"{pdb_id}.clean.pdb"),
        Path(f"{pdb_id}.PDB"),
        Path(base_dir) / pdb_id,
        Path(base_dir) / f"{pdb_id}.pdb",
        Path(base_dir) / f"{pdb_id}.clean.pdb",
        Path(base_dir) / f"{pdb_id}.PDB",
        Path(base_dir) / "pdbs_for_analysis" / f"{pdb_id}.pdb",
        Path(base_dir) / "cleaned_pdbs_for_analysis" / f"{pdb_id}.clean.pdb",
        Path(base_dir) / "PDBs" / f"{pdb_id}.pdb",
    ]
    for candidate in candidates:
        if candidate and candidate.exists():
            return str(candidate)
    return None


def normalize_truthy(value):
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_mutation_chain(row):
    cleaned = str(row.get("Mutation.s._cleaned", "")).strip()
    if len(cleaned) >= 2 and cleaned[1].isalpha():
        return cleaned[1].upper()
    return str(row.get("Chain", "")).strip().upper()


def parse_mutation_residue_number(row):
    residue = str(row.get("Residue", "")).strip()
    if residue:
        try:
            return int(float(residue))
        except ValueError:
            pass

    cleaned = str(row.get("Mutation.s._cleaned", "")).strip()
    match = re.match(r"^[A-Za-z][A-Za-z](-?\d+)[A-Za-z]$", cleaned)
    if match:
        return int(match.group(1))

    return None


def make_manifest_rows(pdb_path, row, nearby_protein_residues):
    chain = parse_mutation_chain(row)
    residue_number = parse_mutation_residue_number(row)
    if residue_number is None:
        return []

    mutation_residue = f"{chain}{residue_number}"
    common = {
        "pdb": pdb_path,
        "mutation_residue": mutation_residue,
        "skempi_pdb": str(row.get("pdb", "")).strip(),
        "skempi_mutation": str(row.get("Mutation.s._cleaned", "")).strip(),
        "wt_aa": str(row.get("WT_AA", "")).strip(),
        "mut_aa": str(row.get("Mut_AA", "")).strip(),
    }

    rows = [
        {
            **common,
            "mode": "baseline-no-flex",
            "protein_flex_residue": "",
            "ligand_flex_residue": "",
        },
        {
            **common,
            "mode": "ligand-flex",
            "protein_flex_residue": "",
            "ligand_flex_residue": mutation_residue,
        },
    ]

    seen = set()
    for residue in nearby_protein_residues:
        try:
            residue_number = int(residue)
        except (TypeError, ValueError):
            continue

        protein_residue = f"A{residue_number}"
        if protein_residue in seen:
            continue
        seen.add(protein_residue)
        rows.append({
            **common,
            "mode": "ligand-and-single-protein-flex",
            "protein_flex_residue": protein_residue,
            "ligand_flex_residue": mutation_residue,
        })

    return rows


def parse_convex_hull_neighbor_map(stdout_text, design_chain):
    pattern = re.compile(
        rf"^Chain\s+{re.escape(design_chain)}\s+Residue\s+(-?\d+)\s+intersects\s+with\s+Chain\s+\S+\s+residue\(s\)\s+(\[.*\])$"
    )

    residue_map = {}
    for line in stdout_text.splitlines():
        line = line.strip()
        match = pattern.match(line)
        if not match:
            continue

        residue_number = int(match.group(1))
        try:
            neighbors = ast.literal_eval(match.group(2))
        except (ValueError, SyntaxError):
            neighbors = []

        if isinstance(neighbors, list):
            residue_map[residue_number] = neighbors

    return residue_map


def mutation_residue_label(chain, residue_number):
    return f"{chain}{residue_number}"


def results_dir_for_entry(base_dir, pdb_path, mutation_residue):
    pdb_base = Path(pdb_path).name
    if pdb_base.lower().endswith(".pdb"):
        pdb_base = pdb_base[:-4]
    safe_mutation = re.sub(r"[^A-Za-z0-9_.-]", "_", mutation_residue)
    return Path(base_dir) / "results" / pdb_base / safe_mutation


def write_convex_hull_logs(base_dir, pdb_path, mutation_residue, stdout_text, stderr_text):
    out_dir = results_dir_for_entry(base_dir, pdb_path, mutation_residue)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "convex_hull.out").write_text(stdout_text or "", encoding="utf-8")
    (out_dir / "convex_hull.err").write_text(stderr_text or "", encoding="utf-8")


def copy_manifest_to_results(base_dir, pdb_path, mutation_residue, manifest_path):
    out_dir = results_dir_for_entry(base_dir, pdb_path, mutation_residue)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(manifest_path)
    shutil.copyfile(manifest_path, out_dir / manifest_path.name)
    shutil.copyfile(manifest_path, out_dir / "source_manifest.csv")


def run_convex_hull_neighbor_map(pdb_path, design_chain):
    script_path = Path(__file__).with_name("convex_hull.py")
    cmd = [sys.executable, str(script_path), str(pdb_path), str(design_chain)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "convex_hull.py failed")

    residue_map = parse_convex_hull_neighbor_map(result.stdout, design_chain)
    if not residue_map:
        raise RuntimeError("No per-residue neighbor lists were parsed from convex_hull.py output")

    return residue_map, result.stdout, result.stderr


def build_all_manifests(csv_path, output_dir, base_dir=None, limit=None, cleaned_pdb_dir=None):
    base_dir = Path(base_dir or os.getcwd())
    csv_path = Path(csv_path)
    cleaned_pdb_dir = Path(cleaned_pdb_dir) if cleaned_pdb_dir else (base_dir / "cleaned_pdbs_for_analysis")

    with open(csv_path, newline="") as handle:
        all_rows = list(csv.DictReader(handle))

    single_mutant_rows = [row for row in all_rows if normalize_truthy(row.get("is_single_mutant", "1"))]
    if not single_mutant_rows:
        single_mutant_rows = all_rows

    if limit is not None:
        single_mutant_rows = single_mutant_rows[:limit]

    os.makedirs(output_dir, exist_ok=True)
    fieldnames = [
        "pdb",
        "mode",
        "mutation_residue",
        "protein_flex_residue",
        "ligand_flex_residue",
        "skempi_pdb",
        "skempi_mutation",
        "wt_aa",
        "mut_aa",
    ]

    written = []
    cache = {}

    for row in single_mutant_rows:
        pdb_id = str(row.get("pdb", "")).strip()
        if not pdb_id:
            continue

        mutation_chain = parse_mutation_chain(row)
        mutation_residue_number = parse_mutation_residue_number(row)
        if mutation_chain != "B" or mutation_residue_number is None:
            # IAS.py currently mutates chain B residues, which is the ligand chain for this workflow.
            continue

        pdb_path = resolve_pdb_path(pdb_id, base_dir, cleaned_pdb_dir=cleaned_pdb_dir)
        if not pdb_path:
            print(f"Skipping {pdb_id}: no local PDB file found")
            continue

        mutation_residue = mutation_residue_label(mutation_chain, mutation_residue_number)
        cache_key = (pdb_path, mutation_chain)
        if cache_key not in cache:
            try:
                cache[cache_key] = run_convex_hull_neighbor_map(pdb_path, mutation_chain)
            except BaseException as exc:
                print(f"Skipping {pdb_id}: convex-hull failed ({type(exc).__name__}: {exc})")
                cache[cache_key] = None

        cache_value = cache.get(cache_key)
        if cache_value is None:
            continue

        residue_map, convex_stdout, convex_stderr = cache_value
        write_convex_hull_logs(base_dir, pdb_path, mutation_residue, convex_stdout, convex_stderr)

        nearby_protein_residues = residue_map.get(mutation_residue_number, [])

        manifest_rows = make_manifest_rows(pdb_path, row, nearby_protein_residues)
        if not manifest_rows:
            continue

        manifest_name = str(row.get("Mutation.s._cleaned", "")).strip() or f"{pdb_id}_{mutation_chain}{mutation_residue_number}"
        safe_manifest_name = re.sub(r"[^A-Za-z0-9_.-]", "_", manifest_name)
        manifest_path = Path(output_dir) / f"{pdb_id}.{safe_manifest_name}.manifest.csv"

        with open(manifest_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for manifest_row in manifest_rows:
                writer.writerow(manifest_row)

        copy_manifest_to_results(base_dir, pdb_path, mutation_residue, manifest_path)

        written.append((pdb_id, manifest_path, len(manifest_rows)))

    return written


def merge_manifest_files(manifest_paths, output_path):
    manifest_paths = [Path(path) for path in manifest_paths if Path(path).exists()]
    fieldnames = [
        "pdb",
        "mode",
        "mutation_residue",
        "protein_flex_residue",
        "ligand_flex_residue",
        "skempi_pdb",
        "skempi_mutation",
        "wt_aa",
        "mut_aa",
    ]

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
    parser.add_argument("--cleaned-pdb-dir", default=None)
    parser.add_argument("--aggregate-manifest", default=None)
    args = parser.parse_args()

    written = build_all_manifests(
        args.csv_path,
        args.output_dir,
        base_dir=args.base_dir or Path(args.csv_path).parent,
        limit=args.limit,
        cleaned_pdb_dir=args.cleaned_pdb_dir,
    )

    if args.aggregate_manifest:
        manifest_paths = [str(path) for _, path, _ in written]
        merged_path = merge_manifest_files(manifest_paths, args.aggregate_manifest)
        print(f"Wrote merged manifest to {merged_path}")


if __name__ == "__main__":
    main()
