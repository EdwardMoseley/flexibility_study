#!/usr/bin/env bash
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci

## Call with:
# sbatch --array=1-4%1 IAS_array.sh manifests/example_manifest.csv
# sbatch --array=1-11%1 IAS_array.sh manifests/2LOB_Model_One.clean_manifest.csv

set -euo pipefail

manifest="$1"

if [[ -z "${manifest:-}" ]]; then
    echo "Usage: sbatch --array=1-N%1 IAS_array.sh <manifest.csv>" >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "This script must be run as a SLURM array job." >&2
    exit 1
fi

PYTHON=/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python

row_idx=$((SLURM_ARRAY_TASK_ID - 1))
selection="$($PYTHON - <<'PY' "$manifest" "$row_idx"
import csv
import sys

manifest_path = sys.argv[1]
row_idx = int(sys.argv[2])

with open(manifest_path, newline='') as handle:
    rows = list(csv.DictReader(handle))

if row_idx < 0 or row_idx >= len(rows):
    raise SystemExit(f"Row index {row_idx} is out of range for {manifest_path}")

row = rows[row_idx]
print(f"{row['pdb']}|{row.get('protein_flex_residue', '')}|{row.get('ligand_flex_residue', '')}|{row.get('mode', '')}")
PY
)"

IFS='|' read -r pdb protein_residue ligand_residue mode <<< "$selection"

base=$(basename "$pdb" .pdb)

mkdir -p logs

stdout_log="logs/${base}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
stderr_log="logs/${base}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

exec >"$stdout_log" 2>"$stderr_log"

echo "Running PDB: $pdb"
echo "Mode: $mode"
echo "Protein flexible residue: $protein_residue"
echo "Ligand flexible residue: $ligand_residue"

bash IAS.sh "$pdb" "$protein_residue" "$ligand_residue"