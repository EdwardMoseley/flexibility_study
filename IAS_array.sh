#!/usr/bin/env bash
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci

## Call with:
# sbatch --array=1-4%1 IAS_array.sh manifests/example_manifest.csv
# sbatch --array=1-11%1 IAS_array.sh manifests/2LOB_Model_One.clean_manifest.csv

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
WORK_ROOT="$SCRIPT_DIR"

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
mutation_residue = row.get('mutation_residue', '') or row.get('ligand_flex_residue', '')
import os
pdb_path = os.path.abspath(row['pdb'])
print(f"{pdb_path}|{mutation_residue}|{row.get('protein_flex_residue', '')}|{row.get('ligand_flex_residue', '')}|{row.get('mode', '')}")
PY
)"

IFS='|' read -r pdb mutation_residue protein_residue ligand_residue mode <<< "$selection"

base=$(basename "$pdb" .pdb)
ligand_folder=$(printf '%s' "$mutation_residue" | tr -cd 'A-Za-z0-9_.-')

mkdir -p "${WORK_ROOT}/logs"

task_dir="${WORK_ROOT}/runs/${base}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$task_dir"
results_dir="${WORK_ROOT}/results/${base}/${ligand_folder}"
mkdir -p "$results_dir"

stdout_log="${WORK_ROOT}/logs/${base}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
stderr_log="${WORK_ROOT}/logs/${base}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

exec >"$stdout_log" 2>"$stderr_log"

echo "Running PDB: $pdb"
echo "Mode: $mode"
echo "Mutation residue: $mutation_residue"
echo "Protein flexible residue: $protein_residue"
echo "Ligand flexible residue: $ligand_residue"
echo "Task working directory: $task_dir"

cd "$task_dir"

bash "$SCRIPT_DIR/IAS.sh" "$pdb" "$mutation_residue" "$protein_residue" "$ligand_residue"

shopt -s nullglob
tsv_files=(*.tsv)
if (( ${#tsv_files[@]} == 0 )); then
    echo "No TSV outputs were produced in $task_dir"
else
    for tsv_file in "${tsv_files[@]}"; do
        mv "$tsv_file" "$results_dir/"
        echo "Moved $tsv_file -> $results_dir/"
    done
fi

seq_files=(seq*.pdb)
if (( ${#seq_files[@]} > 0 )); then
    rm -f -- "${seq_files[@]}"
    echo "Deleted ${#seq_files[@]} seq*.pdb file(s) from $task_dir"
fi

# Keep task logs alongside results for this PDB/mutation to simplify troubleshooting.
cp -f "$stdout_log" "$results_dir/task_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
cp -f "$stderr_log" "$results_dir/task_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"