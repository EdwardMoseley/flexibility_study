#!/usr/bin/env bash
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci

## Call with
## sbatch --array=1-14%1 IAS_array_noproteinflex.sh Boltz_3EIC_MKI_model_0001.pdb

set -euo pipefail

pdb="$1"
protein_flexible_residue="A0"

lig_positions=(B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14)

idx=$((SLURM_ARRAY_TASK_ID - 1))
target_res="${lig_positions[$idx]}"

base=$(basename "$pdb" .pdb)

mkdir -p logs

stdout_log="logs/${base}_${target_res}_flex_${protein_flexible_residue}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
stderr_log="logs/${base}_${target_res}_flex_${protein_flexible_residue}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

exec >"$stdout_log" 2>"$stderr_log"

echo "Running PDB: $pdb"
echo "Ligand mutation position: $target_res"
echo "Protein flexible residue tag: $protein_flexible_residue"
echo "No protein residues will be set flexible."

bash IAS_noproteinflex.sh "$pdb" "$protein_flexible_residue" "$target_res"