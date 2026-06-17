#!/usr/bin/env bash
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci

## Call with:
# sbatch --array=1-4%1 IAS_array.sh Boltz_3EIC_MKI_model_0001.pdb A203 A204 A207 A217
# B4
# sbatch --array=1-4%1 IAS_array.sh 2LOB_Model_One.clean.pdb A50 A89 A90 A91 A92 A93 A95

set -euo pipefail

pdb="$1"
shift

prot_flex_residues=("$@")

idx=$((SLURM_ARRAY_TASK_ID - 1))
residue="${prot_flex_residues[$idx]}"

base=$(basename "$pdb" .pdb)

mkdir -p logs

stdout_log="logs/${base}_${residue}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
stderr_log="logs/${base}_${residue}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

exec >"$stdout_log" 2>"$stderr_log"

echo "Running PDB: $pdb"
echo "Protein flexible residue: $residue"

bash IAS.sh "$pdb" "$residue"