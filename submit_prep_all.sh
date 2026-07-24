#!/bin/bash
#SBATCH --job-name=prep-all-skempi
#SBATCH --output=prep_all_%j.out
#SBATCH --error=prep_all_%j.err
#SBATCH --partition=compsci
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/home/users/etm33/src/OSPREY3/flexibility_study

set -euo pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

cd "$SCRIPT_DIR"

INPUT_DIR="$SCRIPT_DIR/SKEMPI2/pdbs_for_analysis"
OUTPUT_DIR="$SCRIPT_DIR/SKEMPI2/cleaned_pdbs_for_analysis"
LOG_DIR="$SCRIPT_DIR/slurm_logs"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Optional conda activation if this environment is available on the cluster.
if [ -f "$HOME/.bashrc" ]; then
  # shellcheck disable=SC1090
  source "$HOME/.bashrc"
fi

if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate osprey-jdk17 2>/dev/null || true
fi

shopt -s nullglob
pdb_files=("$INPUT_DIR"/*.pdb)

if [ ${#pdb_files[@]} -eq 0 ]; then
  echo "No .pdb files found in $INPUT_DIR"
  exit 0
fi

for input_path in "${pdb_files[@]}"; do
  input_name="$(basename "$input_path")"
  stem="${input_name%.pdb}"
  cleaned_name="${stem}.clean.pdb"
  cleaned_path="$OUTPUT_DIR/$cleaned_name"

  echo "Submitting $input_name"
  sbatch \
    --job-name="prep-${stem}" \
    --partition=compsci \
    --mem=128G \
    --cpus-per-task=4 \
    --output="$LOG_DIR/${stem}_%j.out" \
    --error="$LOG_DIR/${stem}_%j.err" \
    --wrap="cd '$SCRIPT_DIR' && ./prep.sh '$input_path' && mkdir -p '$OUTPUT_DIR' && mv '$INPUT_DIR/$cleaned_name' '$cleaned_path'"
done
