#!/usr/bin/env bash
#SBATCH --job-name=build-manifests
#SBATCH --output=build_manifests_%j.out
#SBATCH --error=build_manifests_%j.err
#SBATCH --partition=compsci
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/home/users/etm33/src/OSPREY3/flexibility_study

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch submit_manifests_batch.sh /path/to/cleaned_pdb_dir [--limit N] [--dry-run]" >&2
  exit 1
fi

PDB_DIR="$1"
shift || true

# ensure conda environment available on compute node
if [ -f "$HOME/.bashrc" ]; then
  # shellcheck disable=SC1090
  source "$HOME/.bashrc"
fi
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate osprey-jdk17 2>/dev/null || true
fi

echo "Building manifests for $PDB_DIR (inside SLURM job)"
./submit_manifests.sh "$PDB_DIR" "$@"

echo "Batch builder finished"
