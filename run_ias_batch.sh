#!/usr/bin/env bash
#SBATCH --job-name=ias-batch
#SBATCH --output=logs/ias-batch-%j.out
#SBATCH --error=logs/ias-batch-%j.err
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci

set -euo pipefail

cd /home/users/etm33/src/OSPREY3/flexibility_study

mkdir -p logs manifests_from_csv

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submitting run_ias_batch.sh through sbatch for SLURM compliance"
    sbatch "$0" "$@"
    exit 0
fi

PYTHON_BIN="${PYTHON_BIN:-/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python}"
CSV_PATH="${1:-SKEMPI2/SKEMPI2_processed18Jun26.csv}"
OUTPUT_DIR="${2:-manifests_from_csv}"
LIMIT="${3:-}"
AGGREGATED_MANIFEST="${OUTPUT_DIR}/all_experiments.csv"
CLEANED_PDB_DIR="${CLEANED_PDB_DIR:-SKEMPI2/cleaned_pdbs_for_analysis}"

if [[ -n "$LIMIT" ]]; then
    "$PYTHON_BIN" build_csv_manifests.py "$CSV_PATH" "$OUTPUT_DIR" --cleaned-pdb-dir "$CLEANED_PDB_DIR" --limit "$LIMIT" --aggregate-manifest "$AGGREGATED_MANIFEST"
else
    "$PYTHON_BIN" build_csv_manifests.py "$CSV_PATH" "$OUTPUT_DIR" --cleaned-pdb-dir "$CLEANED_PDB_DIR" --aggregate-manifest "$AGGREGATED_MANIFEST"
fi

if [[ ! -f "$AGGREGATED_MANIFEST" ]]; then
    echo "No aggregated manifest was generated." >&2
    exit 1
fi

experiment_count=$(($(tail -n +2 "$AGGREGATED_MANIFEST" | wc -l)))
if [[ "$experiment_count" -eq 0 ]]; then
    echo "No experiments were generated." >&2
    exit 1
fi

sbatch --array=1-${experiment_count}%10 IAS_array.sh "$AGGREGATED_MANIFEST"
