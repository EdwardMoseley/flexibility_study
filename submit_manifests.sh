#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 /path/to/cleaned_pdb_dir [--limit N] [--dry-run]

Builds manifests from the given directory and submits an sbatch array to run IAS.
EOF
}

if [ $# -lt 1 ]; then
  usage
  exit 1
fi

PDB_DIR="$1"
shift || true
LIMIT=0
DRY_RUN=0
while [ $# -gt 0 ]; do
  case "$1" in
    --limit)
      LIMIT="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    *)
      echo "Unknown arg: $1" >&2; exit 1
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

OUT_MANIFEST="$SCRIPT_DIR/manifests.csv"

PYTHON=/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python
CMD=("$PYTHON" build_manifests.py --pdb-dir "$PDB_DIR" --out "$OUT_MANIFEST")
if [ "$LIMIT" -gt 0 ]; then
  CMD+=(--limit "$LIMIT")
fi
if [ "$DRY_RUN" -eq 1 ]; then
  CMD+=(--dry-run)
fi

echo "Building manifests with: ${CMD[*]}"
"${CMD[@]}"

if [ ! -f "$OUT_MANIFEST" ]; then
  echo "Manifest not created: $OUT_MANIFEST" >&2
  exit 2
fi

total_lines=$(wc -l < "$OUT_MANIFEST" | tr -d ' ')
# subtract header if present
if [ "$total_lines" -gt 0 ]; then
  data_rows=$((total_lines - 1))
else
  data_rows=0
fi

if [ "$data_rows" -eq 0 ]; then
  echo "No manifest data rows to submit (manifest has $total_lines line(s)); not submitting"; exit 0
fi

if [ "$DRY_RUN" -eq 1 ]; then
  echo "Dry-run: created $OUT_MANIFEST with $data_rows data row(s); not submitting to SLURM."
  exit 0
fi

echo "Submitting SLURM array with $data_rows tasks"
sbatch --job-name=ias-array --array=0-$((data_rows-1)) --export=MANIFEST="$OUT_MANIFEST" run_ias_task.sh
