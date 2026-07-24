#!/usr/bin/env bash
#SBATCH --job-name=ias-task
#SBATCH --output=logs/ias_%A_%a.out
#SBATCH --error=logs/ias_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=compsci

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MANIFEST="${MANIFEST:-$SCRIPT_DIR/manifests.csv}"

if [ ! -f "$MANIFEST" ]; then
  echo "Manifest not found: $MANIFEST" >&2
  exit 2
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
  echo "This script is intended to be run as an sbatch --array job (SLURM_ARRAY_TASK_ID must be set)" >&2
  exit 2
fi

IDX=$SLURM_ARRAY_TASK_ID
PYTHON=/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python

MANIFEST="$MANIFEST" IDX="$IDX" "$PYTHON" - <<'PY'
import csv, json, os, subprocess, sys
manifest = os.environ['MANIFEST']
idx = int(os.environ['IDX'])
with open(manifest, newline='') as fh:
  reader = list(csv.DictReader(fh))
  if idx < 0 or idx >= len(reader):
    print(f'No manifest row for index {idx}', file=sys.stderr)
    sys.exit(3)
  row = reader[idx]
  pdb = row['source_pdb']
  prot = row.get('protein_flex_residue','') or ''
  lig = row.get('ligand_flex_residue','') or ''
  print('Running IAS for', pdb, 'prot_flex=', prot, 'lig_flex=', lig)
  ret = subprocess.run([os.path.join('/home/users/etm33/miniconda3/envs/osprey-jdk17/bin','python'), 'IAS.py', pdb, prot, lig], capture_output=True, text=True)
  print(ret.stdout)
  if ret.returncode != 0:
    print(ret.stderr, file=sys.stderr)
    sys.exit(ret.returncode)
PY

exit 0

exit 0
