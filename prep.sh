#!/usr/bin/env bash
#SBATCH --output=submit.out
#SBATCH --error=submit.err
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4 #48
#SBATCH --partition=compsci

set -euo pipefail

FILENAME="${1:-}"
if [[ -z "$FILENAME" ]]; then
	echo "Usage: sbatch prep.sh /path/to/input.pdb" >&2
	exit 1
fi

PYTHON=/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python
"$PYTHON" prep.py "$FILENAME"

echo "Completed!"
