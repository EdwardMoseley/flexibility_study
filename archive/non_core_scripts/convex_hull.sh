#!/usr/bin/env bash

set -euo pipefail

FILENAME="${1:-}"

if [[ -z "$FILENAME" ]]; then
    echo "Usage: sbatch convex_hull.sh /path/to/input.pdb" >&2
    exit 1
fi

PYTHON=/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python

which "$PYTHON"
"$PYTHON" -V
export PYTHONPATH=${PYTHONPATH:-}
echo "PYTHONPATH=${PYTHONPATH}"
"$PYTHON" -c "import typing_extensions; print('typing_extensions', typing_extensions.__file__, hasattr(typing_extensions,'TypeIs'))"
"$PYTHON" -c "import sys; print('sys.path', sys.path)"
"$PYTHON" -c "import Bio; from Bio.PDB import PDBParser; print('Bio OK', Bio.__file__)"

"$PYTHON" convex_hull.py "$FILENAME"

echo "Completed!"