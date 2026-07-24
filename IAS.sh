#!/usr/bin/env bash

set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <file_name> <mutation_residue> [protein_flexible_residue] [ligand_flexible_residue]"
    exit 1
fi

FILENAME="${1:-}"
MUTATION_RESIDUE="${2:-}"
PROT_FLEX="${3:-}"
LIG_FLEX_RESIDUE="${4:-}"

if [[ -z "$MUTATION_RESIDUE" ]]; then
    echo "Missing mutation residue (example: B9)" >&2
    exit 1
fi

"$PYTHON_BIN" "$SCRIPT_DIR/IAS.py" "$FILENAME" "$MUTATION_RESIDUE" "$PROT_FLEX" "$LIG_FLEX_RESIDUE"

echo "Completed!"