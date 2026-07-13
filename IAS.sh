#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <file_name> [protein_flexible_residue] [ligand_residue]"
    exit 1
fi

FILENAME="${1:-}"
PROT_FLEX="${2:-}"
LIG_RESIDUE="${3:-}"

python3 IAS.py "$FILENAME" "$PROT_FLEX" "$LIG_RESIDUE"

echo "Completed!"