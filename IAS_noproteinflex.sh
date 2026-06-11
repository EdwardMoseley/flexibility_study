#!/usr/bin/env bash

set -euo pipefail

FILENAME="$1"
PROT_FLEX="$2"
TARGET_RES="$3"

python3 IAS_noproteinflex.py "$FILENAME" "$PROT_FLEX" "$TARGET_RES"

echo "Completed!"