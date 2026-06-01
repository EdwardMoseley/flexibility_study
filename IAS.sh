#!/usr/bin/env bash

set -euo pipefail

FILENAME="$1"
PROT_FLEX="$2"

python3 IAS.py "$FILENAME" "$PROT_FLEX"

echo "Completed!"