#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch submit_ias_array.sh <manifest.csv> [array_limit]" >&2
    exit 1
fi

manifest="$1"
array_limit="${2:-4}"

if [[ ! -f "$manifest" ]]; then
    echo "Manifest not found: $manifest" >&2
    exit 1
fi

rows=$(($(wc -l < "$manifest") - 1))
if [[ "$rows" -lt 1 ]]; then
    echo "Manifest has no data rows: $manifest" >&2
    exit 1
fi

script_dir="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

echo "Submitting IAS array for $manifest with $rows row(s)"
sbatch --array=1-${rows}%${array_limit} "$script_dir/IAS_array.sh" "$manifest"
