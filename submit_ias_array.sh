#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch submit_ias_array.sh <manifest.csv> [array_limit] [--dry-run] [--parsable] [--time TIME]" >&2
    exit 1
fi

manifest="$1"
shift

array_limit=4
dry_run=0
parsable=0
time_limit="${IAS_TIME_LIMIT:-4-00:00:00}"

if [[ $# -gt 0 && "${1:-}" =~ ^[0-9]+$ ]]; then
    array_limit="$1"
    shift
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            dry_run=1
            shift
            ;;
        --parsable)
            parsable=1
            shift
            ;;
        --time)
            time_limit="$2"
            shift 2
            ;;
        *)
            echo "Unknown arg: $1" >&2
            exit 1
            ;;
    esac
done

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
array_expr="1-${rows}%${array_limit}"

cmd=(sbatch --time="$time_limit" --array="$array_expr" "$script_dir/IAS_array.sh" "$manifest")

if [[ "$dry_run" -eq 1 ]]; then
    echo "Dry-run: would submit IAS array for $manifest with $rows row(s)"
    echo "Time limit: $time_limit"
    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    echo
    echo "Expected outputs:"
    echo "  - Task logs under: $script_dir/logs"
    echo "  - Results under: $script_dir/results/<pdb_basename>/<mutation_residue>/"
    exit 0
fi

echo "Submitting IAS array for $manifest with $rows row(s)"
submit_output="$(${cmd[@]})"
echo "$submit_output"

if [[ "$parsable" -eq 1 ]]; then
    if [[ "$submit_output" =~ Submitted[[:space:]]batch[[:space:]]job[[:space:]]([0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    fi
fi
