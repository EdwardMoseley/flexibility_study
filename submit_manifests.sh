#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

usage() {
cat <<'EOF'
Usage:
	sbatch submit_manifests.sh <selector> [array_limit] [--dry-run] [--partition NAME]

Selectors:
	1                   First observation row from SKEMPI CSV
	2                   Second observation row
	1:2                 Inclusive range (rows 1 through 2)
	pdb:5XCO           All rows where the first CSV column (pdb) matches 5XCO
	pdb:5XCO:range:2:3 All rows for PDB 5XCO restricted to observation rows 2 through 3
	all                 All observation rows

Examples:
	sbatch submit_manifests.sh 1
	sbatch submit_manifests.sh 1:2 8
	sbatch submit_manifests.sh pdb:5XCO 9
	sbatch submit_manifests.sh pdb:5XCO:range:2:3 4
	sbatch submit_manifests.sh all 16

Notes:
	- Simple mode uses SKEMPI2/SKEMPI2_processed18Jun26.csv by default.
	- This script must be launched with sbatch.
	- Override defaults with env vars:
			SKEMPI_SOURCE_CSV
			SKEMPI_CLEANED_PDB_DIR
			SKEMPI_RUN_ROOT
EOF
}

is_selector() {
	local value="$1"
	[[ "$value" == "all" || "$value" =~ ^[0-9]+$ || "$value" =~ ^[0-9]+:[0-9]+$ || "$value" =~ ^pdb:[A-Za-z0-9_.-]+$ || "$value" =~ ^pdb:[A-Za-z0-9_.-]+:range:[0-9]+:[0-9]+$ ]]
}

SKEMPI_SOURCE_CSV="${SKEMPI_SOURCE_CSV:-$SCRIPT_DIR/SKEMPI2/SKEMPI2_processed18Jun26.csv}"
SKEMPI_CLEANED_PDB_DIR="${SKEMPI_CLEANED_PDB_DIR:-$SCRIPT_DIR/SKEMPI2/cleaned_pdbs_for_analysis}"
SKEMPI_RUN_ROOT="${SKEMPI_RUN_ROOT:-$SCRIPT_DIR/tmp/observation_batches}"

if [[ $# -lt 1 ]]; then
	usage
	exit 1
fi

SELECTOR="$1"
shift

if [[ "$SELECTOR" == "-h" || "$SELECTOR" == "--help" ]]; then
	usage
	exit 0
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
	echo "submit_manifests.sh must be launched via sbatch." >&2
	echo "Example: sbatch submit_manifests.sh all" >&2
	exit 1
fi

if ! is_selector "$SELECTOR"; then
	echo "Invalid selector: $SELECTOR" >&2
	usage
	exit 1
fi

ARRAY_LIMIT=4
if [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; then
	ARRAY_LIMIT="$1"
	shift
fi

if [[ ! -f "$SKEMPI_SOURCE_CSV" ]]; then
	echo "Source CSV not found: $SKEMPI_SOURCE_CSV" >&2
	exit 1
fi

selector_tag="$(printf '%s' "$SELECTOR" | sed 's/[^A-Za-z0-9_.-]/-/g')"
batch_root="$SKEMPI_RUN_ROOT/$selector_tag"
subset_csv="$batch_root/subset.csv"
output_dir="$batch_root/manifests"

mkdir -p "$batch_root"

if [[ "$SELECTOR" == "all" ]]; then
	cp "$SKEMPI_SOURCE_CSV" "$subset_csv"
elif [[ "$SELECTOR" =~ ^pdb:([A-Za-z0-9_.-]+):range:([0-9]+):([0-9]+)$ ]]; then
	pdb_id="${BASH_REMATCH[1]}"
	start="${BASH_REMATCH[2]}"
	end="${BASH_REMATCH[3]}"
	if (( start < 1 || end < 1 || start > end )); then
		echo "Invalid range selector: $SELECTOR" >&2
		exit 1
	fi
	awk -F',' -v pdb="$pdb_id" -v s="$start" -v e="$end" '
		NR==1 { print; next }
		{
			val = $1
			gsub(/^"/, "", val)
			gsub(/"$/, "", val)
			gsub(/^[ \t]+|[ \t]+$/, "", val)
			if (toupper(val) == toupper(pdb)) {
				row_num++
				if (row_num >= s && row_num <= e) {
					print
				}
			}
		}
	' "$SKEMPI_SOURCE_CSV" > "$subset_csv"
elif [[ "$SELECTOR" =~ ^pdb:([A-Za-z0-9_.-]+)$ ]]; then
	pdb_id="${BASH_REMATCH[1]}"
	awk -F',' -v pdb="$pdb_id" '
		NR==1 { print; next }
		{
			val = $1
			gsub(/^"/, "", val)
			gsub(/"$/, "", val)
			gsub(/^[ \t]+|[ \t]+$/, "", val)
			if (toupper(val) == toupper(pdb)) {
				print
			}
		}
	' "$SKEMPI_SOURCE_CSV" > "$subset_csv"
elif [[ "$SELECTOR" =~ ^([0-9]+):([0-9]+)$ ]]; then
	start="${BASH_REMATCH[1]}"
	end="${BASH_REMATCH[2]}"
	if (( start < 1 || end < 1 || start > end )); then
		echo "Invalid range selector: $SELECTOR" >&2
		exit 1
	fi
	awk -v s="$start" -v e="$end" 'NR==1 || (NR>=s+1 && NR<=e+1)' "$SKEMPI_SOURCE_CSV" > "$subset_csv"
else
	obs="$SELECTOR"
	if (( obs < 1 )); then
		echo "Observation selector must be >= 1" >&2
		exit 1
	fi
	awk -v n="$obs" 'NR==1 || NR==n+1' "$SKEMPI_SOURCE_CSV" > "$subset_csv"
fi

rows=$(( $(wc -l < "$subset_csv") - 1 ))
if (( rows < 1 )); then
	echo "Selector $SELECTOR matched no observation rows in $SKEMPI_SOURCE_CSV" >&2
	exit 1
fi

echo "[submit_manifests.sh] Selected observations: $SELECTOR ($rows rows)"
echo "[submit_manifests.sh] Subset CSV: $subset_csv"
echo "[submit_manifests.sh] Output dir: $output_dir"

pipeline_args=(
	"$subset_csv"
	"$output_dir"
	--cleaned-pdb-dir "$SKEMPI_CLEANED_PDB_DIR"
	--aggregate-manifest "$output_dir/all_experiments.csv"
	--array-limit "$ARRAY_LIMIT"
)

if [[ $# -gt 0 ]]; then
	pipeline_args+=("$@")
fi

exec bash "$SCRIPT_DIR/run_manifest_pipeline.sh" "${pipeline_args[@]}"
