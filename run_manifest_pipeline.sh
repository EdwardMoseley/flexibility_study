#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_manifest_pipeline.sh <skempi_csv> <output_dir> [options]

Builds manifests with build_csv_manifests.py, then submits IAS_array.sh through submit_ias_array.sh.
All execution is launched through sbatch.

Options:
  --cleaned-pdb-dir DIR      Directory with cleaned PDBs (default: SKEMPI2/cleaned_pdbs_for_analysis)
  --aggregate-manifest PATH  Aggregated manifest path (default: <output_dir>/all_experiments.csv)
  --limit N                  Limit number of CSV rows used during manifest build
  --array-limit N            Max concurrent SLURM array tasks (default: 4)
  --dry-run                  Print exact commands and expected outputs without executing
  --python-bin PATH          Python interpreter to use
  --partition NAME           SLURM partition for helper jobs (default: compsci)
EOF
}

quote_cmd() {
  printf '%q ' "$@"
}

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

CSV_PATH="$1"
OUTPUT_DIR="$2"
shift 2

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CLEANED_PDB_DIR="$SCRIPT_DIR/SKEMPI2/cleaned_pdbs_for_analysis"
AGGREGATE_MANIFEST=""
LIMIT=""
ARRAY_LIMIT=4
DRY_RUN=0
PYTHON_BIN="${PYTHON_BIN:-/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python}"
PARTITION="${PARTITION:-compsci}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cleaned-pdb-dir)
      CLEANED_PDB_DIR="$2"
      shift 2
      ;;
    --aggregate-manifest)
      AGGREGATE_MANIFEST="$2"
      shift 2
      ;;
    --limit)
      LIMIT="$2"
      shift 2
      ;;
    --array-limit)
      ARRAY_LIMIT="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --python-bin)
      PYTHON_BIN="$2"
      shift 2
      ;;
    --partition)
      PARTITION="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ ! -f "$CSV_PATH" ]]; then
  echo "Input CSV not found: $CSV_PATH" >&2
  exit 1
fi

if [[ -z "$AGGREGATE_MANIFEST" ]]; then
  AGGREGATE_MANIFEST="$OUTPUT_DIR/all_experiments.csv"
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$SCRIPT_DIR/logs"

build_cmd=(
  "$PYTHON_BIN"
  "$SCRIPT_DIR/build_csv_manifests.py"
  "$CSV_PATH"
  "$OUTPUT_DIR"
  --base-dir "$SCRIPT_DIR"
  --cleaned-pdb-dir "$CLEANED_PDB_DIR"
  --aggregate-manifest "$AGGREGATE_MANIFEST"
)

if [[ -n "$LIMIT" ]]; then
  build_cmd+=(--limit "$LIMIT")
fi

submit_cmd=(bash "$SCRIPT_DIR/submit_ias_array.sh" "$AGGREGATE_MANIFEST" "$ARRAY_LIMIT")
if [[ "$DRY_RUN" -eq 1 ]]; then
  submit_cmd+=(--dry-run)
fi

build_wrap="$(quote_cmd "${build_cmd[@]}")"
build_sbatch_cmd=(
  sbatch
  --wait
  --partition="$PARTITION"
  --mem=8G
  --cpus-per-task=1
  --output="$SCRIPT_DIR/logs/manifest_build_%j.out"
  --error="$SCRIPT_DIR/logs/manifest_build_%j.err"
  --wrap "$build_wrap"
)

submit_wrap="cd $(quote_cmd "$SCRIPT_DIR") && $(quote_cmd "${submit_cmd[@]}")"
submit_sbatch_cmd=(
  sbatch
  --wait
  --partition="$PARTITION"
  --mem=2G
  --cpus-per-task=1
  --output="$SCRIPT_DIR/logs/manifest_submit_%j.out"
  --error="$SCRIPT_DIR/logs/manifest_submit_%j.err"
  --wrap "$submit_wrap"
)

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry-run mode enabled"
  printf 'Build command (sbatch): '
  printf '%q ' "${build_sbatch_cmd[@]}"
  echo
  printf 'Submit command (sbatch): '
  printf '%q ' "${submit_sbatch_cmd[@]}"
  echo
  echo "Expected outputs:"
  echo "  - Aggregated manifest: $AGGREGATE_MANIFEST"
  echo "  - Build logs: $SCRIPT_DIR/logs/manifest_build_<jobid>.out|.err"
  echo "  - Submit helper logs: $SCRIPT_DIR/logs/manifest_submit_<jobid>.out|.err"
  echo "  - Task logs: $SCRIPT_DIR/logs/<pdb_basename>_<jobid>_<array_task>.out|.err"
  echo "  - Results: $SCRIPT_DIR/results/<pdb_basename>/<mutation_residue>/*.tsv"
  exit 0
fi

echo "Building manifests via sbatch..."
build_submit_output="$("${build_sbatch_cmd[@]}")"
echo "$build_submit_output"

build_job_id=""
if [[ "$build_submit_output" =~ Submitted[[:space:]]batch[[:space:]]job[[:space:]]([0-9]+) ]]; then
  build_job_id="${BASH_REMATCH[1]}"
fi

if [[ ! -f "$AGGREGATE_MANIFEST" ]]; then
  echo "Aggregated manifest was not created: $AGGREGATE_MANIFEST" >&2
  exit 2
fi

rows=$(( $(wc -l < "$AGGREGATE_MANIFEST") - 1 ))
if [[ "$rows" -lt 1 ]]; then
  echo "Aggregated manifest has no data rows: $AGGREGATE_MANIFEST" >&2
  exit 2
fi

echo "Submitting IAS array..."
submit_output="$("${submit_sbatch_cmd[@]}")"
echo "$submit_output"

submit_helper_job_id=""
if [[ "$submit_output" =~ Submitted[[:space:]]batch[[:space:]]job[[:space:]]([0-9]+) ]]; then
  submit_helper_job_id="${BASH_REMATCH[1]}"
fi

job_id=""
if [[ -n "$submit_helper_job_id" ]]; then
  submit_log="$SCRIPT_DIR/logs/manifest_submit_${submit_helper_job_id}.out"
  if [[ -f "$submit_log" ]]; then
    if grep -Eo 'Submitted batch job [0-9]+' "$submit_log" >/dev/null 2>&1; then
      job_id="$(grep -Eo 'Submitted batch job [0-9]+' "$submit_log" | tail -n1 | awk '{print $4}')"
    fi
  fi
fi

echo "Submission summary"
echo "  Manifest: $AGGREGATE_MANIFEST"
echo "  Rows submitted: $rows"
if [[ -n "$build_job_id" ]]; then
  echo "  Build helper job id: $build_job_id"
fi
if [[ -n "$submit_helper_job_id" ]]; then
  echo "  Submit helper job id: $submit_helper_job_id"
fi
if [[ -n "$job_id" ]]; then
  echo "  IAS array job id: $job_id"
  echo "  Check status: squeue -j $job_id"
  echo "  Logs directory: $SCRIPT_DIR/logs"
else
  echo "  IAS array job id: not detected (check manifest_submit_<jobid>.out)"
  echo "  Logs directory: $SCRIPT_DIR/logs"
fi
echo "  Results directory: $SCRIPT_DIR/results"
