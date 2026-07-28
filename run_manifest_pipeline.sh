#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_manifest_pipeline.sh <skempi_csv> <output_dir> [options]

Builds manifests with build_csv_manifests.py, then submits IAS_array.sh through submit_ias_array.sh.
All execution is launched through sbatch.

When run inside a SLURM allocation, this script submits dependent helper jobs and exits
without waiting, so the parent allocation does not time out while children keep running.

Options:
  --cleaned-pdb-dir DIR      Directory with cleaned PDBs (default: SKEMPI2/cleaned_pdbs_for_analysis)
  --aggregate-manifest PATH  Aggregated manifest path (default: <output_dir>/all_experiments.csv)
  --limit N                  Limit number of CSV rows used during manifest build
  --array-limit N            Max concurrent SLURM array tasks (default: 4)
  --build-time TIME          SLURM time limit for manifest build helper (default: 4-00:00:00)
  --submit-time TIME         SLURM time limit for submit helper (default: 00:30:00)
  --ias-time TIME            SLURM time limit for IAS array tasks (default: 4-00:00:00)
  --dry-run                  Print exact commands and expected outputs without executing
  --python-bin PATH          Python interpreter to use
  --partition NAME           SLURM partition for helper jobs (default: compsci)
  --wait-pipeline            Force waiting for helper jobs even when running in SLURM

This script must run inside a SLURM job submitted via sbatch.
EOF
}

quote_cmd() {
  printf '%q ' "$@"
}

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "run_manifest_pipeline.sh must be launched from a SLURM job via sbatch." >&2
  echo "Example: sbatch submit_manifests.sh all" >&2
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
BUILD_TIME="4-00:00:00"
SUBMIT_TIME="00:30:00"
IAS_TIME="4-00:00:00"
DRY_RUN=0
PYTHON_BIN="${PYTHON_BIN:-/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python}"
PARTITION="${PARTITION:-compsci}"
FORCE_WAIT=0

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
    --build-time)
      BUILD_TIME="$2"
      shift 2
      ;;
    --submit-time)
      SUBMIT_TIME="$2"
      shift 2
      ;;
    --ias-time)
      IAS_TIME="$2"
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
    --wait-pipeline)
      FORCE_WAIT=1
      shift
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

submit_cmd=(bash "$SCRIPT_DIR/submit_ias_array.sh" "$AGGREGATE_MANIFEST" "$ARRAY_LIMIT" --time "$IAS_TIME")
if [[ "$DRY_RUN" -eq 1 ]]; then
  submit_cmd+=(--dry-run)
fi

build_wrap="$(quote_cmd "${build_cmd[@]}")"
build_sbatch_cmd=(
  sbatch
  --wait
  --partition="$PARTITION"
  --time="$BUILD_TIME"
  --mem=8G
  --cpus-per-task=1
  --output="$SCRIPT_DIR/logs/manifest_build_%j.out"
  --error="$SCRIPT_DIR/logs/manifest_build_%j.err"
  --wrap "$build_wrap"
)

submit_wrap="cd $(quote_cmd "$SCRIPT_DIR") && $(quote_cmd "${submit_cmd[@]}")"
submit_sbatch_cmd=(
  sbatch
  --partition="$PARTITION"
  --time="$SUBMIT_TIME"
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
wait_for_children=1
if [[ -n "${SLURM_JOB_ID:-}" && "$FORCE_WAIT" -eq 0 ]]; then
  wait_for_children=0
fi

build_sbatch_cmd_effective=("${build_sbatch_cmd[@]}")
submit_sbatch_cmd_effective=("${submit_sbatch_cmd[@]}")

if [[ "$wait_for_children" -eq 1 ]]; then
  build_sbatch_cmd_effective=(sbatch --wait "${build_sbatch_cmd_effective[@]:1}")
  submit_sbatch_cmd_effective=(sbatch --wait "${submit_sbatch_cmd_effective[@]:1}")
fi

build_submit_output="$("${build_sbatch_cmd_effective[@]}")"
echo "$build_submit_output"

build_job_id=""
if [[ "$build_submit_output" =~ Submitted[[:space:]]batch[[:space:]]job[[:space:]]([0-9]+) ]]; then
  build_job_id="${BASH_REMATCH[1]}"
fi
if [[ -z "$build_job_id" && "$build_submit_output" =~ ^([0-9]+) ]]; then
  build_job_id="${BASH_REMATCH[1]}"
fi

if [[ "$wait_for_children" -eq 0 ]]; then
  if [[ -z "$build_job_id" ]]; then
    echo "Could not parse build helper job id from sbatch output." >&2
    exit 2
  fi

  echo "Submitting IAS helper via dependency afterok:$build_job_id ..."
  submit_with_dep=("${submit_sbatch_cmd_effective[@]:0:1}" --dependency="afterok:$build_job_id" "${submit_sbatch_cmd_effective[@]:1}")
  submit_output="$("${submit_with_dep[@]}")"
  echo "$submit_output"

  submit_helper_job_id=""
  if [[ "$submit_output" =~ Submitted[[:space:]]batch[[:space:]]job[[:space:]]([0-9]+) ]]; then
    submit_helper_job_id="${BASH_REMATCH[1]}"
  fi
  if [[ -z "$submit_helper_job_id" && "$submit_output" =~ ^([0-9]+) ]]; then
    submit_helper_job_id="${BASH_REMATCH[1]}"
  fi

  echo "Submission summary"
  echo "  Manifest target: $AGGREGATE_MANIFEST"
  echo "  Build helper job id: $build_job_id"
  if [[ -n "$submit_helper_job_id" ]]; then
    echo "  Submit helper job id: $submit_helper_job_id"
    echo "  Check status: squeue -j $build_job_id,$submit_helper_job_id"
  else
    echo "  Submit helper job id: not detected"
  fi
  echo "  Build logs: $SCRIPT_DIR/logs/manifest_build_<jobid>.out|.err"
  echo "  Submit logs: $SCRIPT_DIR/logs/manifest_submit_<jobid>.out|.err"
  echo "  Results directory: $SCRIPT_DIR/results"
  exit 0
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
submit_output="$("${submit_sbatch_cmd_effective[@]}")"
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
