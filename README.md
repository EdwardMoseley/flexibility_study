# flexibility_study

This folder contains a SLURM-first workflow for running IAS/K* calculations from SKEMPI mutation data, using convex hull neighborhood information to define flexibility regimes.

## Overview

Pipeline intent:

1. Read mutation rows from SKEMPI processed CSV.
2. Identify mutated ligand residues from Mutation.s._cleaned (chain letter + residue position).
3. Use convex_hull.py output (via SCOPE) to map each mutated ligand residue to nearby protein residues.
4. Generate IAS experiments for each mutation with three regimes:
	 1. baseline-no-flex
	 2. ligand-flex
	 3. ligand-and-single-protein-flex (one nearby protein residue at a time)
5. Submit SLURM array jobs to run IAS.py for each manifest row.
6. Move resulting TSV files to a single per-PDB/per-ligand results destination.

All compute execution is done through sbatch.

## Key Inputs

Primary SKEMPI input:

- SKEMPI2/SKEMPI2_processed18Jun26.csv

Primary cleaned PDB source:

- SKEMPI2/cleaned_pdbs_for_analysis

Convex hull entrypoint:

- convex_hull.py

## Scripts and Responsibilities

Mutation and manifest generation:

- build_csv_manifests.py
	- Reads a CSV of SKEMPI rows.
	- Resolves PDB paths from cleaned_pdbs_for_analysis.
	- Calls convex_hull.py and parses per-residue neighbor lines.
	- Emits per-mutation manifests and aggregated manifest.

Pilot selection helper:

- select_skempi_pilot.py
	- Selects one PDB with at least N chain-specific mutations.
	- Writes a subset CSV for pilot execution.
	- Supports explicit --pdb selection.

IAS execution:

- IAS.py
	- Runs OSPREY K* for one experiment row.
	- Accepts explicit mutation residue plus optional protein/ligand flexible residues.
	- Writes one TSV per run with a unique tag that includes ligand/protein flexibility state.

Wrappers:

- IAS.sh
	- Thin wrapper for IAS.py.
	- Resolves IAS.py via SLURM_SUBMIT_DIR-aware script location.

- IAS_array.sh
	- Array-task runner for a manifest row.
	- Creates isolated per-task run directory under runs/.
	- Executes IAS.sh.
	- Moves TSV outputs to results/<pdb_basename>/<mutation_residue>/.
	- Deletes seq*.pdb files after run completion to save space.

- submit_ias_array.sh
	- Friendly submit wrapper.
	- Computes row count from manifest and submits IAS_array.sh with array bounds.

## Output Layout

Logs:

- logs/<pdb_basename>_<jobid>_<arraytask>.out
- logs/<pdb_basename>_<jobid>_<arraytask>.err

Per-task scratch/runtime:

- runs/<pdb_basename>_<jobid>_<arraytask>/

Collected TSV results:

- results/<pdb_basename>/<mutation_residue>/*.tsv

Examples:

- results/1EMV.clean/B86/
- results/1EMV.clean/B97/

Notes:

- TSV basenames are unique by mutation residue plus ligand/protein flexibility labels.
- seq*.pdb files are deleted at the end of each task by IAS_array.sh.

## Recommended Invocation (Pilot on One PDB)

Run these from this folder.

1) Create one-PDB pilot subset CSV (auto-select, excluding known entries):

```bash
sbatch --wait \
	--output=logs/pilot_select_%j.out \
	--error=logs/pilot_select_%j.err \
	--partition=compsci --mem=4G --cpus-per-task=1 \
	--wrap="/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python select_skempi_pilot.py \
		SKEMPI2/SKEMPI2_processed18Jun26.csv \
		tmp/pilot_check/subset.csv \
		--chain B --min-mutations 2 --mutations-per-pdb 2 \
		--exclude-pdb 5XCO --exclude-pdb 1C1Y --exclude-pdb 1EMV"
```

2) Or force a specific construct:

```bash
sbatch --wait \
	--output=logs/pilot_select_1emv_%j.out \
	--error=logs/pilot_select_1emv_%j.err \
	--partition=compsci --mem=4G --cpus-per-task=1 \
	--wrap="/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python select_skempi_pilot.py \
		SKEMPI2/SKEMPI2_processed18Jun26.csv \
		tmp/pilot_1emv/subset.csv \
		--pdb 1EMV --chain B --min-mutations 2 --mutations-per-pdb 2"
```

3) Build manifests from the subset CSV:

```bash
sbatch --wait \
	--output=logs/pilot_manifest_%j.out \
	--error=logs/pilot_manifest_%j.err \
	--partition=compsci --mem=8G --cpus-per-task=1 \
	--wrap="/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python build_csv_manifests.py \
		tmp/pilot_check/subset.csv \
		tmp/pilot_check/manifests \
		--cleaned-pdb-dir SKEMPI2/cleaned_pdbs_for_analysis \
		--aggregate-manifest tmp/pilot_check/manifests/all_experiments.csv"
```

4) Submit IAS array using the wrapper script:

```bash
sbatch submit_ias_array.sh tmp/pilot_check/manifests/all_experiments.csv 4
```

Where argument 2 is max concurrent array tasks.

## Manual Convex Hull Example

Single PDB convex hull check:

```bash
sbatch --wait \
	--output=logs/convex_simple_%j.out \
	--error=logs/convex_simple_%j.err \
	--partition=compsci --mem=8G --cpus-per-task=1 \
	--wrap="/home/users/etm33/miniconda3/envs/osprey-jdk17/bin/python convex_hull.py SKEMPI2/cleaned_pdbs_for_analysis/5XCO.clean.pdb"
```

Read output in the generated log under logs/convex_simple_<jobid>.out.

## What Was Changed

Functional updates made in this workflow:

1. Added mutation-aware manifest generation from SKEMPI rows.
2. Added explicit parsing of convex_hull.py residue-neighbor output per design residue.
3. Updated IAS argument flow so mutation residue is explicit and independent from flexibility flags.
4. Added per-task isolated run directories to avoid runtime collisions.
5. Fixed SLURM path resolution issues using SLURM_SUBMIT_DIR-aware script roots.
6. Added clean array submit helper script (submit_ias_array.sh).
7. Added TSV collection into centralized results folders.
8. Added post-run deletion of seq*.pdb files.
9. Updated run-tag naming to avoid duplicate TSV basenames.
10. Added pilot selection script (select_skempi_pilot.py).

## Troubleshooting Notes

1. If convex_hull.py crashes on a specific structure, skip that PDB for pilot and revisit with higher resources.
2. If expected TSVs are not in results/, check:
	 - task logs in logs/
	 - whether run happened before absolute results-path fix
	 - whether files remained in runs/<task>/ and need one-time backfill move
3. If wrapper jobs fail to find scripts, ensure submit scripts are launched from this folder and not copied elsewhere.

## Typical Scale-Up Path

After pilot validation:

1. Increase mutations-per-pdb in select_skempi_pilot.py or feed broader CSV input.
2. Build aggregated manifest.
3. Submit with submit_ias_array.sh using a tuned concurrency limit.
