#!/usr/bin/env python3
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest


class SubmissionScriptsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.script_dir = pathlib.Path(__file__).resolve().parents[1]
        cls.submit_ias_array = cls.script_dir / "submit_ias_array.sh"
        cls.submit_manifests = cls.script_dir / "submit_manifests.sh"
        cls.pipeline_script = cls.script_dir / "run_manifest_pipeline.sh"

    def test_submit_ias_array_dry_run_prints_exact_sbatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = pathlib.Path(tmpdir) / "manifest.csv"
            manifest.write_text(
                "pdb,mode,mutation_residue,protein_flex_residue,ligand_flex_residue,skempi_pdb,skempi_mutation,wt_aa,mut_aa\n"
                "/tmp/a.pdb,baseline-no-flex,B10,,,,,,\n"
                "/tmp/a.pdb,ligand-flex,B10,,B10,,,,\n",
                encoding="utf-8",
            )

            proc = subprocess.run(
                ["bash", str(self.submit_ias_array), str(manifest), "7", "--dry-run"],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Dry-run: would submit IAS array", proc.stdout)
            self.assertIn("--array=1-2%7", proc.stdout)
            self.assertIn("IAS_array.sh", proc.stdout)

    def test_pipeline_dry_run_prints_build_and_submit(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = pathlib.Path(tmpdir) / "subset.csv"
            out_dir = pathlib.Path(tmpdir) / "manifests"
            csv_path.write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "1ABC,IB9A,B,9,1,I,A\n",
                encoding="utf-8",
            )

            proc = subprocess.run(
                [
                    "bash",
                    str(self.pipeline_script),
                    str(csv_path),
                    str(out_dir),
                    "--dry-run",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Build command (sbatch):", proc.stdout)
            self.assertIn("build_csv_manifests.py", proc.stdout)
            self.assertIn("Submit command (sbatch):", proc.stdout)
            self.assertIn("submit_ias_array.sh", proc.stdout)
            self.assertIn("Expected outputs:", proc.stdout)

    def test_submit_manifests_rejects_passthrough_csv_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = pathlib.Path(tmpdir) / "subset.csv"
            out_dir = pathlib.Path(tmpdir) / "manifests"
            csv_path.write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "1ABC,IB9A,B,9,1,I,A\n",
                encoding="utf-8",
            )

            proc = subprocess.run(
                [
                    "bash",
                    str(self.submit_manifests),
                    str(csv_path),
                    str(out_dir),
                    "--dry-run",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(proc.returncode, 0, msg=proc.stdout)
            self.assertIn("Invalid selector", proc.stderr)

    def test_submit_manifests_selector_range_builds_subset_and_dry_runs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = pathlib.Path(tmpdir) / "skempi.csv"
            cleaned_dir = pathlib.Path(tmpdir) / "cleaned"
            run_root = pathlib.Path(tmpdir) / "runs"
            cleaned_dir.mkdir(parents=True, exist_ok=True)
            run_root.mkdir(parents=True, exist_ok=True)

            csv_path.write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "1ABC,IB9A,B,9,1,I,A\n"
                "2DEF,KB12L,B,12,1,K,L\n"
                "3GHI,MB20T,B,20,1,M,T\n",
                encoding="utf-8",
            )

            env = os.environ.copy()
            env["SKEMPI_SOURCE_CSV"] = str(csv_path)
            env["SKEMPI_CLEANED_PDB_DIR"] = str(cleaned_dir)
            env["SKEMPI_RUN_ROOT"] = str(run_root)

            proc = subprocess.run(
                ["bash", str(self.submit_manifests), "1:2", "6", "--dry-run"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Selected observations: 1:2 (2 rows)", proc.stdout)
            self.assertIn("all_experiments.csv\\ 6\\ --dry-run", proc.stdout)

            subset_csv = run_root / "1-2" / "subset.csv"
            self.assertTrue(subset_csv.exists())
            lines = subset_csv.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(3, len(lines))

    def test_submit_manifests_selector_all_works(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = pathlib.Path(tmpdir) / "skempi.csv"
            cleaned_dir = pathlib.Path(tmpdir) / "cleaned"
            run_root = pathlib.Path(tmpdir) / "runs"
            cleaned_dir.mkdir(parents=True, exist_ok=True)
            run_root.mkdir(parents=True, exist_ok=True)

            csv_path.write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "1ABC,IB9A,B,9,1,I,A\n"
                "2DEF,KB12L,B,12,1,K,L\n",
                encoding="utf-8",
            )

            env = os.environ.copy()
            env["SKEMPI_SOURCE_CSV"] = str(csv_path)
            env["SKEMPI_CLEANED_PDB_DIR"] = str(cleaned_dir)
            env["SKEMPI_RUN_ROOT"] = str(run_root)

            proc = subprocess.run(
                ["bash", str(self.submit_manifests), "all", "--dry-run"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Selected observations: all (2 rows)", proc.stdout)

    def test_submit_manifests_selector_pdb_filters_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = pathlib.Path(tmpdir) / "skempi.csv"
            cleaned_dir = pathlib.Path(tmpdir) / "cleaned"
            run_root = pathlib.Path(tmpdir) / "runs"
            cleaned_dir.mkdir(parents=True, exist_ok=True)
            run_root.mkdir(parents=True, exist_ok=True)

            csv_path.write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "5XCO,IB9A,B,9,1,I,A\n"
                "5XCO,DB12A,B,12,1,D,A\n"
                "1ABC,KB12L,B,12,1,K,L\n",
                encoding="utf-8",
            )

            env = os.environ.copy()
            env["SKEMPI_SOURCE_CSV"] = str(csv_path)
            env["SKEMPI_CLEANED_PDB_DIR"] = str(cleaned_dir)
            env["SKEMPI_RUN_ROOT"] = str(run_root)

            proc = subprocess.run(
                ["bash", str(self.submit_manifests), "pdb:5XCO", "5", "--dry-run"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Selected observations: pdb:5XCO (2 rows)", proc.stdout)
            self.assertIn("all_experiments.csv\\ 5\\ --dry-run", proc.stdout)

            subset_csv = run_root / "pdb-5XCO" / "subset.csv"
            self.assertTrue(subset_csv.exists())
            lines = subset_csv.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(3, len(lines))
            self.assertTrue(all(line.startswith("5XCO") or line.startswith("pdb,") for line in lines))

    def test_submit_manifests_selector_pdb_range_filters_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = pathlib.Path(tmpdir) / "skempi.csv"
            cleaned_dir = pathlib.Path(tmpdir) / "cleaned"
            run_root = pathlib.Path(tmpdir) / "runs"
            cleaned_dir.mkdir(parents=True, exist_ok=True)
            run_root.mkdir(parents=True, exist_ok=True)

            csv_path.write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "5XCO,IB9A,B,9,1,I,A\n"
                "5XCO,DB12A,B,12,1,D,A\n"
                "5XCO,LB7A,B,7,1,L,A\n"
                "1ABC,KB12L,B,12,1,K,L\n",
                encoding="utf-8",
            )

            env = os.environ.copy()
            env["SKEMPI_SOURCE_CSV"] = str(csv_path)
            env["SKEMPI_CLEANED_PDB_DIR"] = str(cleaned_dir)
            env["SKEMPI_RUN_ROOT"] = str(run_root)

            proc = subprocess.run(
                ["bash", str(self.submit_manifests), "pdb:5XCO:range:2:3", "4", "--dry-run"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Selected observations: pdb:5XCO:range:2:3 (2 rows)", proc.stdout)
            subset_csv = run_root / "pdb-5XCO-range-2-3" / "subset.csv"
            self.assertTrue(subset_csv.exists())
            lines = subset_csv.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(3, len(lines))
            self.assertIn("5XCO,DB12A", lines[1])
            self.assertIn("5XCO,LB7A", lines[2])

    def test_convex_hull_can_import_archived_find_doublets_helper(self):
        module_path = self.script_dir / "convex_hull.py"
        code = (
            "import importlib.util; "
            f"spec = importlib.util.spec_from_file_location('convex_hull', {str(module_path)!r}); "
            "module = importlib.util.module_from_spec(spec); "
            "spec.loader.exec_module(module); "
            "print('imported')"
        )

        proc = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        self.assertIn("imported", proc.stdout)

    def test_submit_manifests_uses_slurm_submit_dir_for_default_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_root = pathlib.Path(tmpdir)
            script_dir = temp_root / "script-copy"
            submit_dir = temp_root / "submit-dir"
            script_dir.mkdir(parents=True, exist_ok=True)
            submit_dir.mkdir(parents=True, exist_ok=True)

            (script_dir / "submit_manifests.sh").write_text(
                (self.script_dir / "submit_manifests.sh").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            (submit_dir / "run_manifest_pipeline.sh").write_text(
                "#!/usr/bin/env bash\nexit 0\n",
                encoding="utf-8",
            )
            (submit_dir / "SKEMPI2").mkdir(parents=True, exist_ok=True)
            (submit_dir / "SKEMPI2" / "SKEMPI2_processed18Jun26.csv").write_text(
                "pdb,Mutation.s._cleaned,Chain,Residue,is_single_mutant,WT_AA,Mut_AA\n"
                "5M2O,IB9A,B,9,1,I,A\n",
                encoding="utf-8",
            )

            env = os.environ.copy()
            env["SLURM_SUBMIT_DIR"] = str(submit_dir)

            proc = subprocess.run(
                ["bash", str(script_dir / "submit_manifests.sh"), "pdb:5M2O", "--dry-run"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertIn("Selected observations: pdb:5M2O (1 rows)", proc.stdout)


if __name__ == "__main__":
    unittest.main()
