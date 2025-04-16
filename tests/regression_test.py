import os
import subprocess
import tempfile
from pathlib import Path
from unittest import TestCase

tests_dir = Path(__file__).parent
project_dir = tests_dir.parent


class RegressionTest(TestCase):
    n_snakemake_cores = 4

    def test_create_hyperpackages(self):
        requested_files = [
            "data/hyperpackages/ec_2.7.10.1.refpkg.tar.gz",
            "data/hyperpackages/rhea_10596.refpkg.tar.gz",
            "data/hyperpackages/ec_1.1.3.refpkg.tar.gz",
        ]

        original_wd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)

            os.symlink(project_dir / "workflow", "workflow")
            os.symlink(project_dir / "utils", "utils")
            os.symlink(tests_dir / "test_config.yaml", "config.yaml")

            result = subprocess.run(
                ["snakemake", "--use-conda", "--jobs", str(self.n_snakemake_cores)]
                + requested_files
            )

            self.assertEqual(
                0,
                result.returncode,
                f"Snakemake returned non-zero code: {result.returncode}",
            )

            for file in requested_files:
                file_path = Path(file)
                self.assertTrue(file_path.exists(), f"{file} was not created")
                self.assertTrue(file_path.stat().st_size > 0, f"{file} is empty")

        os.chdir(original_wd)
    
    def test_assign_hyperpackage(self):
        requested_files = [
            "data/assigned_packages/rhea_10596.refpkg.tar.gz",
        ]

        original_wd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)

            os.symlink(project_dir / "workflow", "workflow")
            os.symlink(project_dir / "utils", "utils")
            os.symlink(tests_dir / "test_config.yaml", "config.yaml")
            os.symlink(tests_dir / "geneX.fasta", "geneX.fasta")

            result = subprocess.run(
                ["snakemake", "--use-conda", "--jobs", str(self.n_snakemake_cores)]
                + requested_files
            )

            self.assertEqual(
            
                result.returncode,
                f"Snakemake returned non-zero code: {result.returncode}",
            )

            for file in requested_files:
                file_path = Path(file)
                self.assertTrue(file_path.exists(), f"{file} was not created")
                self.assertTrue(file_path.stat().st_size > 0, f"{file} is empty")

        os.chdir(original_wd)
