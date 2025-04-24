import os
import subprocess
import tempfile
from pathlib import Path
from unittest import TestCase
import sys

tests_dir = Path(__file__).parent
project_dir = tests_dir.parent


class RegressionTest(TestCase):
    n_snakemake_cores = 4

    # these are all together for convenience since it's a simple test and a lot of the
    # files created are the same for both
    def test_create_and_assign_hyperpackages(self):
        requested_hyperpackages = [
            "results/hyperpackages/ec_2.7.10.1.refpkg.tar.gz",
            "results/hyperpackages/rhea_10596.refpkg.tar.gz",
        ]
        requested_assignments = [
            "results/assigned_hyperpackages/geneX/rhea_10596.refpkg.tar.gz",
        ]

        # make and switch to a temporary directory
        original_wd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                # link directories
                os.symlink(project_dir / "workflow", "workflow")
                os.symlink(project_dir / "utils", "utils")
                os.symlink(tests_dir / "test_config.yaml", "config.yaml")
                os.mkdir("data")
                os.symlink(tests_dir / "data" / "geneX.fasta", "data/geneX.fasta")

                # make hyperpackages
                hyperpackage_result = subprocess.run(
                    ["snakemake", "--use-conda", "--jobs", str(self.n_snakemake_cores)]
                    + requested_hyperpackages,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    text=True
                )

                # check hyperpackages result
                self.assertEqual(
                    0,
                    hyperpackage_result.returncode,
                    f"Snakemake returned non-zero code for TreeSAPP create: {hyperpackage_result.returncode}\n"
                )
                for file in requested_hyperpackages:
                    file_path = Path(file)
                    self.assertTrue(file_path.exists(), f"{file} was not created")
                    self.assertTrue(file_path.stat().st_size > 0, f"{file} is empty")

                # make assignments
                assign_result = subprocess.run(
                    ["snakemake", "--use-conda", "--jobs", str(self.n_snakemake_cores)]
                    + requested_assignments,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    text=True
                )

                # check assignment results
                self.assertEqual(
                    0,
                    assign_result.returncode,
                    f"Snakemake returned non-zero code for TreeSAPP assign: {assign_result.returncode}\n"
                )
                for file in requested_assignments:
                    file_path = Path(file)
                    self.assertTrue(file_path.exists(), f"{file} was not created")
                    self.assertTrue(file_path.stat().st_size > 0, f"{file} is empty")
            finally:
                os.chdir(original_wd)
