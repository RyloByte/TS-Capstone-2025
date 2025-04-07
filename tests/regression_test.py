import subprocess
from pathlib import Path
from unittest import TestCase

from unittest_scenarios import IsolatedWorkingDirMixin

tests_dir = Path(__file__).parent
project_dir = tests_dir.parent


class RegressionTest(IsolatedWorkingDirMixin, TestCase):
    external_connections = [
        IsolatedWorkingDirMixin.ExternalConnection(
            external_path=str(project_dir / "utils")
        ),
        IsolatedWorkingDirMixin.ExternalConnection(
            external_path=str(project_dir / "workflow")
        ),
        IsolatedWorkingDirMixin.ExternalConnection(
            external_path=str(tests_dir / "test_config"), internal_path="config"
        ),
    ]

    def test_create_hyperpackages(self):
        requested_files = [
            "data/hyperpackages/ec_2.7.10.1.refpkg.tar.gz",
            "data/hyperpackages/rhea_10596.refpkg.tar.gz",
        ]

        result = subprocess.run(
            ["snakemake", "--use-conda", "--cores", "4"] + requested_files
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
