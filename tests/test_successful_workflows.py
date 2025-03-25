import os
import subprocess
import unittest
from pathlib import Path

from tests.workflow_test_utils import copy_config, gather_files

from unittest_scenarios import ScenarioTestCaseMixin, IsolatedWorkingDirMixin


class SuccessfulWorkflowsTestCase(ScenarioTestCaseMixin, unittest.TestCase):
    scenarios_dir = Path(__file__).parent / "successful_scenarios"
    external_connections = [
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "utils")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "workflow")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "config"), strategy=copy_config)
    ]
    match_final_state_exactly = False


    def run_scenario(self, scenario_name: str, scenario_path: str) -> None:
        items_to_request = gather_files(os.path.join(scenario_path, "final_state"))
        result = subprocess.run(["snakemake", "--use-conda"] + list(items_to_request))
        self.assertEqual(0, result.returncode, f"Snakemake returned non-zero code: {result.returncode}")
