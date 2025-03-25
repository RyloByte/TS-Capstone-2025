import os
import subprocess
import unittest
from pathlib import Path

from unittest_scenarios import IsolatedWorkingDirMixin, ScenarioTestCaseMixin

from tests.workflow_test_utils import copy_config, gather_files


class FailingWorkflowsTestCase(ScenarioTestCaseMixin, unittest.TestCase):
    scenarios_dir = Path(__file__).parent / "failing_scenarios"
    check_strategy = ScenarioTestCaseMixin.OutputChecking.NONE
    external_connections = [
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "utils")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "workflow")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "config"),
                                                   strategy=copy_config)
    ]

    def run_scenario(self, scenario_name: str, expected_state_path: str) -> None:
        items_to_request = gather_files(os.path.join(expected_state_path, "final_state"))
        result = subprocess.run(["snakemake", "--use-conda"] + list(items_to_request))
        self.assertNotEqual(0, result.returncode, "Snakemake completed normally")
