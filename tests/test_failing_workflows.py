import os
import subprocess
import unittest
from pathlib import Path

from unittest_scenarios import IsolatedWorkingDirMixin, ScenarioTestCaseMixin

from tests.workflow_test_utils import copy_config, gather_files

tests_dir = Path(__file__).parent
project_dir = tests_dir.parent

class FailingWorkflowsTestCase(ScenarioTestCaseMixin, unittest.TestCase):
    scenarios_dir = tests_dir / "failing_scenarios"
    check_strategy = ScenarioTestCaseMixin.OutputChecking.NONE
    external_connections = [
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(project_dir / "utils")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(project_dir / "workflow")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(project_dir / "config"),
                                                   strategy=copy_config)
    ]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        utils_dir = project_dir / "utils"
        if not utils_dir.exists():
            utils_dir.mkdir(exist_ok=True)

    def run_scenario(self, scenario_name: str, expected_state_path: str) -> None:
        items_to_request = gather_files(os.path.join(expected_state_path, "final_state"))
        result = subprocess.run(["snakemake", "--use-conda"] + list(items_to_request))
        self.assertNotEqual(0, result.returncode, "Snakemake completed normally")
