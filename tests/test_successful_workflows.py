import os
import subprocess
import unittest
from pathlib import Path

from tests.workflow_test_utils import copy_config, gather_files
from tests.workflow_test_utils.scenario_testcase import ScenarioTestCase

from unittest_scenarios import ScenarioTestCaseMixin, IsolatedWorkingDirMixin


class SuccessfulWorkflowsTestCase(ScenarioTestCaseMixin, unittest.TestCase):
    scenarios_dir = Path(__file__).parent / "successful_scenarios"
    external_connections = [
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "utils")),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "workflow"), strategy="copy"),
        IsolatedWorkingDirMixin.ExternalConnection(external_path=str(Path(__file__).parent.parent / "config"), strategy=copy_config)
    ]
    extra_final_items_allowed = True


    def run_scenario(self, scenario_name: str, scenario_path: str) -> None:
        # if os.path.exists(os.path.join(scenario_path, "initial_state")):
        #     input_items = gather_files(os.path.join(scenario_path, "initial_state"))
        # else:
        #     input_items = set()
        # output_items = gather_files(os.path.join(scenario_path, "final_state"))
        # items_to_request = output_items - input_items

        items_to_request = gather_files(os.path.join(scenario_path, "final_state"))
        result = subprocess.run(["snakemake", "--use-conda"] + list(items_to_request))
        self.assertEqual(0, result.returncode, f"Snakemake returned non-zero code: {result.returncode}")


# class TestSuccessfulWorkflows(ScenarioTestCase, unittest.TestCase):
#     scenarios_dir = Path(__file__).parent / "successful_scenarios"
#     check_strategy = ScenarioTestCase.OutputCheckOptions.CHECK_FILE_CONTENTS
#     link_items = {
#         Path(__file__).parent.parent / "utils": "link-make",
#         Path(__file__).parent.parent / "config": copy_config,
#         Path(__file__).parent.parent / "workflow": "link"
#     }
#
#     def run_scenario(self, scenario_name, expected_state_path):
#         expected_outputs = gather_files(expected_state_path)
#         result = subprocess.run(["snakemake", "--use-conda"] + expected_outputs)
#         self.assertEqual(0, result.returncode, f"Snakemake returned non-zero code: {result.returncode}")
