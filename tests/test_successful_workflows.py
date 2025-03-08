import subprocess
import unittest
from pathlib import Path

from tests.workflow_test_utils import copy_config, gather_expected_outputs
from tests.workflow_test_utils.scenario_testcase import ScenarioTestCase


class TestSuccessfulWorkflows(ScenarioTestCase, unittest.TestCase):
    scenarios_dir = Path(__file__).parent / "successful_scenarios"
    check_strategy = ScenarioTestCase.OutputCheckOptions.CHECK_FILE_CONTENTS
    link_items = {
        Path(__file__).parent.parent / "utils": "link",
        Path(__file__).parent.parent / "config": copy_config,
        Path(__file__).parent.parent / "workflow": "link"
    }

    def run_scenario(self, scenario_name, expected_state_path):
        expected_outputs = gather_expected_outputs(expected_state_path)
        result = subprocess.run(["snakemake", "--use-conda"] + expected_outputs)
        self.assertEqual(0, result.returncode, f"Snakemake returned non-zero code: {result.returncode}")
