import subprocess
import unittest
from pathlib import Path

from tests.workflow_test_utils import copy_config, gather_expected_outputs
from tests.workflow_test_utils.scenario_testcase import ScenarioTestCase


class TestFailingWorkflows(ScenarioTestCase, unittest.TestCase):
    scenarios_dir = Path(__file__).parent / "failing_scenarios"
    check_strategy = ScenarioTestCase.OutputCheckOptions.NO_CHECK
    link_items = {
        "utils": "link",
        "config": copy_config,
        "workflow": "link"
    }

    def run_scenario(self, scenario_name, expected_state_path):
        expected_outputs = gather_expected_outputs(expected_state_path)
        result = subprocess.run(["snakemake", "--use-conda"] + expected_outputs)
        self.assertNotEqual(0, result.returncode, "Snakemake completed normally")
