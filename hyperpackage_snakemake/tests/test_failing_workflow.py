import subprocess
import unittest
import os

from tests.workflow_runner import WorkflowRunner


class TestFailingWorkflows(WorkflowRunner, unittest.TestCase):
    scenarios_dir = "tests/failure_scenarios"

    @classmethod
    def generate_test(cls, scenario_path, expected_output_path):
        def test(self):
            self.expected_output_path = expected_output_path
            self.setup_workflow_dirs(scenario_path)
            expected_outputs = self.gather_outputs(expected_output_path)
            self.run_workflow(expected_outputs)
            self.check_for_expected_outputs(expected_outputs)

        return test

    def run_workflow(self, expected_outputs):
        result = subprocess.run(["snakemake", "--use-conda"] + expected_outputs)
        self.assertNotEqual(
            0,
            result.returncode,
            f"snakemake returned zero return code",
        )

    def check_for_expected_outputs(self, expected_outputs):
        for expected_output in expected_outputs:
            self.assertFalse(os.path.exists(expected_output), f"{expected_output} was created")
