import os
import subprocess
import unittest
import hashlib
import difflib

from tests.workflow_runner import WorkflowRunner


class TestSuccessfulWorkflows(WorkflowRunner, unittest.TestCase):
    scenarios_dir = "tests/successful_scenarios"

    @classmethod
    def generate_test(cls, scenario_path, expected_output_path):
        def test(self):
            self.expected_output_path = expected_output_path
            self.setup_workflow_dirs(scenario_path)
            expected_outputs = self.gather_outputs(expected_output_path)
            self.run_workflow(expected_outputs)
            self.compare_expected_outputs(expected_outputs)

        return test

    def run_workflow(self, expected_outputs):
        result = subprocess.run(["snakemake", "--use-conda"] + expected_outputs)
        self.assertEqual(
            0,
            result.returncode,
            f"snakemake returned non-zero return code: {result.returncode}",
        )

    def compare_expected_outputs(self, expected_outputs):
        for expected_output in expected_outputs:
            correct_output_path = os.path.join(
                self.expected_output_path, expected_output
            )
            # raise Exception(correct_output_path)
            # check the file was produced
            self.assertTrue(
                os.path.exists(expected_output),
                f"Workflow failed to produce: {expected_output}",
            )
            # check the hash matches
            if self.get_hash(expected_output) != self.get_hash(correct_output_path):
                diff = self.get_diff(correct_output_path, expected_output)
                self.fail(
                    f"Output file {expected_output} does not match:\n{''.join(diff)}"
                )

    def get_hash(self, file_path):
        with open(file_path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()

    def get_diff(self, file_path1, file_path2):
        with open(file_path1, "r") as f1, open(file_path2, "r") as f2:
            return list(
                difflib.unified_diff(
                    f1.readlines(),
                    f2.readlines(),
                    fromfile=file_path1,
                    tofile=file_path2,
                )
            )
