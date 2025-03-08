import os
import subprocess
import unittest
import hashlib
import itertools

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
            # check the file was produced
            self.assertTrue(
                os.path.exists(expected_output),
                f"Workflow failed to produce: {expected_output}",
            )
            if self.is_text_file(correct_output_path):
                print(f"Comparing {correct_output_path} by text")
                self.assertFilesTextEqual(correct_output_path, expected_output)
            else:
                print(f"Comparing {correct_output_path} by hash")
                self.assertFilesHashEqual(correct_output_path, expected_output)

    def assertFilesHashEqual(self, correct_path, actual_path):
        self.assertTrue(self.get_hash(correct_path) == self.get_hash(actual_path), f"Hash of {actual_path} does not match {correct_path}")

    def assertFilesTextEqual(self, correct_output_path, actual_output_path):
        with open(correct_output_path, "r", newline=None) as f_correct, open(actual_output_path, "r", newline=None) as f_actual:
            for i, (line_correct, line_actual) in enumerate(itertools.zip_longest(f_correct, f_actual, fillvalue=None)):
                self.assertIsNotNone(line_correct, f"{correct_output_path} has {i} lines, {actual_output_path} is longer")
                self.assertIsNotNone(line_actual, f"{actual_output_path} ends on line {i}, {correct_output_path} continues")
                self.assertEqual(line_correct, line_actual, f"{actual_output_path} does not match {correct_output_path} on line {i}: '{line_correct}' vs. '{line_actual}'")

    def get_hash(self, file_path):
        with open(file_path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()

    def is_text_file(self, file_path):
        try:
            with open(file_path, "r") as f:
                _ = f.read(1)
        except UnicodeDecodeError:
            return False
        return True
