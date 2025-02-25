import os
import subprocess
import unittest
import tempfile
import shutil
import hashlib
import difflib


class SuccessfulWorkflowTests(unittest.TestCase):
    scenarios_dir = "tests/successful_scenarios"
    outputs_dir = "tests/expected_outputs"
    link_if_missing = ["config"]
    link = ["workflow", "utils"]

    def __new__(cls, *args, **kwargs):
        if not os.path.exists(cls.scenarios_dir):
            raise Exception(f"Cannot find a scenarios dir {cls.scenarios_dir}, no tests will be collected.")

        # create tests from scenarios
        for scenario in os.listdir(cls.scenarios_dir):
            scenario_path = os.path.join(cls.scenarios_dir, scenario)
            if os.path.isdir(scenario_path):
                test_name = f"test_{scenario}"
                setattr(
                    cls,
                    test_name,
                    cls.generate_test(
                        os.path.join(os.getcwd(), scenario_path),
                        os.path.join(os.getcwd(), cls.outputs_dir, scenario),
                    ),
                )
        return super().__new__(cls)

    @classmethod
    def generate_test(cls, scenario_path, expected_output_path):
        def test(self):
            self.expected_output_path = expected_output_path
            self.setup_workflow_dirs(scenario_path)
            expected_outputs = self.gather_outputs(expected_output_path)
            self.run_workflow(expected_outputs)
            self.compare_expected_outputs(expected_outputs)

        return test

    def setUp(self):
        # create a temporary directory to run the workflow in
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_wd = self.temp_dir.name
        # preserve the original working directory
        self.original_wd = os.getcwd()
        # change to the temporary directory for running workflow
        os.chdir(self.temp_wd)

    def setup_workflow_dirs(self, scenario_path):
        # copy scenario into temp directory
        shutil.copytree(scenario_path, self.temp_wd, dirs_exist_ok=True)
        # copy things like config from main directory if not present in scenario
        for item in self.link_if_missing:
            if not os.path.exists(os.path.join(self.temp_wd, item)):
                os.symlink(
                    os.path.join(self.original_wd, item),
                    os.path.join(self.temp_wd, item),
                )
        # link workflow dirs
        for item in self.link:
            os.symlink(
                os.path.join(self.original_wd, item), os.path.join(self.temp_wd, item)
            )

    def gather_outputs(self, expected_output_path):
        expected_outputs = []
        for root, _, files in os.walk(expected_output_path):
            for file in files:
                full_output_path = os.path.join(root, file)
                relative_output_path = os.path.relpath(
                    full_output_path, start=expected_output_path
                )
                expected_outputs.append(relative_output_path)
        return expected_outputs

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
