from unittest import TestCase
import os
import shutil
import subprocess


class TestPassFail(TestCase):
    # list of input setups and output files that should work
    passing_cases: list[tuple[str, str]] = [
        ("good_rhea_id", "data/test-uniprot_mapped.faa")
    ]
    # list of input setups and output files that should not work
    failing_cases: list[tuple[str, str]] = []
    # input prefix - where to find the inputs
    input_prefix = "tests/test_inputs"
    # relevant input / intermediate directories, note: utils is not messed with
    dirs_to_reset: list[str] = ["data", "input", "results"]
    # where existing data/, input/, results/ get put and then restored from
    temp_dir: str = "testing_temp_"

    def setUp(self):
        if not os.getcwd().endswith("hyperpackage_snakemake"):
            raise Exception(f"This test is meant to be run from `hyperpackage_snakemake`. You are in {os.getcwd()}.")

        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
        os.mkdir(self.temp_dir)
        for directory in self.dirs_to_reset:
            if os.path.exists(directory):
                shutil.move(directory, os.path.join(self.temp_dir, directory))

    def tearDown(self):
        for directory in self.dirs_to_reset:
            if os.path.exists(directory):
                shutil.rmtree(directory)
            if os.path.exists(os.path.join(self.temp_dir, directory)):
                shutil.move(os.path.join(self.temp_dir, directory), os.getcwd())
        shutil.rmtree(self.temp_dir)

    def setup_dirs(self, input_dir: str):
        for directory in self.dirs_to_reset:
            if os.path.isdir(directory):
                shutil.rmtree(directory)
        for item in os.listdir(os.path.join(self.input_prefix, input_dir)):
            src = os.path.join(self.input_prefix, input_dir, item)
            dest = item
            if os.path.isdir(src):
                shutil.copytree(src, dest)
            else:
                shutil.copy(src, dest)

    def test_passing(self):
        for input_dir, output_file in self.passing_cases:
            self.setup_dirs(input_dir)
            result = subprocess.run(["snakemake", "--use-conda", output_file])
            self.assertEqual(0, result.returncode, f"Workflow returned status code {result.returncode}")
            self.assertTrue(os.path.exists(output_file), f"Workflow did not produce output file {output_file}")

    def test_failing(self):
        for input_dir, output_file in self.failing_cases:
            self.setup_dirs(input_dir)
            result = subprocess.run(["snakemake", "--use-conda", output_file])
            self.assertNotEqual(result.returncode, 0, "Workflow returned status code 0")
            self.assertFalse(os.path.exists(output_file), f"Workflow produced output file {output_file}")
