import enum
import hashlib
import itertools
import os
import shutil
import tarfile
import tempfile
import zipfile

from tests.workflow_test_utils.isolated_working_dir_testcase import IsolatedWorkingDirTestCase


class ScenarioTestCase(IsolatedWorkingDirTestCase):
    """
    Runs scenarios from the given directory. Each scenario should contain directories `initial_state` and `final_state`.

    When checking outputs, extra files produced by the scenario that are not present in the `final_state` will not be
    checked or penalized.
    """
    class OutputCheckOptions(enum.Enum):
        NO_CHECK = 0
        CHECK_FILE_NAMES = 1
        CHECK_FILE_CONTENTS = 2

    scenarios_dir: str = None
    check_strategy: OutputCheckOptions = OutputCheckOptions.CHECK_FILE_CONTENTS
    link_items: dict = {}

    def run_scenario(self, scenario_name, expected_state_path):
        """
        Run a scenario. A scenario run should depend on up to 3 things:
        - scenario name, provided by arg
        - expected output, provided by absolute path arg
        - initial state, provided in working directory at runtime
        """
        raise NotImplementedError("Please provide a function for running a scenario")

    def __new__(cls, *args, **kwargs):
        # check for valid scenarios directory
        if cls.scenarios_dir is None:
            raise Exception("Please provide a scenarios_dir")
        if not os.path.exists(cls.scenarios_dir):
            raise Exception(f"Could not find scenarios_dir {cls.scenarios_dir}")

        for scenario in os.listdir(cls.scenarios_dir):
            scenario_path = os.path.join(cls.scenarios_dir, scenario)
            if os.path.isdir(scenario_path):
                test_name = f"test_{scenario}"
                initial_state_path = os.path.abspath(os.path.join(scenario_path, "initial_state"))
                expected_state_path = os.path.abspath(os.path.join(scenario_path, "final_state"))
                setattr(cls, test_name, cls.generate_test(scenario, initial_state_path, expected_state_path))
        return super().__new__(cls)

    @classmethod
    def generate_test(cls, scenario_name, initial_state_path, expected_state_path):
        def test_func(self):
            self.copy_initial_state(initial_state_path)
            self.perform_linking()
            self.run_scenario(scenario_name, expected_state_path)
            self.check_final_state(expected_state_path)
        return test_func

    def copy_initial_state(self, scenario_path):
        shutil.copytree(scenario_path, self.test_path, dirs_exist_ok=True)

    def perform_linking(self):
        for item, strategy in self.link_items.items():
            external_item_path = item if os.path.isabs(item) else os.path.join(self.original_wd, item)
            internal_item_path = os.path.join(self.test_path, os.path.basename(item))
            if not os.path.exists(internal_item_path):
                if strategy == "link-make" or strategy == "copy-make":
                    if not os.path.exists(external_item_path):
                        os.makedirs(external_item_path, exist_ok=True)
                if strategy == "copy" or strategy == "copy-make":
                    if os.path.isdir(external_item_path):
                        shutil.copytree(external_item_path, internal_item_path, dirs_exist_ok=True)
                    else:
                        shutil.copy(external_item_path, internal_item_path)
                elif strategy == "link" or strategy == "link-make":
                    os.symlink(external_item_path, internal_item_path, target_is_directory=os.path.isdir(external_item_path))
                elif callable(strategy):
                    strategy(external_item_path, internal_item_path)
                else:
                    raise Exception(f"Unknown linking strategy {strategy} for {item}")

    def check_final_state(self, expected_state_path):
        if self.check_strategy == self.OutputCheckOptions.NO_CHECK:
            return

        if self.check_strategy == self.OutputCheckOptions.CHECK_FILE_NAMES:
            for root, _, files in os.walk(expected_state_path):
                for file in files:
                    expected_item_path = os.path.join(root, file)
                    self.assertTrue(os.path.exists(expected_item_path), f"Did not find expected item {expected_item_path} in scenario end state")

        if self.check_strategy == self.OutputCheckOptions.CHECK_FILE_CONTENTS:
            self.compare_paths(expected_state_path, os.getcwd())

    def compare_paths(self, expected_path, actual_path):
        self.assertTrue(os.path.exists(actual_path), f"{actual_path} does not exist")

        # directory
        if os.path.isdir(expected_path):
            self.compare_directories(expected_path, actual_path)

        # archive
        elif self.is_archive_file(expected_path):
            self.compare_archives(expected_path, actual_path)

        # text file
        elif self.is_text_file(expected_path):
            self.compare_text_files(expected_path, actual_path)

        # binary file
        else:
            self.compare_file_hashes(expected_path, actual_path)

    def compare_directories(self, expected_directory, actual_directory):
        self.assertTrue(os.path.isdir(actual_directory), f"{actual_directory} is not a directory")
        for item in os.listdir(expected_directory):
            self.compare_paths(os.path.join(expected_directory, item), os.path.join(actual_directory, item))

    def is_archive_file(self, file_path):
        archive_extensions = [".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tar.xz", ".txz"]
        return any(file_path.endswith(extension) for extension in archive_extensions)

    def compare_archives(self, expected_archive, actual_archive):
        expected_extracted = tempfile.TemporaryDirectory()
        actual_extracted = tempfile.TemporaryDirectory()

        if expected_archive.endswith(".zip"):
            with zipfile.ZipFile(expected_archive, "r") as f:
                f.extractall(expected_extracted.name)
            with zipfile.ZipFile(actual_archive, "r") as f:
                f.extractall(actual_extracted.name)

        elif expected_archive.endswith(".tar"):
            with tarfile.open(expected_archive) as f:
                f.extractall(expected_extracted.name)
            with tarfile.open(actual_archive) as f:
                f.extractall(actual_extracted.name)

        elif expected_archive.endswith(".tar.gz") or expected_archive.endswith(".tgz"):
            with tarfile.open(expected_archive, "r:gz") as f:
                f.extractall(expected_extracted.name)
            with tarfile.open(actual_archive, "r:gz") as f:
                f.extractall(actual_extracted.name)

        elif expected_archive.endswith(".tar.bz2") or expected_archive.endswith(".tbz2"):
            with tarfile.open(expected_archive, "r:bz2") as f:
                f.extractall(expected_extracted.name)
            with tarfile.open(actual_archive, "r:bz2") as f:
                f.extractall(actual_extracted.name)

        elif expected_archive.endswith(".tar.xz") or expected_archive.endswith(".txz"):
            with tarfile.open(expected_archive, "r:xz") as f:
                f.extractall(expected_extracted.name)
            with tarfile.open(actual_archive, "r:xz") as f:
                f.extractall(actual_extracted.name)

        self.compare_paths(expected_extracted.name, actual_extracted.name)

    def is_text_file(self, file_path):
        try:
            with open(file_path, "r") as f:
                _ = f.read(1)
        except UnicodeDecodeError:
            return False
        return True

    def compare_text_files(self, expected_text_file, actual_text_file):
        with open(expected_text_file, "r", newline=None) as f_expected, open(actual_text_file, "r", newline=None) as f_actual:
            for i, (line_expected, line_actual) in enumerate(itertools.zip_longest(f_expected, f_actual, fillvalue=None)):
                self.assertIsNotNone(line_actual, f"{actual_text_file} ends on line {i + 1}, expected to continue")
                self.assertIsNotNone(line_expected, f"{actual_text_file} continues past line {i}, expected to end")
                self.assertEqual(line_expected, line_actual, f"{actual_text_file} does not match {expected_text_file} on line {i + 1}")

    def compare_file_hashes(self, expected_file, actual_file):
        with open(expected_file, "rb") as f:
            expected_hash = hashlib.sha256(f.read()).hexdigest()
        with open(actual_file, "rb") as f:
            actual_hash = hashlib.sha256(f.read()).hexdigest()
        self.assertEqual(expected_hash, actual_hash, f"Hash of {actual_file} does not match {expected_file}")
