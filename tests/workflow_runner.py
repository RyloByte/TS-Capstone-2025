import os
import shutil
import tempfile


class WorkflowRunner:
    scenarios_dir = ""
    outputs_dir = "tests/expected_outputs"
    link_if_missing = ["config", "utils"]
    link = ["workflow"]

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

    def setUp(self):
        # create a temporary directory to run the workflow in
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_wd = self.temp_dir.name
        # preserve the original working directory
        self.original_wd = os.getcwd()
        # change to the temporary directory for running workflow
        os.chdir(self.temp_wd)

    def tearDown(self):
        os.chdir(self.original_wd)

    def copy_example_config(self, real_config_path, temp_config_path, ext=".example"):
        os.makedirs(temp_config_path, exist_ok=True)
        for root, _, files in os.walk(real_config_path):
            relative_path = os.path.relpath(root, real_config_path)
            target_path = os.path.join(temp_config_path, relative_path)

            if not os.path.exists(target_path):
                os.makedirs(target_path)

            for file in files:
                if file.endswith(ext):
                    old_path = os.path.join(root, file)
                    new_name = os.path.splitext(file)[0]
                    new_path = os.path.join(target_path, new_name)
                    os.symlink(old_path, new_path)

    def setup_workflow_dirs(self, scenario_path):
        # copy scenario into temp directory
        shutil.copytree(scenario_path, self.temp_wd, dirs_exist_ok=True)
        # copy things like config from main directory if not present in scenario
        for item in self.link_if_missing:
            temp_item_path = os.path.join(self.temp_wd, item)
            outside_item_path = os.path.join(self.original_wd, item)
            if not os.path.exists(temp_item_path):
                if item == "config":
                    # use example config values
                    self.copy_example_config(outside_item_path, temp_item_path)
                else:
                    # just link the directory
                    # need to create it if there isn't one at the top level
                    if not os.path.exists(outside_item_path):
                        print(f"{outside_item_path} not found, creating...")
                        os.makedirs(outside_item_path, exist_ok=True)
                    os.symlink(
                        outside_item_path,
                        temp_item_path,
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
