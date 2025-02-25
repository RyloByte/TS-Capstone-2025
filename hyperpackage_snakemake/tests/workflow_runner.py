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
