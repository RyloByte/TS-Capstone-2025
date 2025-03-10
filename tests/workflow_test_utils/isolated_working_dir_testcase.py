import os
import tempfile


class IsolatedWorkingDirTestCase:

    def setUp(self):
        self.test_temp_dir = tempfile.TemporaryDirectory()
        self.test_path = self.test_temp_dir.name
        self.original_wd = os.getcwd()
        os.chdir(self.test_path)

    def tearDown(self):
        os.chdir(self.original_wd)
