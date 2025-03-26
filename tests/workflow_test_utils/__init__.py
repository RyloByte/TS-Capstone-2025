import os
import shutil
from os import PathLike


def copy_config(config_path, test_config_path, ext=".example"):
    os.makedirs(test_config_path, exist_ok=True)
    for root, _, files in os.walk(config_path):
        config_dir = os.path.relpath(root, config_path)
        test_config_dir = os.path.join(test_config_path, config_dir)

        if not os.path.exists(test_config_dir):
            os.makedirs(test_config_dir, exist_ok=True)

        for file in files:
            if file.endswith(ext):
                old_path = os.path.join(root, file)
                new_name = os.path.splitext(file)[0]
                new_path = os.path.join(test_config_dir, new_name)
                shutil.copy(old_path, new_path)


def gather_files(expected_state_path: str | PathLike[str]) -> set[str]:
    expected_outputs = set()
    for root, _, files in os.walk(expected_state_path):
        for file in files:
            full_output_path = os.path.join(root, file)
            relative_output_path = os.path.relpath(full_output_path, start=expected_state_path)
            expected_outputs.add(relative_output_path)
    return expected_outputs