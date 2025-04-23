import shutil
import tarfile
import subprocess
import os
import tempfile
from tqdm import tqdm

config = snakemake.config["treesapp_assign"]
MUTE_TREESAPP = config["mute_treesapp"]
EXTRA_ARGS = []
for item in config["extra_args"]:
    EXTRA_ARGS += item.split()

def run_treesapp_assign(input_fasta: str, output_dir: str, refpkg_path: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    result = subprocess.run(
        [
            "treesapp", "assign", 
             "-i", input_fasta,
             "-o", output_dir,
             "--refpkg_dir", refpkg_path,
        ]
        + EXTRA_ARGS,
        capture_output=MUTE_TREESAPP,
        text=True
    )

    if result.returncode != 0:
        if MUTE_TREESAPP:
            print(result.stderr)
        raise RuntimeError(f"Got non-zero return code from TreeSAPP: {result.returncode}")


def create_tar_gz(output_path: str, input_dir: str):
    with tarfile.open(output_path, "w:gz") as tar:
        tar.add(input_dir, arcname=os.path.basename(input_dir))


if __name__ == "__main__":
    # Inputs
    hyperpackage = snakemake.input[0]  # .tar.gz of reference packages
    fasta = snakemake.input[1]

    # Outputs
    assigned_tar_output = snakemake.output[0]

    # Extract tar.gz into temp dir
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(hyperpackage, "r:gz") as tar:
            members = tar.getmembers()
            if len(members) == 0:
                raise RuntimeError(f"No reference packages found in {hyperpackage}")
            
            tar.extractall(path=tmp)
        
        # Create dir in temp dir which holds packages that have completed treesapp assign
        assigned_dir = os.path.join(tmp, "assigned_packages")
        shutil.copy(os.path.join(tmp, "manifest.json"), os.path.join(assigned_dir, "manifest.json"))
        refpkg_dirs = {member.name.split("/")[0] for member in members if "/" in member.name}

        # Perform treesapp assign on each refpkg
        for refpkg in tqdm(sorted(refpkg_dirs), desc="Assigning", unit="refpkg"):
            output_pkg_dir = os.path.join(assigned_dir, refpkg)
            input_pkg_dir = os.path.join(tmp, refpkg, "final_outputs")
            if os.path.isdir(input_pkg_dir):
                run_treesapp_assign(fasta, output_pkg_dir, input_pkg_dir)

        # Create output
        create_tar_gz(assigned_tar_output, assigned_dir)