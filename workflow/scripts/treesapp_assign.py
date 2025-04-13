import tarfile
import subprocess
import os
import tempfile
from tqdm import tqdm

def run_treesapp_assign(input_fasta: str, output_dir: str, refpkg_path: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    result = subprocess.run(
        [
            "treesapp", "assign", 
             "-i", input_fasta, 
             "-m", "prot", 
             "--trim_align", 
             "-o", output_dir,
             "--refpkg_dir", refpkg_path
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"[ERROR] Failed on {refpkg_path}")
        print(result.stderr)

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
        
        assigned_dir = os.path.join(tmp, "assigned_packages")
        # # os.makedirs(assigned_dir, exist_ok=True)

        refpkg_dirs = {member.name.split("/")[0] for member in members if "/" in member.name}

        for i, refpkg in enumerate(tqdm(sorted(refpkg_dirs), desc="Assigning")):
            # if i >= 3:
            #     break
            output_pkg_dir = os.path.join(assigned_dir, refpkg)
            input_pkg_dir = os.path.join(tmp, refpkg, "final_outputs")
            if os.path.isdir(input_pkg_dir):
                run_treesapp_assign(fasta, output_pkg_dir, input_pkg_dir)

        create_tar_gz(assigned_tar_output, assigned_dir)