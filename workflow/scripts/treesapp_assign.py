import tarfile
import subprocess
import os
import tempfile
from tqdm import tqdm



def run_treesapp_assign(input_fasta: str, output_dir: str, refpkg_path: str) -> None:
    os.makedirs(output_path, exist_ok=True)
    result = subprocess.run(
        [
            "treesapp", "assign", 
             "-i", input_fasta, 
             "-m", "prot", 
             "--trim_align", 
             "-o", output_dir,
             "--refpkg", refpkg_path
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
    # os.makedirs(assigned_hyperpackage, exist_ok=True)

    # Extract tar.gz into temp dir
    # with tempfile.TemporaryDirectory() as tmp:
    # tmp = tempfile.mkdtemp(prefix="treesapp_debug_")
    tmp = os.path.abspath("treesapp_debug")
    print(f"[DEBUG] Temporary working directory: {tmp}")
    print(f"{tmp} this is the tmp pathhhh")
    with tarfile.open(hyperpackage, "r:gz") as tar:
        members = tar.getmembers()
        if len(members) == 0:
            raise RuntimeError(f"No reference packages found in {hyperpackage}")
        
        tar.extractall(path=tmp)
    
    assigned_dir = os.path.join(tmp, "assigned_packages")
    os.makedirs(assigned_dir, exist_ok=True)

    refpkg_dirs = {member.name.split("/")[0] for member in members if "/" in member.name}

    for i, refpkg in enumerate(tqdm(sorted(refpkg_dirs), desc="Assigning")):
        if i >= 10:
            break
        refpkg_path = os.path.join(tmp, refpkg)
        output_path = os.path.join(assigned_dir, refpkg)
        if os.path.isdir(refpkg_path):
            run_treesapp_assign(fasta, refpkg_path, output_path)

    create_tar_gz(assigned_tar_output, assigned_dir)






    # with tarfile.open(hyperpackage, "r:gz") as tar:
        

    #     with tempfile.TemporaryDirectory() as tmp_dir:
            
    #         tmp_pkg_path = os.path.join(tmp_dir, "assigned_packages")
    #         os.makedirs(tmp_pkg_path, exist_ok=True)            

    #         # Run treesapp assign on each reference package
    #         for i, member in enumerate(tqdm(members, desc="Running treesapp assign", total=len(members))):
    #             if i >= 10:
    #                 break

    #             # We only want to run on top-level .refpkg directories
    #             top_level_dir = member.name.split("/")[0]
    #             refpkg_path = os.path.join(tmp_dir, top_level_dir)

    #             if not os.path.exists(refpkg_path):
    #                 continue
    #             output_path = os.path.join(tmp_pkg_path, top_level_dir)
    #             os.makedirs(output_path, exist_ok=True)
    #             run_treesapp_assign(fasta, output_path, refpkg_path)

    #         create_tar_gz(assigned_hyperpackage, tmp_pkg_path)

    #         print("Done")

            # Bundle into assigned_hyperpackage
            # curr_wd = os.getcwd()
            # assigned_hyperpackage_dir = os.path.join(curr_wd, "data/assigned_hyperpackages")
            # assigned_fasta_dir = os.path.join(assigned_hyperpackage_dir, "data/assigned_hyperpackages")
            # os.makedirs(assigned_fasta_dir, exist_ok=True)


# if __name__ == "__main__":
#     #inputs
#     hyperpackage = snakemake.input[0]
#     fasta = snakemake.input[1]

#     #outputs
#     assigned_hyperpackage = snakemake.output[0]
#     print(os.getcwd())

#     with tarfile.open(hyperpackage, "r:gz") as tar:
#         members = tar.getmembers()
#         if len(members) == 0:
#             raise RuntimeError(f"No reference packages found in {hyperpackage}")

#         original_wd = os.getcwd()
#         tmp_dir = original_wd

#         for member in tqdm(members, desc="Running treesapp assign", total=len(members), unit="reference package"):
#             run_treesapp_assign(fasta, )

