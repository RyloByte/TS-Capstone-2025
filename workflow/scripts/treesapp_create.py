from snakemake.script import snakemake
import subprocess
import tarfile
import tempfile
import os
from pathlib import Path


config = snakemake.config["treesapp_create"]
extra_args = config["extra_args"]


def run_treesapp(input_faa: str, refpkg_name: str, output_dir: str):
    print(f"Running TreeSAPP, input: {input_faa} output: {output_dir}")

    result = subprocess.run(
        [
            "treesapp",
            "create",
            "-i",
            input_faa,
            "-c",
            refpkg_name,
            "-o",
            output_dir,
            "--headless",
        ],
        capture_output=True,
        text=True,
    )

    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())

    if result.returncode != 0:
        raise RuntimeError(f"Got non-zero return code from TreeSAPP: {result.returncode}")


if __name__ == "__main__":
    # extract the sequences from input
    with tempfile.TemporaryDirectory() as extracted_input:
        with tarfile.open(snakemake.input[1], "r:gz") as tar:
            tar.extractall(path=extracted_input)

        # create a temporary directory for the output
        with tempfile.TemporaryDirectory() as output_dir:

            # iterate through each faa file
            for cluster_file in os.listdir(extracted_input):
                # expecting clustering_file = {sample}-{exemplar uniprot ID}-sequence_cluster.faa
                exemplar_id = cluster_file.split("-")[-2]

    # # create temp output
    # with tempfile.TemporaryDirectory() as temp_output:
    #     # open input tarfile
    #     with tarfile.open(snakemake.input[0], "r:gz") as tar:
    #         # go through each .faa file in tar
    #         for member in tar.getmembers():
    #             # write that tar to a tempfile
    #             with tempfile.NamedTemporaryFile() as temp_fasta:
    #                 temp_fasta.write(tar.extractfile(member).read())
    #                 temp_fasta.flush()
    #
    #                 # create output dir
    #                 refpkg_output = os.path.join(temp_output, Path(member).stem)
    #                 os.makedirs(refpkg_output, exist_ok=True)
    #
    #                 # run treesapp
    #                 refpkg_name = generate_refpkg_name(6)
    #                 run_treesapp(temp_fasta.name, refpkg_name, refpkg_output)
    #
    #     # gather the outputs into a tarfile
    #     with tarfile.open(snakemake.output[0], "w:gz") as tar:
    #         tar.add(temp_output, arcname=".")
