import itertools
import json
import string
import sys

from snakemake.script import snakemake
import subprocess
import tarfile
import tempfile
import os
from pathlib import Path


config = snakemake.config["treesapp_create"]
extra_args: str = config["extra_args"]
ref_pkg_name_length: int = config["ref_pkg_name_length"]
show_treesapp_output: bool = config["show_treesapp_output"]


def ref_pkg_name_generator(n: int):
    characters = string.ascii_uppercase + string.digits
    for combo in itertools.product(characters, repeat=n):
        yield "".join(combo)


def run_treesapp(input_faa: str, ref_pkg_name: str, output_dir: str, taxonomy_file: str) -> None:
    print(f"Running TreeSAPP, input: {input_faa} output: {output_dir}")

    result = subprocess.run(
        [
            "treesapp",
            "create",
            "-i",
            input_faa,
            "-c",
            ref_pkg_name,
            "-o",
            output_dir,
            "--headless",
            "--seqs2lin",
            taxonomy_file,
            extra_args
        ],
        stdout=sys.stdout if show_treesapp_output else None,
        stderr=sys.stderr if show_treesapp_output else None,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Got non-zero return code from TreeSAPP: {result.returncode}")


if __name__ == "__main__":
    swissprot_data_file = Path(snakemake.input[0])
    swissprot_taxonomy_file = Path(snakemake.input[1])
    clusters_archive = Path(snakemake.input[2])
    hyperpackage_output = Path(snakemake.output[0])

    assert swissprot_data_file.suffixes == [".tsv", ".gz"], f"Expected {swissprot_data_file} to be a .tsv.gz swissprot data file"
    assert swissprot_taxonomy_file.suffixes == [".tsv", ".gz"], f"Expected {swissprot_taxonomy_file} to be a .tsv.gz swissprot taxonomy file"
    assert clusters_archive.suffixes == [".tar", ".gz"], f"Expected {clusters_archive} to be a .tar.gz clusters archive"
    assert hyperpackage_output.suffixes == [".tar", ".gz"], f"Expected {hyperpackage_output} to be a .tar.gz hyperpackage output"


    name_generator = ref_pkg_name_generator(ref_pkg_name_length)
    hyper_package_manifest = {"ref_pkgs": []}
    expected_ending = "-hyperpackage.tar.gz"
    if not hyperpackage_output.name.endswith(expected_ending):
        print("Got an unexpected hyperpackage output name, can't determine activity number used to create it")
    else:
        hyper_package_manifest["activity_number"] = hyperpackage_output.name[:-len(expected_ending)]

    # create a temporary directory for the output
    with tempfile.TemporaryDirectory() as output_dir:

        # extract the sequences from input
        with tempfile.TemporaryDirectory() as extracted_input:
            with tarfile.open(clusters_archive, "r:gz") as tar:
                tar.extractall(path=extracted_input)

            # iterate through each faa file
            for cluster_file in os.listdir(extracted_input):
                cluster_file = Path(cluster_file)
                if cluster_file.suffix != ".faa":
                    print(f"Skipping non-faa file: {cluster_file}")
                    continue
                cluster_name = cluster_file.stem
                pkg_name = next(name_generator)
                hyper_package_manifest["ref_pkgs"].append({"pkg_name": pkg_name, "exemplar_sequence": cluster_name})

                run_treesapp(input_faa=cluster_file.name, ref_pkg_name=pkg_name, output_dir=output_dir, taxonomy_file=swissprot_taxonomy_file.name)

        # write manifest file
        with open(os.path.join(output_dir, "hyper_package_manifest.json"), "w") as manifest_f:
            json.dump(hyper_package_manifest, manifest_f)

        # put everything in the output archive
        with tarfile.open(hyperpackage_output, "w:gz") as output_tar:
            for file in os.listdir(output_dir):
                tar.add(os.path.join(os.path.join(output_dir, file)), arcname=file)
