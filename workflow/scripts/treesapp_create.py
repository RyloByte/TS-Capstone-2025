import itertools
import json
import os
import string
import subprocess
import tarfile
import tempfile
from pathlib import Path

from cluster_utils import extract_input, n_input_fastas
from snakemake.script import snakemake
from tqdm.auto import tqdm
from uniprot_utils import UniprotFastaParser

config = snakemake.config["treesapp_create"]
MUTE_TREESAPP = config["mute_treesapp"]
EXTRA_ARGS = []
for item in config["extra_args"]:
    EXTRA_ARGS += item.split()


def ref_pkg_name_generator(n: int = 6):
    characters = string.ascii_uppercase + string.digits
    for combo in itertools.product(characters, repeat=n):
        yield "".join(combo)


def run_treesapp(input_fasta: str, ref_pkg_name: str, output_dir: str) -> None:
    result = subprocess.run(
        [
            "treesapp",
            "create",
            "-i",
            input_fasta,
            "-c",
            ref_pkg_name,
            "-o",
            output_dir,
            "--headless",
        ]
        + EXTRA_ARGS,
        stdout=subprocess.PIPE if MUTE_TREESAPP else None,
        stderr=subprocess.PIPE if MUTE_TREESAPP else None,
    )

    # check result
    if result.returncode != 0:
        if MUTE_TREESAPP:
            print(result.stderr)
        raise RuntimeError(
            f"Got non-zero return code from TreeSAPP: {result.returncode}"
        )


if __name__ == "__main__":
    # inputs
    sequence_clusters_file = snakemake.input[0]

    # outputs
    hyperpackage_output = snakemake.output[0]

    name_generator = ref_pkg_name_generator()
    hyperpackage_manifest = {"ref_pkgs": []}

    expected_input_type = ".fasta.tar.gz"
    basename = os.path.basename(sequence_clusters_file)
    if basename.endswith(expected_input_type):
        hyperpackage_manifest["origin"] = basename[: -len(expected_input_type)]
    else:
        hyperpackage_manifest["origin"] = basename

    original_wd = os.getcwd()

    # create a temporary directories for output and running
    with tempfile.TemporaryDirectory() as output_dir, tempfile.TemporaryDirectory() as working_dir:

        print(
            "Note: TreeSAPP may take a long time to retrieve NCBI lineage info for large packages"
        )

        n_sequence_clusters = n_input_fastas(sequence_clusters_file)
        if n_sequence_clusters == 0:
            raise RuntimeError(f"No sequence clusters found in {n_sequence_clusters}")

        os.chdir(working_dir)
        for seq_cluster_filename, seq_cluster_file in tqdm(
            extract_input(os.path.join(original_wd, sequence_clusters_file)),
            desc="Running TreeSAPP create",
            total=n_sequence_clusters,
            unit="cluster",
        ):

            # write fasta file to disk
            with open(seq_cluster_filename, "wb") as f:
                f.write(seq_cluster_file.read())

            # get properties of hyperpackage
            cluster_name = Path(seq_cluster_filename).stem
            ref_pkg_name = next(name_generator)
            accession_ids = [
                accession for accession, _ in UniprotFastaParser(seq_cluster_filename)
            ]

            hyperpackage_manifest["ref_pkgs"].append(
                {
                    "ref_pkg_name": ref_pkg_name,
                    "exemplar_sequence": cluster_name,
                    "sequences": accession_ids,
                }
            )

            run_treesapp(
                seq_cluster_filename,
                ref_pkg_name,
                os.path.join(output_dir, ref_pkg_name),
            )

        os.chdir(original_wd)

        # write manifest file
        with open(os.path.join(output_dir, "manifest.json"), "w") as manifest_f:
            json.dump(hyperpackage_manifest, manifest_f)

        # put everything in the output archive
        with tarfile.open(hyperpackage_output, "w:gz") as output_tar:
            for file in os.listdir(output_dir):
                output_tar.add(
                    os.path.join(os.path.join(output_dir, file)), arcname=file
                )
