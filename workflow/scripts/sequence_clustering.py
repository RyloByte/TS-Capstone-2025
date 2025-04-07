import os
import subprocess
import tarfile
import tempfile

import pandas as pd
from Bio import SeqIO
from cluster_utils import filter_by_size, print_cluster_sizes, save_fasta_clusters
from snakemake.script import snakemake
from tqdm.auto import tqdm

config = snakemake.config["sequence_clustering"]
MUTE_MMSEQS = config["mute_mmseqs"]
MIN_CLUSTER_SIZE = config.get("min_cluster_size")
MAX_CLUSTER_SIZE = config.get("max_cluster_size")
MMSEQS_ARGS = []
for item in config["mmseqs_args"]:
    MMSEQS_ARGS += item.split()

mmseqs_output = "output"
mmseqs_tmp = "tmp"


def n_input_fastas(archive_path: str) -> int:
    i = 0
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar:
            if member.isfile() and (
                member.name.endswith(".fasta") or member.name.endswith(".faa")
            ):
                i += 1
    return i


def extract_input(archive_path: str):
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar:
            if member.isfile() and (
                member.name.endswith(".fasta") or member.name.endswith(".faa")
            ):
                yield member.name, tar.extractfile(member)


if __name__ == "__main__":
    # inputs
    structure_clusters = snakemake.input[0]

    # outputs
    output_archive = snakemake.output[0]

    clusters = {}
    original_wd = os.getcwd()

    # iterate through input
    n_structure_clusters = n_input_fastas(structure_clusters)
    for struct_cluster_filename, struct_cluster_file in tqdm(
        extract_input(structure_clusters),
        desc="Running mmseqs2 easy-linclust",
        total=n_structure_clusters,
        unit="cluster",
    ):

        # create temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)

            # write the fasta file
            with open(struct_cluster_filename, "wb") as f:
                f.write(struct_cluster_file.read())

            # get the records from that file
            records = {}
            for record in SeqIO.parse(struct_cluster_filename, "fasta"):
                uniprot_accession = record.id.split("|")[1]
                records[uniprot_accession] = record

            # run mmseqs
            result = subprocess.run(
                [
                    "mmseqs",
                    "easy-linclust",
                    struct_cluster_filename,
                    mmseqs_output,
                    mmseqs_tmp,
                ]
                + MMSEQS_ARGS,
                stdout=subprocess.PIPE if MUTE_MMSEQS else None,
                stderr=subprocess.PIPE if MUTE_MMSEQS else None,
                text=True,
            )

            # check result
            if result.returncode != 0:
                if MUTE_MMSEQS:
                    print(result.stderr)
                raise RuntimeError(
                    f"mmseqs easy-linclust failed with code {result.returncode}"
                )

            # save clusters
            result = pd.read_csv(
                os.path.join(tmpdir, f"{mmseqs_output}_cluster.tsv"),
                names=["rep", "mem"],
                sep="\t",
            )
            for cluster_name, cluster_df in result.groupby("rep"):
                clusters[cluster_name] = [
                    records[accession] for accession in cluster_df["mem"]
                ]

    os.chdir(original_wd)

    # print cluster sizes
    print_cluster_sizes(clusters)

    # filter by size
    clusters = filter_by_size(
        clusters, min_cluster_size=MIN_CLUSTER_SIZE, max_cluster_size=MAX_CLUSTER_SIZE
    )

    # create output
    save_fasta_clusters(clusters, output_archive)
