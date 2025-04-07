import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from cluster_utils import (
    extract_input,
    filter_by_size,
    n_input_fastas,
    print_cluster_sizes,
    save_fasta_clusters,
)
from snakemake.script import snakemake
from tqdm.auto import tqdm

config = snakemake.config["sequence_clustering"]
MUTE_MMSEQS = config["mute_mmseqs"]
MIN_CLUSTER_SIZE = config["min_cluster_size"]
MAX_CLUSTER_SIZE = config["max_cluster_size"]
MMSEQS_ARGS = []
for item in config["mmseqs_args"]:
    MMSEQS_ARGS += item.split()

mmseqs_output = "output"
mmseqs_tmp = "tmp"


def run_mmseqs(input_fasta: str, output_name: str, temp_name: str) -> None:
    result = subprocess.run(
        [
            "mmseqs",
            "easy-linclust",
            input_fasta,
            output_name,
            temp_name,
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
        raise RuntimeError(f"mmseqs easy-linclust failed with code {result.returncode}")


if __name__ == "__main__":
    # inputs
    structure_clusters = snakemake.input[0]

    # outputs
    output_archive = snakemake.output[0]

    clusters = {}
    original_wd = os.getcwd()

    # iterate through input
    n_structure_clusters = n_input_fastas(structure_clusters)
    if n_structure_clusters == 0:
        raise RuntimeError(f"No structure clusters found in {structure_clusters}")
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
            run_mmseqs(struct_cluster_filename, mmseqs_output, mmseqs_tmp)

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
    print_cluster_sizes(
        clusters, min_cluster_size=MIN_CLUSTER_SIZE, max_cluster_size=MAX_CLUSTER_SIZE
    )

    # filter by size
    clusters = filter_by_size(
        clusters, min_cluster_size=MIN_CLUSTER_SIZE, max_cluster_size=MAX_CLUSTER_SIZE
    )

    # create output
    save_fasta_clusters(clusters, output_archive)
