import os
import subprocess
import tarfile
from collections import defaultdict

from Bio import SeqIO
from snakemake.script import snakemake

MIN_SEQ_ID = snakemake.config["sequence_clustering"]["min_seq_id"]
MIN_SEQ_COV = snakemake.config["sequence_clustering"]["min_seq_cov"]
COV_MODE = snakemake.config["sequence_clustering"]["cov_mode"]
K_LENGTH = snakemake.config["sequence_clustering"]["k_length"]
KMER_PER_SEQ = snakemake.config["sequence_clustering"]["kmer_per_seq"]
THREADS = snakemake.config["sequence_clustering"]["threads"]
SHUFFLE = snakemake.config["sequence_clustering"]["shuffle"]
HASH_SHIFT = snakemake.config["sequence_clustering"]["hash_shift"]
REMOVE_TMP_FILES = snakemake.config["sequence_clustering"]["remove_tmp_files"]
FORCE_REUSE_TMP = snakemake.config["sequence_clustering"]["force_reuse_tmp"]
ALIGNMENT_MODE = snakemake.config["sequence_clustering"]["alignment_mode"]
SIMILARITY_TYPE = snakemake.config["sequence_clustering"]["similarity_type"]
REALIGN = snakemake.config["sequence_clustering"]["realign"]
SPACED_KMER_MODE = snakemake.config["sequence_clustering"]["spaced_kmer_mode"]
REMOVE_OUTPUT_FILES = snakemake.config["sequence_clustering"]["remove_output_files"]


def parse_cluster_tsv(cluster_file):
    cluster_map = defaultdict(list)

    with open(cluster_file, "r") as file:
        for line in file:
            cluster_id, sequence_id = line.strip().split("\t")
            cluster_map[cluster_id].append(sequence_id)

    return cluster_map


def extract_sequences(fasta_file):
    sequences = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = record

    return sequences


def write_cluster_files(cluster_map, sequences, output_dir, sample_name):
    os.makedirs(output_dir, exist_ok=True)

    for cluster_id, seq_ids in cluster_map.items():
        output_file = os.path.join(
            output_dir, f"{sample_name}-{cluster_id}-sequence_cluster.faa"
        )

        with open(output_file, "w") as f:
            for seq_id in seq_ids:
                if seq_id in sequences:
                    SeqIO.write(sequences[seq_id], f, "fasta")


def extract_tar_gz(archive_path, extract_to="./data/inputs"):
    with tarfile.open(archive_path, "r:gz") as tar:
        os.makedirs(extract_to, exist_ok=True)
        tar.extractall(path=extract_to)


def compress_cluster_files(output_tar_path, cluster_dir):
    with tarfile.open(output_tar_path, "w:gz") as tar:
        for file in os.listdir(cluster_dir):
            if file.endswith(".faa"):
                file_path = os.path.join(cluster_dir, file)
                tar.add(file_path, arcname=file)  # Add file to archive


# snakemake.input[0] = "data/{sample}-structure_clusters.tar.gz" , fill with {0..n}.faa
# snakemake.output[0] = "data/{sample}-sequence_clusters.tar.gz", fill with {0..n}.faa

if __name__ == "__main__":
    extract_path = "./data/inputs"
    extract_tar_gz(snakemake.input[0], extract_path)

    output_path = "./data/outputs"
    tmp_dir = "./data/tmp"
    tmp_cluster_dir = "./data/cluster_splits"
    final_output = f"./data"

    for file in os.listdir(extract_path):
        os.makedirs(output_path, exist_ok=True)
        os.makedirs(tmp_dir, exist_ok=True)
        os.makedirs(tmp_cluster_dir, exist_ok=True)

        if file.startswith("._"):
            continue
        filepath = os.path.join(extract_path, file)
        filename = file.rsplit("-", 1)[0]
        # subprocess.run(["mmseqs", "easy-linclust", os.path.join(extract_path, file), f"{output_path}/{filename}-sequence_clusters", tmp_dir, "--min-seq-id", f"{MIN_SEQ_ID}", "--min-seq-cov", f"{MIN_SEQ_COV}", "--cov-mode", f"{COV_MODE}", "-k", f"{K_LENGTH}", "--threads", f"{THREADS}", "shuffle", f"{SHUFFLE}"])
        subprocess.run(
            [
                "mmseqs",
                "easy-linclust",
                os.path.join(extract_path, file),
                f"{output_path}/{filename}-sequence_clusters",
                tmp_dir,
                "--min-seq-id",
                f"{MIN_SEQ_ID}",
                "-c",
                f"{MIN_SEQ_COV}",
                "--cov-mode",
                f"{COV_MODE}",
                "-k",
                f"{K_LENGTH}",
                "--kmer-per-seq",
                f"{KMER_PER_SEQ}",
                "--threads",
                f"{THREADS}",
                "--shuffle",
                f"{SHUFFLE}",
                # "--hash-shift", f"{HASH_SHIFT}",
                "--remove-tmp-files",
                f"{REMOVE_TMP_FILES}",
                "--force-reuse",
                f"{FORCE_REUSE_TMP}",
                "--alignment-mode",
                f"{ALIGNMENT_MODE}",
                "--similarity-type",
                f"{SIMILARITY_TYPE}",
                "--realign",
                f"{REALIGN}",
                "--spaced-kmer-mode",
                f"{SPACED_KMER_MODE}",
            ]
        )

        cluster_map = parse_cluster_tsv(
            f"{output_path}/{filename}-sequence_clusters_cluster.tsv"
        )
        sequences = extract_sequences(
            f"{output_path}/{filename}-sequence_clusters_all_seqs.fasta"
        )
        write_cluster_files(cluster_map, sequences, tmp_cluster_dir, filename)

    filename = snakemake.input[0].rsplit("/")[-1]
    filename = filename.rsplit("-")[0]
    final_output = os.path.join(final_output, f"{filename}-sequence_clusters.tar.gz")
    compress_cluster_files(final_output, tmp_cluster_dir)

    subprocess.run(
        f"rm -rf {tmp_dir} {tmp_cluster_dir} {extract_path}",
        shell=True,
        executable="/bin/bash",
    )
    if REMOVE_OUTPUT_FILES:
        subprocess.run(f"rm -rf {output_path}", shell=True, executable="/bin/bash")
