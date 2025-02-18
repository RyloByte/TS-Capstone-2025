from snakemake.script import snakemake
from Bio import SeqIO
import os


CLUSTER_THRESHOLD = snakemake.config["sequence_cluster_threshold"]


# snakemake.input[0] = "data/{sample}-structure_clusters", fill with {0..n}.faa
# snakemake.output[0] = "data/{sample}-sequence_clusters", fill with {0..n}.faa
if __name__ == "__main__":
    for file in os.listdir(snakemake.input[0]):
        sequences = SeqIO.parse(os.path.join(snakemake.input[0], file), "fasta")
