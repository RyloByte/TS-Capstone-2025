from snakemake.script import snakemake
from Bio import SeqIO


# snakemake.input = ["data/structure_clusters/0.faa", ... "data/structure_clusters/{n}.faa"]
# output = "data/sequence_clusters/{0..n}.faa"
if __name__ == "__main__":
    for file in snakemake.input:
        sequences = SeqIO.parse(file, "fasta")
