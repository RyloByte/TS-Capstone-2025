from snakemake.script import snakemake
from Bio import SeqIO


# snakemake.input[0] = "data/input.faa"
# output: "data/uniprot_mapped.faa"
if __name__ == "__main__":
    sequences = SeqIO.parse(snakemake.input[0], "fasta")
