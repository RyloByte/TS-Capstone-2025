from snakemake.script import snakemake
from Bio import SeqIO
import requests

UNIPROT_API = "https://"

# snakemake.input[0] = "data/uniprot_mapped.faa"
# output = "data/homologous_proteins.faa"
if __name__ == "__main__":
    # load sequences from input
    sequences = SeqIO.parse(snakemake.input[0], "fasta")

    # query uniprot for each of them

    # write them to a file
