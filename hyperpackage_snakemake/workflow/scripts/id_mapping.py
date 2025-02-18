from Bio import SeqIO
from snakemake.script import snakemake


# snakemake.input[0] = "data/{sample}.faa"
# snakemake.output[0] = "data/{sample}-uniprot_mapped.faa"
if __name__ == "__main__":
    sequences = SeqIO.parse(snakemake.input[0], "fasta")
