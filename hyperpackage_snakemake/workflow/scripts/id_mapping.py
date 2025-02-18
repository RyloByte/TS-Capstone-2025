from Bio import SeqIO
from snakemake.script import snakemake
import re


id_patterns = {
    "Uniprot/Swiss-Prot": re.compile(r"[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d([A-Z][A-Z0-9]{2}\d){1,2}"),
    "NCBI/RefSeq": re.compile(r"[A-Z]{2}_[A-Z0-9]+\.\d+"),
    "GenBank": re.compile(r"[A-Z]{2}\d{5}\d*"),
    "PDB": "",
    "TrEMBLE": ""
}



# snakemake.input[0] = "data/{sample}.faa"
# snakemake.output[0] = "data/{sample}-uniprot_mapped.faa"
if __name__ == "__main__":
    sequences = SeqIO.parse(snakemake.input[0], "fasta")
