import pandas as pd
from snakemake.script import snakemake
import re
from Bio import SeqIO
import logging
import gzip

logger = logging.getLogger()

ec_pattern = re.compile(r"(?:EC:|ec:|)((?:\d+)(?:\.(?:\d+)){0,2}(?:\.(?:n?\d+))?)")
MAX_HOMOLOGOUS_PROTEINS = snakemake.config["max_homologous_proteins"]

# snakemake.input[0] = ["data/{sample}-ecnumber.txt", "utils/swissprot_data.tsv.gz", "utils/swissprot_sequences.fasta.gz"]
# snakemake.output[0]= "data/{sample}-uniprot_mapped.faa"
if __name__ == "__main__":
    # load the EC Number
    with open(snakemake.input[0], "r") as f:
        ec_number = re.search(ec_pattern, f.read())
    if ec_number is None:
        raise RuntimeError(f"Did not find EC number in {snakemake.input[0]}, ex. EC:1.2.3.4 or 4.5 etc.")
    ec_number = ec_number.group()
    logger.info(f"Found EC number: {ec_number}.")

    # get the matching proteins
    swissprot_df = pd.read_csv(snakemake.input[1], sep="\t", compression="gzip")
    number_and_subnumber_matcher = re.compile(rf"(?<!\d|\.){ec_number}")
    matching_sequences = swissprot_df[swissprot_df["EC number"].str.contains(number_and_subnumber_matcher)]
    logger.info(f"Found {len(matching_sequences)} sequences with compatible EC number.")

    # check if none were found
    if len(matching_sequences) == 0:
        raise RuntimeError(f"Did not find any sequences with EC number: {ec_number}.")

    # clip to max allowed
    logger.info(f"Clipping matching sequences to {MAX_HOMOLOGOUS_PROTEINS}.")
    matching_sequences = matching_sequences.head(MAX_HOMOLOGOUS_PROTEINS)

    # find the sequence for each accession ID
    # painfully inefficient but... ¯\_(ツ)_/¯
    output_sequences = []
    with gzip.open(snakemake.input[2], "r") as f:
        for accession_id in matching_sequences["Entry"]:
            found_sequence = False
            for record in SeqIO.parse(accession_id, "fasta"):
                if accession_id in record.id:
                    output_sequences.append((accession_id, record.seq))
                    found_sequence = True
                    break
            if not found_sequence:
                logger.error(f"Did not find sequence for accession ID {accession_id}.")

    # write list of sequences to file
    SeqIO.write(output_sequences, snakemake.output[0], "fasta")
