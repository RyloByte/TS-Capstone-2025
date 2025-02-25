import pandas as pd
from snakemake.script import snakemake
import re
from Bio import SeqIO
import logging
import gzip

logger = logging.getLogger()

reah_pattern = re.compile(r"(?:REAH:|reah:|)(\d{5})")
MAX_HOMOLOGOUS_PROTEINS = snakemake.config["max_homologous_proteins"]

# snakemake.input = ["data/{sample}-reahid.txt", "utils/swissprot_data.tsv.gz", "utils/swissprot_sequences.fasta.gz"]
# snakemake.output[0]= "data/{sample}-uniprot_mapped.faa"
if __name__ == "__main__":
    # load the ReahID
    with open(snakemake.input[0], "r") as f:
        reah_id = re.search(reah_pattern, f.read())
    if reah_id is None:
        logger.critical(f"Did not find Reah ID in {snakemake.input[0]}, ex. REAH:12345 or 67890")
        raise RuntimeError(f"Did not find Reah ID in {snakemake.input[0]}, ex. REAH:12345 or 67890")
    reah_id = reah_id.group(0)
    logger.info(f"Found Reah ID: {reah_id}.")

    # get the matching proteins
    swissprot_df = pd.read_csv(snakemake.input[1], sep="\t", compression="gzip")
    matching_sequences = swissprot_df[swissprot_df["Reah ID"].str.contains(reah_id)]
    logger.info(f"Found {len(matching_sequences)} sequences with matching REAH ID.")

    # check if none were found
    if len(matching_sequences) == 0:
        logger.critical(f"Did not find any sequences that match Reah ID {reah_id}.")
        raise RuntimeError(f"Did not find any sequences that match Reah ID {reah_id}.")

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
