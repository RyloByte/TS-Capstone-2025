import pandas as pd
from snakemake.script import snakemake
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

logger = logging.getLogger()

reah_pattern = re.compile(r"(?:REAH:|reah:|)(\d{5})")
MAX_HOMOLOGOUS_PROTEINS = snakemake.config["max_homologous_proteins"]

# snakemake.input[0] = "data/{sample}-reahid.txt"
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
    swissprot_df = pd.read_csv(snakemake.input[1])  # TODO update with swissprot download step
    matching_sequences = swissprot_df[swissprot_df["REAH ID"] == reah_id]  # TODO update with swissprot download step
    logger.info(f"Found {len(matching_sequences)} sequences with matching REAH ID.")

    # check if none were found
    if len(matching_sequences) == 0:
        logger.critical(f"Did not find any sequences that match Reah ID {reah_id}.")
        raise RuntimeError(f"Did not find any sequences that match Reah ID {reah_id}.")

    # clip to max allowed
    logger.info(f"Clipping matching sequences to {MAX_HOMOLOGOUS_PROTEINS}.")
    matching_sequences = matching_sequences.head(MAX_HOMOLOGOUS_PROTEINS)

    # create list of sequences
    sequences = []
    for _, row in matching_sequences.iterrows():
        sequences.append(SeqRecord(Seq(row["sequence"]), id=row["uniprot_id"]))  # TODO update with swissprot download step

    # write list of sequences to file
    SeqIO.write(sequences, snakemake.output[0], "fasta")
