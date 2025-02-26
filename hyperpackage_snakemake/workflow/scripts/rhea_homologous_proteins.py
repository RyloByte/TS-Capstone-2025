import pandas as pd
from snakemake.script import snakemake
import re
from Bio import SeqIO
import logging
import gzip

logger = logging.getLogger()

rhea_pattern = re.compile(r"(?<!\d)(?:RHEA:|rhea:|)(\d{5})(?!\d)")
MAX_HOMOLOGOUS_PROTEINS = snakemake.config["max_homologous_proteins"]

# snakemake.input = ["data/{sample}-rheaid.txt", "utils/swissprot_data.tsv.gz", "utils/swissprot_sequences.fasta.gz"]
# snakemake.output[0]= "data/{sample}-uniprot_mapped.faa"
if __name__ == "__main__":
    # load the Rhea ID
    with open(snakemake.input[0], "r") as f:
        rhea_id = re.search(rhea_pattern, f.read())
    if rhea_id is None:
        raise RuntimeError(
            f"Did not find Rhea ID in {snakemake.input[0]}, ex. RHEA:12345 or 67890"
        )
    rhea_id = rhea_id.group()
    logger.info(f"Found Rhea ID: {rhea_id}")

    # get the matching proteins
    swissprot_df = pd.read_csv(snakemake.input[1], sep="\t", compression="gzip")
    matching_sequences = swissprot_df[swissprot_df["Rhea ID"].str.contains(rhea_id, na=False, case=False)]
    logger.info(f"Found {len(matching_sequences)} sequences with matching Rhea ID")

    # check if none were found
    if len(matching_sequences) == 0:
        raise RuntimeError(f"Did not find any sequences that match Rhea ID {rhea_id}")

    # clip to max allowed
    logger.info(f"Clipping matching sequences to {MAX_HOMOLOGOUS_PROTEINS}")
    matching_sequences = matching_sequences.head(MAX_HOMOLOGOUS_PROTEINS)

    # find the sequence for each accession ID
    # painfully inefficient but... ¯\_(ツ)_/¯
    output_sequences = []
    with gzip.open(snakemake.input[2], "rt") as f:
        records = list(SeqIO.parse(f, "fasta"))
    for accession_id in matching_sequences["Entry"]:
        found_sequence = False
        for record in records:
            if accession_id in record.id:
                output_sequences.append(record)
                found_sequence = True
                break
        if not found_sequence:
            logger.error(f"Did not find sequence for accession ID {accession_id}")

    # write list of sequences to file
    SeqIO.write(output_sequences, snakemake.output[0], "fasta")
