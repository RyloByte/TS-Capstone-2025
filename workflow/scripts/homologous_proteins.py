import gzip
import os

import pandas as pd
from Bio import SeqIO
from snakemake.script import snakemake
import re


rhea_pattern = re.compile(r"\d{5}")
ec_pattern = re.compile(r"(\d+)(?:\.(\d+)){0,2}(?:\.(n?\d+))?")

config = snakemake.config["homologous_lookup"]
MAX_HOMOLOGOUS_PROTEINS = config["max_proteins"]

def get_rhea_matching_proteins(rhea_id: str, swissprot_df: pd.DataFrame) -> pd.DataFrame:
    return swissprot_df[swissprot_df["Rhea ID"].str.contains(rhea_id, na=False, case=False)]

def get_ec_matching_proteins(ec_num: str, swissprot_df: pd.DataFrame) -> pd.DataFrame:
    return swissprot_df[
        swissprot_df["EC number"].str.contains(rf"(?<!\d|\.){re.escape(ec_num)}", na=False)
    ]


if __name__ == "__main__":
    activity_number = snakemake.wildcards.number

    swissprot_df = pd.read_csv(snakemake.input[0], sep="\t", compression="gzip")

    output_file_name = os.path.basename(snakemake.output[0]).lower()

    if output_file_name.startswith("rhea_"):
        if re.fullmatch(rhea_pattern, activity_number) is None:
            raise RuntimeError(f"{activity_number} is not a valid rhea activity number, ex. 12345")
        matching_proteins = get_rhea_matching_proteins(activity_number, swissprot_df)
        print(f"Found {len(matching_proteins)} proteins matching rhea ID '{activity_number}'")
    elif output_file_name.startswith("ec_"):
        if re.fullmatch(ec_pattern, activity_number) is None:
            raise RuntimeError(f"{activity_number} is not a valid ec activity number, ex. 1.2.3.4 or 5.6.7 etc.")
        matching_proteins = get_ec_matching_proteins(activity_number, swissprot_df)
        print(f"Found {len(matching_proteins)} proteins matching EC number '{activity_number}'")
    else:
        raise RuntimeError(f"{output_file_name} should have started with 'rhea_' or 'ec_'")

    if len(matching_proteins) == 0:
        raise RuntimeError(f"Did not find any homologous proteins for activity number {output_file_name.split('-')[0]}")

    print(f"Clipping matching proteins to {MAX_HOMOLOGOUS_PROTEINS}")
    matching_proteins = matching_proteins.head(MAX_HOMOLOGOUS_PROTEINS)

    output_sequences = []
    with gzip.open(snakemake.input[1], "rt") as f:
        records = list(SeqIO.parse(f, "fasta"))
    for accession_id in matching_proteins["Entry"]:
        found_sequence = False
        for record in records:
            if accession_id in record.id:
                output_sequences.append(record)
                found_sequence = True
                break
        if not found_sequence:
            print(f"Did not find sequence for accession ID {accession_id}")

    # write list of sequences to file
    with open(snakemake.output[0], "w", newline="\n") as f:
        SeqIO.write(output_sequences, f, "fasta")
