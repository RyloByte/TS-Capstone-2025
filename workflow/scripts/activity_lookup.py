import gzip
import json
import random

import pandas as pd
from Bio import SeqIO, SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from snakemake.script import snakemake
from tqdm import tqdm
import re


# ec_num_pattern = re.compile(r"EC[ |:=]*((\d+\.)*\d+)", flags=re.IGNORECASE)
# rhea_pattern = re.compile(r"RHEA[ |:=]*(\d+)", flags=re.IGNORECASE)


# config
config = snakemake.config["activity_lookup"]
MAX_HOMOLOGOUS_PROTEINS = config["max_proteins"]


if __name__ == "__main__":
    # wildcards
    activity_number_type = snakemake.wildcards.number_type.lower()
    activity_number = snakemake.wildcards.number.lower()

    # inputs
    annotated_sequences_path = snakemake.input[0]

    # outputs
    fasta_output_file = snakemake.output[0]
    taxonomy_output_file = snakemake.output[1]

    # create pattern
    pattern = re.compile(rf"{re.escape(activity_number_type)}[ |:=]*{re.escape(activity_number)}(?!\d)", flags=re.IGNORECASE)

    print(pattern)

    # iterate through records
    records = []
    taxonomy_items = []
    with gzip.open(annotated_sequences_path, "rt") as handle:
        for record in tqdm(SwissProt.parse(handle), desc=f"Searching SwissProt for {activity_number_type.upper()} {activity_number}", unit="entries", total=572970):
            # search description
            use_record = bool(re.search(pattern, record.description))

            # search comments
            use_record = use_record or any(bool(re.search(pattern, comment)) for comment in record.comments)

            # check cross-references
            if not use_record:
                for cross_ref in record.cross_references:
                    if cross_ref[0].lower() == activity_number_type and re.match(pattern, activity_number_type + cross_ref[1]):
                        use_record = True
                        break

            # append record
            if use_record:
                records.append(SeqRecord(
                    Seq(record.sequence),
                    id=record.entry_name,
                    description="buns",
                ))

                taxonomy_items.append((
                    record.entry_name,
                    "Root; " + "; ".join(record.organism_classification)
                ))

    if len(records) == 0:
        raise ValueError(f"No sequences matching {activity_number_type} {activity_number} found")
    print(f"Found {len(records)} sequences")

    # trim if necessary
    if len(records) > MAX_HOMOLOGOUS_PROTEINS:
        print(f"Trimming {len(records)} sequences to {MAX_HOMOLOGOUS_PROTEINS}")
        indices = random.sample(range(len(records)), MAX_HOMOLOGOUS_PROTEINS)
        records = [records[i] for i in indices]
        taxonomy_items = [taxonomy_items[i] for i in indices]

    # write sequences
    with open(fasta_output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")

    # write taxonomy
    df = pd.DataFrame(taxonomy_items, columns=["SeqID", "Lineage"])
    df.to_csv(taxonomy_output_file, sep="\t", index=False)
