from snakemake.script import snakemake
import pandas as pd
from tqdm import tqdm
import gzip
from Bio import SwissProt


prefix_ranks = ("r__Root",)
ranks = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
columns = ("SeqID", "Lineage")


if __name__ == "__main__":
    # inputs
    annotated_sequences_path = snakemake.input[0]

    # outputs
    taxonomies_file = snakemake.output[0]

    # extract taxonomy info
    taxonomy_items = []
    with gzip.open(annotated_sequences_path, "rt") as handle:
        for record in tqdm(SwissProt.parse(handle), desc=f"Extracting SwissProt taxonomies", unit="entries", total=572970):
            entry_taxonomy = list(prefix_ranks)
            for prefix, rank in zip(ranks, record.organism_classification):
                entry_taxonomy.append(f"{prefix}{rank}")
            taxonomy_items.append((record.entry_name, "; ".join(entry_taxonomy)))

    # create and save table
    df = pd.DataFrame(taxonomy_items, columns=columns)
    df.to_csv(taxonomies_file, sep="\t", index=False, compression="gzip")
