import pandas as pd
import sys
import os
import gzip


def convert_lineage(lineage_str):
    rank_map = {
        "no rank": None, "clade": None,  # Omit these
        "cellular organisms": "r__",
        "superkingdom": "d__",
        "kingdom": None,
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__"
    }
    
    lineage_parts = lineage_str.split(", ")
    formatted_ranks = []
    
    for part in lineage_parts:
        if "(" not in part or ")" not in part:
            continue  # Skip malformed entries
        
        try:
            taxon, rank_name = part.rsplit(" (", 1)
            taxon = taxon.strip()
            rank_name = rank_name.rstrip(")").strip()

        except ValueError:
            continue  # Skip if it doesn't match expected format

        if rank_name in rank_map:
            prefix = rank_map[rank_name]
            if prefix is not None:  # Ignore omitted ranks
                formatted_ranks.append(f"{prefix}{taxon}")

    # Ensure all ranks exist in order, leaving blanks if missing
    rank_order = ["r__", "d__", "p__", "c__", "o__", "f__", "g__", "s__"]
    rank_dict = {r: "r__Root" if r == "r__" else r for r in rank_order}
    
    for entry in formatted_ranks:
        key = entry[:3]  # Extract prefix (e.g., "d__")
        rank_dict[key] = entry
    
    return ";".join(rank_dict.values()) + ";"

def process_tsv(input_file, dataframe):
    df = dataframe
    
    last_column = df.columns[-1]
    df["formatted_lineage"] = df[last_column].apply(lambda x: convert_lineage(str(x)))
    df = df[["Entry", "formatted_lineage"]]
    
    output_name = input_file.split("/", maxsplit=1)[1]
    output_name = f"utils/{output_name.split("_")[0]}"

    output_file = f"{output_name}_taxonomies.tsv"
    df.to_csv(output_file, sep="\t", index=False, header=None)
    print(f"Processed file saved as: {output_file}")

def process_gz_tsv(gz_file):
    with gzip.open(gz_file, "rt") as f:  # Open as text mode
        df = pd.read_csv(f, sep="\t")
        return df

if __name__ == "__main__":
    df = process_gz_tsv(snakemake.input[0])
    process_tsv(snakemake.input[0], df)