import gzip
import os
import sqlite3

import pandas as pd
from Bio import SeqIO
from snakemake.script import snakemake
from tqdm.auto import tqdm

config = snakemake.config["cluster_db"]
FILTER_BY_SPROT = config.get("filter_by_sprot", True)
CHUNK_SIZE = config.get("chunk_size", 10000000)


if __name__ == "__main__":
    # inputs
    foldseek_clusters_file = snakemake.input[0]
    sprot_fasta = snakemake.input[1]

    # outputs
    db_file = snakemake.output[0]

    # load swissprot accessions for filtering
    if FILTER_BY_SPROT:
        print("Loading SwissProt accessions...")
        with gzip.open(sprot_fasta, "rt") as f:
            swissprot_accessions = set(
                record.id.split("|")[1] for record in SeqIO.parse(f, "fasta")
            )
    else:
        print(
            "Building database with ALL FoldSeek entries... (this will take a long time)"
        )

    # remove existing database
    if os.path.exists(db_file):
        os.remove(db_file)

    conn = sqlite3.connect(db_file)

    # create table - built in index for primary key
    conn.execute(
        """
    CREATE TABLE IF NOT EXISTS clusters
    (repId TEXT, memId TEXT PRIMARY KEY)
    """
    )

    # known numbers from the specific file we're using
    n_lines = 214684310
    n_chunks = n_lines // CHUNK_SIZE + 1

    # read and insert chunks
    for chunk in tqdm(
        pd.read_csv(
            foldseek_clusters_file,
            names=["repId", "memId", "cluFlag", "taxId"],
            sep="\t",
            chunksize=CHUNK_SIZE,
            compression="gzip",
        ),
        total=n_chunks,
        desc="Building cluster database",
        unit="chunk",
    ):
        chunk = chunk.drop(columns=["cluFlag", "taxId"])
        if FILTER_BY_SPROT:
            chunk = chunk[chunk["memId"].isin(swissprot_accessions)]
        chunk.to_sql("clusters", conn, if_exists="append", index=False)

    row_count = conn.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]

    if FILTER_BY_SPROT:
        print(
            f"Found clusters for {row_count}/{len(swissprot_accessions)} ({row_count / len(swissprot_accessions):.1%}) SwissProt accessions"
        )
    else:
        print(f"Found clusters for {row_count} proteins")

    # commit changes
    conn.commit()
    conn.close()
