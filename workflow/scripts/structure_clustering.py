import sqlite3
from collections import defaultdict

from Bio import SeqIO
from cluster_utils import filter_by_size, print_cluster_sizes, save_fasta_clusters
from snakemake.script import snakemake

config = snakemake.config["structure_clustering"]
MIN_CLUSTER_SIZE = config["min_cluster_size"]
MAX_CLUSTER_SIZE = config["max_cluster_size"]


if __name__ == "__main__":
    # inputs
    fasta_file = snakemake.input[0]
    cluster_db = snakemake.input[1]

    # outputs
    output_archive = snakemake.output[0]

    # get uniprot accessions from sequences
    records = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        uniprot_accession = record.id.split("|")[1]
        records[uniprot_accession] = record
    if len(records) == 0:
        raise RuntimeError(f"No sequence records found in {fasta_file}")

    # query cluster db for accessions
    print(f"Querying structure cluster db for {len(records)} sequences...")
    conn = sqlite3.connect(cluster_db)
    placeholders = ",".join(["?"] * len(records))
    query = f"SELECT repId, memId FROM clusters WHERE memId in ({placeholders})"

    results = conn.execute(query, list(records.keys())).fetchall()
    conn.close()

    print(f"Found clusters for {len(results)}/{len(records)} sequences")

    # group clusters
    clusters = defaultdict(list)
    for repId, memId in results:
        clusters[repId].append(records[memId])

    # print cluster sizes
    print_cluster_sizes(
        clusters, min_cluster_size=MIN_CLUSTER_SIZE, max_cluster_size=MAX_CLUSTER_SIZE
    )

    # filter by size
    clusters = filter_by_size(
        clusters, min_cluster_size=MIN_CLUSTER_SIZE, max_cluster_size=MAX_CLUSTER_SIZE
    )

    # create output
    save_fasta_clusters(clusters, output_archive)
