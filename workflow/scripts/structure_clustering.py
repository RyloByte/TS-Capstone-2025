import sqlite3
from collections import defaultdict

from cluster_utils import filter_by_size, print_cluster_sizes, save_fasta_clusters
from snakemake.script import snakemake
from tqdm.auto import tqdm
from uniprot_utils import UniprotFastaParser

config = snakemake.config["structure_clustering"]
MIN_CLUSTER_SIZE = config["min_cluster_size"]
MAX_CLUSTER_SIZE = config["max_cluster_size"]


if __name__ == "__main__":
    # inputs
    fasta_file = snakemake.input[0]
    cluster_db = snakemake.input[1]
    sprot_fasta = snakemake.input[2]

    # outputs
    output_archive = snakemake.output[0]

    # get uniprot accessions from sequences
    records = {
        accession: record for accession, record in UniprotFastaParser(fasta_file)
    }
    if len(records) == 0:
        raise RuntimeError(f"No sequence records found in {fasta_file}")

    # query cluster db for accessions
    print(f"Querying structure cluster db for {len(records)} sequences in {fasta_file}...")
    conn = sqlite3.connect(cluster_db)
    placeholders = ",".join(["?"] * len(records))
    query = f"SELECT DISTINCT repId FROM clusters WHERE memId in ({placeholders})"

    cluster_ids = [cluster_id[0] for cluster_id in conn.execute(query, list(records.keys())).fetchall()]

    placeholders = ",".join(["?"] * len(cluster_ids))
    query = f"SELECT repId, memId FROM clusters WHERE repId in ({placeholders})"
    seqs_with_clusters = conn.execute(query, cluster_ids).fetchall()

    conn.close()

    print(f"Found {len(cluster_ids)} clusters with {len(seqs_with_clusters)} sequences")

    # group clusters
    clusters = defaultdict(list)
    reverse_map = {memId: repId for repId, memId in seqs_with_clusters}

    # iterate through all uniprot sequences
    for accession, record in tqdm(
        UniprotFastaParser(sprot_fasta),
        desc="Retrieving Cluster Items from SwissProt",
        unit="records",
    ):
        if accession in reverse_map:
            cluster_id = reverse_map.pop(accession)
            clusters[cluster_id].append(record)
            if len(reverse_map) == 0:
                break

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
