import os
import sqlite3
import tarfile
import tempfile
from collections import defaultdict, Counter

from Bio import SeqIO
from snakemake.script import snakemake


config = snakemake.config["structure_clustering"]
MIN_CLUSTER_SIZE = config["min_cluster_size"]


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
    cluster_sizes = Counter(len(cluster) for cluster in clusters.values())
    max_width = 80
    most_popular_cluster = max(cluster_sizes.values())
    multiplier = 1 if most_popular_cluster <= max_width else max_width / most_popular_cluster
    print("Cluster size / count")
    for cluster_size, cluster_count in sorted(cluster_sizes.items()):
        print(f"{cluster_size:3}: {'█' * max(1, int(cluster_count * multiplier))} {cluster_count}")

    # filter by size
    clusters = {k: v for k, v in clusters.items() if len(v) >= MIN_CLUSTER_SIZE}
    print(f"Keeping {len(clusters)} clusters larger than {MIN_CLUSTER_SIZE}")

    # create output
    with tarfile.open(output_archive, "w:gz") as tar:
        for cluster_id, cluster in clusters.items():
            with tempfile.NamedTemporaryFile(suffix="fasta", delete=False) as tmp:
                try:
                    SeqIO.write(cluster, tmp.name, "fasta")
                    tmp.close()
                    tar.add(tmp.name, arcname=f"{cluster_id}.fasta")
                finally:
                    if os.path.exists(tmp.name):
                        os.unlink(tmp.name)
