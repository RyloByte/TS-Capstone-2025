import os
import tarfile
import tempfile
from collections import Counter

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def print_cluster_sizes(clusters: dict, max_width: int = 80) -> None:
    cluster_sizes = Counter(len(cluster) for cluster in clusters.values())
    most_popular_cluster = max(cluster_sizes.values())
    multiplier = (
        1 if most_popular_cluster <= max_width else max_width / most_popular_cluster
    )
    print("Cluster size / count")
    for cluster_size, cluster_count in sorted(cluster_sizes.items()):
        print(
            f"{cluster_size:3}: {'█' * max(1, int(cluster_count * multiplier))} {cluster_count}"
        )


def filter_by_size(
    clusters: dict,
    min_cluster_size: int | None = None,
    max_cluster_size: int | None = None,
) -> dict:
    if min_cluster_size is not None:
        clusters = {k: v for k, v in clusters.items() if len(v) >= min_cluster_size}
    if max_cluster_size is not None:
        clusters = {k: v for k, v in clusters.items() if len(v) <= max_cluster_size}

    print(f"Keeping {len(clusters)} clusters", end="")
    if min_cluster_size is not None:
        print(f" >= {min_cluster_size}", end="")
    if max_cluster_size is not None:
        print(f" <= {max_cluster_size}", end="")
    print()

    return clusters


def save_fasta_clusters(clusters: dict[str, SeqRecord], output_archive: str) -> None:
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
