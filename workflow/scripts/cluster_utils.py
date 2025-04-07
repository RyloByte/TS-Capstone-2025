import os
import tarfile
import tempfile
from collections import Counter

from Bio import SeqIO


def print_cluster_sizes(
    clusters: dict, max_width: int = 80, min_cluster_size=None, max_cluster_size=None
) -> None:
    cluster_sizes = Counter(len(cluster) for cluster in clusters.values())
    most_popular_cluster = max(cluster_sizes.values())
    multiplier = (
        1 if most_popular_cluster <= max_width else max_width / most_popular_cluster
    )
    print("Cluster size / count")
    printed_min, printed_max = False, False
    for cluster_size, cluster_count in sorted(cluster_sizes.items()):
        if (
            not printed_min
            and min_cluster_size is not None
            and cluster_size >= min_cluster_size
        ):
            print(f"=========== MIN CLUSTER SIZE: {min_cluster_size} ===========")
            printed_min = True
        if (
            not printed_max
            and max_cluster_size is not None
            and cluster_size > max_cluster_size
        ):
            print(f"=========== MAX CLUSTER SIZE: {max_cluster_size} ===========")
            printed_max = True
        print(
            f"{cluster_size:3}: {'█' * max(1, int(cluster_count * multiplier))} {cluster_count}"
        )


def filter_by_size(
    clusters: dict,
    min_cluster_size=None,
    max_cluster_size=None,
) -> dict:
    if min_cluster_size is not None:
        clusters = {k: v for k, v in clusters.items() if len(v) >= min_cluster_size}
    if max_cluster_size is not None:
        clusters = {k: v for k, v in clusters.items() if len(v) <= max_cluster_size}

    print(f"Keeping {len(clusters)} clusters", end="")
    p_and = False
    if min_cluster_size is not None:
        print(f" >= {min_cluster_size}", end="")
        p_and = True
    if max_cluster_size is not None:
        if p_and:
            print(" AND", end="")
        print(f" <= {max_cluster_size}", end="")
    print(f" ({sum(len(cluster) for cluster in clusters.values())} total sequences)")

    return clusters


def save_fasta_clusters(clusters: dict, output_archive: str) -> None:
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


def n_input_fastas(archive_path: str) -> int:
    i = 0
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar:
            if member.isfile() and (
                member.name.endswith(".fasta") or member.name.endswith(".faa")
            ):
                i += 1
    return i


def extract_input(archive_path: str):
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar:
            if member.isfile() and (
                member.name.endswith(".fasta") or member.name.endswith(".faa")
            ):
                yield member.name, tar.extractfile(member)
