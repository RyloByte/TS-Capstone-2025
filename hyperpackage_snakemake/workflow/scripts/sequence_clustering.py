from snakemake.script import snakemake
from Bio import SeqIO
import os
import tarfile
import subprocess


CLUSTER_THRESHOLD = snakemake.config["sequence_cluster_threshold"]

def extract_tar_gz(archive_path, extract_to="./data/inputs"):
    with tarfile.open(archive_path, "r:gz") as tar:
        os.makedirs(extract_to, exist_ok=True)
        tar.extractall(path=extract_to)

# snakemake.input[0] = "data/{sample}-structure_clusters.tar.gz" , fill with {0..n}.faa
# snakemake.output[0] = "data/{sample}-sequence_clusters.tar.gz", fill with {0..n}.faa
if __name__ == "__main__":
    extract_path = "./data/inputs"
    extract_tar_gz(snakemake.input[0], extract_path)

    output_path = "./data/outputs"
    os.makedirs(output_path, exist_ok=True)

    tmp_dir = "./data/tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    for file in os.listdir(extract_path):
        if file.startswith("._"):
            continue
        filepath = os.path.join(extract_path, file)
        filename = file.rsplit("-", 1)[0]
        subprocess.run(["mmseqs", "easy-linclust", os.path.join(extract_path, file), f"{output_path}/{filename}-sequence_clusters.faa", tmp_dir])

    # for file in os.listdirs()

    # for file in 
        # file = tar.extractfile(f"{faa_file}-structure_clusters.faa")
        # content = file.read().decode()
        # print(content)
        # paths = tar.getnames()  # Lists all files inside
    
    # for file in paths:
    #     if file.split(
    # for file in os.listdir(f"data/{filename}"):
    #     print(file)

