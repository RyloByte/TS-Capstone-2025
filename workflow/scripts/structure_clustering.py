from snakemake.script import snakemake
from Bio import SeqIO
import pandas as pd
import os
import gzip
import subprocess
import shutil
import tarfile
import requests

RES_DIR = "data/foldseek_results"
CLUST_DIR = "data/foldseek_cluster"

MAX_SEQS = snakemake.config["foldseek_max_seqs"]
SP_ONLY = snakemake.config["foldseek_swissprot_only"]

def create_cluster_output(cluster_df, homologous_seqs, output_tar_file):    
    # Write each cluster to a separate .faa file
    os.makedirs(os.path.join("data", "foldseek_cluster"), exist_ok=True)

    faa_files = []
    for i, cluster_id in enumerate(cluster_df["repID"].unique()):
        # Define the output file name
        output_file = os.path.join("data", "foldseek_cluster", f"cluster_{cluster_id}.faa")
        faa_files.append(output_file)
        
        # Write sequences to the file
        with open(output_file, "w", newline="\n") as f:

            # first, save fastas for original homologous sequences 
            for h_seq in homologous_seqs:
                if h_seq in cluster_df[cluster_df["repID"] == cluster_id]["memID"].unique():
                    SeqIO.write(homologous_seqs[h_seq], f, "fasta")

            # then get sequences for other clusters
            n_members = 0
            for member_ids in cluster_df[cluster_df["repID"] == cluster_id]["memID"].unique():
                if n_members > MAX_SEQS - 1:
                    continue
                if SP_ONLY:
                    found_sequence, record =  (member_ids)
                    if found_sequence:
                        SeqIO.write(record, f, "fasta")
                        n_members += 1
                else:
                    member_fasta = get_fasta_from_uniprot(member_ids)
                    f.write(member_fasta)
                    n_members += 1    
            print(f"Found {n_members} qualifying sequences for cluster {cluster_id}")            
                
    with tarfile.open(output_tar_file, "w:gz") as tar:
        for file in faa_files:
            tar.add(file, arcname=os.path.basename(file))

    print(f"Saved structure clusters to {output_tar_file}!")

def get_fasta_from_uniprot(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error: Unable to fetch data for {uniprot_id}. Status code: {response.status_code}")
        return None

def get_sp_status(uniprot_id):
    found_sequence = False
    with gzip.open(snakemake.input[2], "rt") as f:
        records = list(SeqIO.parse(f, "fasta"))
        for record in records:
            if uniprot_id in record.id:
                found_sequence = True
                return found_sequence, record
    return found_sequence, ""

# snakemake.input[0] = "data/{sample}-homologous_proteins.faa"
# snakemake.input[1] = "utils/foldseek_clusters/foldseek_clusters.tsv.gz"
# snakemake.input[2] = "utils/swissprot_sequences.fasta.gz"
# snakemake.output[0] = "data/{sample}-structure_clusters.tar.gz", fill with {0..n}.faa
if __name__ == "__main__":
    # load sequences
    homologous_seqs = {x.id.split("|")[1]: x for x in list(SeqIO.parse(snakemake.input[0], "fasta"))}
    print(f"Detected the following Rhea homologous proteins: {list(homologous_seqs.keys())}")

    print("Loading Foldseek cluster database... (this will take a moment)")

    foldseek_clusters = pd.read_csv(snakemake.input[1], sep="\t", names=["repID", "memID", "cluFlag", "taxID"], compression="gzip")

    # filter to remove singletons
    foldseek_clusters = foldseek_clusters[(~foldseek_clusters["cluFlag"] != 3) & (~foldseek_clusters["cluFlag"] != 4)]
    print(f"Loaded Foldseek cluster database containing {len(foldseek_clusters)} proteins in {foldseek_clusters['repID'].nunique()} non-singleton clusters!")
    
    # get all cluster ids matching homologous sequences\
    cluster_ids = set()

    for h_seq in homologous_seqs:
        h_df = foldseek_clusters[foldseek_clusters["memID"] == h_seq]
        if len(h_df) == 0:
            print(f"ERROR: {h_seq} not found in non-singleton cluster, continuing...")
        else:
            cluster_ids.add(h_df["repID"].unique()[0])
    
    print(f"Identified {len(cluster_ids)} non-singleton clusters, extracting cluster information...")

    # get all proteins in matching cluster_ids
    seq_clusters = foldseek_clusters[foldseek_clusters["repID"].isin(cluster_ids)]
    seq_clusters = seq_clusters.drop(columns=["cluFlag"])
    
    print(f"Found {len(seq_clusters)} structually similar proteins! Retrieving fasta sequences...")

    if SP_ONLY:
        print("NOTE: Skipping all non-swiss-prot proteins in clusters!")

    create_cluster_output(seq_clusters, homologous_seqs, snakemake.output[0])