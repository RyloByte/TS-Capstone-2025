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

# def get_experimental_pdb(uniprot_map: pd.DataFrame, uniprot_id: str):
#     pdb_ids = uniprot_map[(uniprot_map["Entry"] == uniprot_id) & (~uniprot_map["PDB"].isna())] 
#     if len(pdb_ids) == 0:
#         return None
#     pdb_ids = uniprot_map["PDB"].unique()[0]
#     pdb_ids = [x for x in pdb_ids.strip().split(";") if x != ""]
#     print(f"EXP: {pdb_ids}")
#     print(f"Found {len(pdb_ids)} experimental PDBs, selecting first structure...")
#     return pdb_ids[0]

# def get_alphafold_pdb(uniprot_map: pd.DataFrame, uniprot_id: str):
#     pdb_ids = uniprot_map[(uniprot_map["Entry"] == uniprot_id) & (~uniprot_map["AlphaFoldDB"].isna())] 
#     if len(pdb_ids) == 0:
#         return None
#     pdb_ids = pdb_ids["AlphaFoldDB"].unique()[0]
#     pdb_ids = [x for x in pdb_ids.strip().split(";") if x != ""]
#     print(f"AF: {pdb_ids}")
#     print(f"Found {len(pdb_ids)} AlphaFoldDB PDBs, selecting first structure...")
#     return pdb_ids[0]

def create_cluster_output(cluster_df, output_tar_file):    
    # Write each cluster to a separate .faa file
    os.makedirs(os.path.join("data", "foldseek_cluster"), exist_ok=True)

    faa_files = []
    for i, cluster_id in enumerate(cluster_df["cluster_id"].unique()):
        # Define the output file name
        output_file = os.path.join("data", "foldseek_cluster", f"cluster_{i}.faa")
        faa_files.append(output_file)
        
        # Write sequences to the file
        with open(output_file, "w", newline="\n") as f:
            n_members = 0
            for member_ids in cluster_df[cluster_df["cluster_id"] == cluster_id]["member_id"].unique():
                if n_members > MAX_SEQS:
                    continue
                if SP_ONLY:
                    found_sequence, record = get_sp_status(member_ids)
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
    homologous_seqs = [x.id.split("|")[1] for x in list(SeqIO.parse(snakemake.input[0], "fasta"))]
    print(homologous_seqs)

    foldseek_clusters = pd.read_csv(snakemake.input[1], sep="\t", names=["memberID", "repID", "taxID"], compression="gzip")

    seq_clusters = {}
    cluster_ids = set()

    # find all clusters
    with gzip.open(snakemake.input[1], "rt") as infile:
        for line in infile:
            repID, memID, cluFlag, taxID = line.strip("\n").split("\t")
            if memID in homologous_seqs:
                if cluFlag == 3 or cluFlag == 4:
                    print(f"{memID} identified as singleton, skipping...")
                else:
                    print(f"{memID} identified in {repID} cluster!")
                    if memID not in seq_clusters:
                        seq_clusters[memID] = {}
                        seq_clusters[memID]["repID"] = []

                    seq_clusters[memID]["taxID"] = taxID
                    seq_clusters[memID]["repID"].append(repID)
                    cluster_ids.add(repID)

    print(f"Identified {len(seq_clusters)}/{len(homologous_seqs)} as non-singletons comprising {len(cluster_ids)} clusters, extracting cluster information...")

    # iter through file again and get all ids 
    cluster_lst = []
    member_lst = []
    tax_lst = []

    with gzip.open(snakemake.input[1], "rt") as infile:
        for line in infile:
            repID, memID, cluFlag, taxID = line.strip("\n").split("\t")                
            if repID in cluster_ids:
                cluster_lst.append(repID)
                member_lst.append(memID)
                tax_lst.append(taxID)

    cluster_out = pd.DataFrame({"cluster_id": cluster_lst, "member_id": member_lst, "taxonomy_id": tax_lst})
    
    print(f"Found {len(cluster_out)} structually similar proteins! Retrieving fasta sequences...")

    create_cluster_output(cluster_out, snakemake.output[0])
        
