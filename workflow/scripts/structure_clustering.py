from snakemake.script import snakemake
from Bio import SeqIO
import pandas as pd
import os
import subprocess
import shutil
import tarfile

PDB_DIR = "data/pdb_structures"
RES_DIR = "data/foldseek_results"
CLUST_DIR = "data/foldseek_cluster"

CLUSTER_THRESHOLD = snakemake.config["structure_cluster_threshold"]
E_THRESHOLD = snakemake.config["foldseek_eval_threshold"]
N_ITERS = snakemake.config["foldseek_num_iterations"]
MAX_SEQS = snakemake.config["foldseek_max_seqs"]

def get_experimental_pdb(uniprot_map: pd.DataFrame, uniprot_id: str):
    pdb_ids = uniprot_map[(uniprot_map["Entry"] == uniprot_id) & (~uniprot_map["PDB"].isna())] 
    if len(pdb_ids) == 0:
        return None
    pdb_ids = uniprot_map["PDB"].unique()[0]
    pdb_ids = [x for x in pdb_ids.strip().split(";") if x != ""]
    print(f"EXP: {pdb_ids}")
    print(f"Found {len(pdb_ids)} experimental PDBs, selecting first structure...")
    return pdb_ids[0]

def get_alphafold_pdb(uniprot_map: pd.DataFrame, uniprot_id: str):
    pdb_ids = uniprot_map[(uniprot_map["Entry"] == uniprot_id) & (~uniprot_map["AlphaFoldDB"].isna())] 
    if len(pdb_ids) == 0:
        return None
    pdb_ids = pdb_ids["AlphaFoldDB"].unique()[0]
    pdb_ids = [x for x in pdb_ids.strip().split(";") if x != ""]
    print(f"AF: {pdb_ids}")
    print(f"Found {len(pdb_ids)} AlphaFoldDB PDBs, selecting first structure...")
    return pdb_ids[0]

def create_cluster_output(all_seqs_fasta, output_tar_file):
    sequences = list(SeqIO.parse(all_seqs_fasta, "fasta"))
    clusters = {}
    current_cluster = None
    for seq in sequences:
        if len(str(seq.seq)) == 0:
            current_cluster = seq.id
            clusters[current_cluster] = []
        else:
            clusters[current_cluster].append(seq)
    
    # Write each cluster to a separate .faa file
    faa_files = []
    for i, cluster_id in enumerate(clusters):
        # Define the output file name
        output_file = os.path.join("data", "foldseek_cluster", f"cluster_{i}.faa")
        faa_files.append(output_file)
        # Write sequences to the file
        seqs = clusters[cluster_id]
        with open(output_file, "w", newline="\n") as f:
            SeqIO.write(seqs, f, "fasta")
    
    with tarfile.open(output_tar_file, "w:gz") as tar:
        for file in faa_files:
            tar.add(file, arcname=os.path.basename(file))


# snakemake.input[0] = "data/{sample}-homologous_proteins.faa"
# snakemake.input[1] = "utils/foldseek_spdb/foldseek_spdb"
# snakemake.input[2] = "utils/uniprot_spdb/uniprotkb_AND_reviewed_true.tsv"
# snakemake.output[0] = "data/{sample}-structure_clusters.tar.gz", fill with {0..n}.faa
if __name__ == "__main__":
    # load sequences
    homologous_seqs = list(SeqIO.parse(snakemake.input[0], "fasta"))

    uniprot_map = pd.read_csv(snakemake.input[2], sep="\t")
    print(uniprot_map.columns)
    print(uniprot_map.head())

    # download pdb files
    os.makedirs(PDB_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)
    os.makedirs(CLUST_DIR, exist_ok=True)

    for h_seq in homologous_seqs:
        h_id = h_seq.id.split("|")[1]
        print(h_id)
        # check if there is already a pdb file for this protein
        if not os.path.exists(os.path.join(PDB_DIR, f"{h_id}.pdb")):
            if h_id not in uniprot_map["Entry"].unique():
                print(f"Warning: {h_id} not found in Swiss-Prot database, skipping...")
                continue
            pdb_id = get_experimental_pdb(uniprot_map, h_id)
            
            if pdb_id is None:
                pdb_id = get_alphafold_pdb(uniprot_map, h_id)
                if pdb_id is None:
                    print(f"Warning: {h_id} does not have any structures, skipping...")
                    continue
                else:
                    subprocess.run(["wget", f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb", "-O", f"{PDB_DIR}/{h_id}.pdb"])
            else:
                subprocess.run(["wget", f"https://files.rcsb.org/download/{pdb_id}.pdb", "-O", f"{PDB_DIR}/{h_id}.pdb"])
        else:
            print(f"Found PDB_DIR/{h_id}.pdb, skipping download")

    for pdb in os.listdir(PDB_DIR):
        if os.path.basename(pdb).split(".")[0] in [h_seq.id.split("|")[1] for h_seq in homologous_seqs]:
            subprocess.run([
                "foldseek", 
                "easy-search", 
                f"{PDB_DIR}/{pdb}", 
                snakemake.input[1], 
                f"{RES_DIR}/{os.path.splitext(pdb)[0]}", 
                "tmp",
                "--format-output", "query,target,alntmscore,lddt,prob",
                "-e", str(E_THRESHOLD),
                "--num-iterations", str(N_ITERS),
                "--max-seqs", str(MAX_SEQS),
                ])
    
    for pdb_res in os.listdir(RES_DIR):
        res_df = pd.read_csv(os.path.join(RES_DIR, pdb_res), sep="\t", names=["query", "target", "alntmscore", "lddt", "prob"])

        for i, af_pdb in enumerate(res_df["target"].unique()):
            af_id = af_pdb.split("-")[1]
            if i > MAX_SEQS:
                continue
            if not os.path.exists(os.path.join(PDB_DIR, f"{af_id}.pdb")):
                subprocess.run(["wget", f"https://alphafold.ebi.ac.uk/files/{af_pdb}.pdb", "-O", f"{PDB_DIR}/{af_id}.pdb"])


    print(f"FOUND {len(os.listdir(PDB_DIR))} STRUCTURALLY SIMILAR PROTEINS, NOW CLUSTERING...")

    subprocess.run([
        "foldseek", "easy-cluster", PDB_DIR, os.path.join(CLUST_DIR, "fs_cluster_results"), "tmp",
        "--alignment-type", "1",
        "-e", str(E_THRESHOLD),])
 
    create_cluster_output(os.path.join(CLUST_DIR, "fs_cluster_results_all_seqs.fasta"), snakemake.output[0])

    shutil.rmtree("tmp")

        
