from Bio import SeqIO
import itertools
import os
import subprocess
from snakemake.script import snakemake

PDB_DIR = "data/pdb_structures"
CLUSTER_THRESHOLD = snakemake.config["structure_cluster_threshold"]


def get_rcsb_id(uniprot_id: str) -> str | None:
    pass

def get_alphafold_pdb(uniprot_id: str):
    pass

def get_rcsb_pdb(uniprot_id: str):
    pass

def pdb_pairs() -> list[tuple[str, str]]:
    # TODO unknown how the files come in
    return [(os.path.join(PDB_DIR, f1), os.path.join(PDB_DIR, f2)) for f1, f2 in itertools.combinations(snakemake.input, 2)]

def tm_score(pdb_1: str, pdb_2: str) -> float:
    result = subprocess.run(["tmalign", pdb_1, pdb_2], capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(":(")
    output = result.stdout
    # TODO regex the result and convert it to float
    

# snakemake.input[0] = "data/{sample}-homologous_proteins.faa"
# snakemake.output[0] = "data/{sample}-structure_clusters.tar.gz", fill with {0..n}.faa
if __name__ == "__main__":
    # load sequences
    uniprot_ids = [seq.name for seq in SeqIO.parse(snakemake.input[0], "fasta")]

    # download pdb files
    for uniprot_id in uniprot_ids:
        # check if there is already a pdb file for this protein
        if not os.path.exists(os.path.join(PDB_DIR, f"{uniprot_id}.pdb")):
            rcsb_id = get_rcsb_id(uniprot_id)
            if rcsb_id is None:
                get_alphafold_pdb(uniprot_id)
            else:
                get_rcsb_pdb(rcsb_id)

    # get pairwise pdb scores

    # create clusters

    # write to fasta files in data/structure_clusters
