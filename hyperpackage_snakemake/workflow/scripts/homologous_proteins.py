from snakemake.script import snakemake
from Bio import SeqIO
from typing import Optional
import pandas as pd
import re

MAX_HOMOLOGOUS_THRESHOLD = snakemake.config["max_homologous_sequences"]

def parse_uniprot_mapper(mapper_path: str) -> tuple[dict, dict]:
    uniprot_map = pd.read_csv(mapper_path, sep="\t")

    uniprot2rhea = {}
    rhea2uniprot = {}

    for i, row in uniprot_map.iterrows():
        uniprot_id = row["Entry"]
        rhea_ids = set(row["Rhea ID"].split(" ")) if isinstance(row["Rhea ID"], str) else None
        uniprot2rhea[uniprot_id] = rhea_ids
        if rhea_ids: # if uniprot has rhea id(s), add to dict
            for r_id in rhea_ids:
                if r_id not in rhea2uniprot:
                    rhea2uniprot[r_id] = set()
                rhea2uniprot[r_id].add(uniprot_id)
    
    return uniprot2rhea, rhea2uniprot

uniprot_id_pattern = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

def extract_uniprot_id(seq_id: str) -> Optional[str]:
    id_match = uniprot_id_pattern.search(seq_id)    
    if id_match is not None:
        id_match = id_match.group()
    return id_match


# snakemake.input[0] = "data/{sample}-uniprot_mapped.faa"
# snakemake.input[1] = "utils/uniprot_spdb/uniprot_sprot.fasta"
# snakemake.input[2] = "utils/uniprot_spdb/uniprotkb_AND_reviewed_true_2025_02_22.tsv"
# snakemake.output[0] = "data/{sample}-homologous_proteins.faa"
if __name__ == "__main__":
    # load sequences from input
    sequences = list(SeqIO.parse(snakemake.input[0], "fasta"))

    # load uniprot-rhea dicts
    uniprot2rhea, rhea2uniprot = parse_uniprot_mapper(snakemake.input[2])

    homologous_rhea = set()
    for seq in sequences:
        u_id = extract_uniprot_id(seq.id)
        r_ids = uniprot2rhea.get(u_id)
        if r_ids is not None:
            for r_id in r_ids:
                homologous_rhea.add(r_id)
        
    homologous_ids = set()
    for r_id in homologous_rhea:
        u_ids = rhea2uniprot.get(r_id)
        if u_ids is not None:
            for u_id in u_ids:
                homologous_ids.add(u_id)

    if len(homologous_ids) > MAX_HOMOLOGOUS_THRESHOLD:
        raise RuntimeError(f"Input faa produces {len(homologous_ids)} homologous proteins, max threshold is {MAX_HOMOLOGOUS_THRESHOLD}")

    # load all swiss-prot sequences from swiss-prot database
    sp_sequences = list(SeqIO.parse(snakemake.input[1], "fasta"))
    filtered_sp_sequences = []
    for sp_seq in sp_sequences:
        sp_id = extract_uniprot_id(sp_seq.id)
        if sp_id in homologous_ids:
            sp_seq.id = sp_id
            sp_seq.description = ""
            filtered_sp_sequences.append(sp_seq)

    print(f"Found total of {len(filtered_sp_sequences)} homologous proteins from {len(sequences)} input proteins (mapped to {len(homologous_rhea)} total Rhea IDs)")
    # write homologous protein sequences to file
    SeqIO.write(filtered_sp_sequences, snakemake.output[0], "fasta")
