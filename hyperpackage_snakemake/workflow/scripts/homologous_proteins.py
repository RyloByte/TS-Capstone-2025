from snakemake.script import snakemake
from Bio import SeqIO
import pandas as pd

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
                rhea2uniprot[r_id].append(uniprot_id)
    
    return uniprot2rhea, rhea2uniprot


# snakemake.input[0] = "data/{sample}-uniprot_mapped.faa"
# snakemake.input[1] = "utils/uniprot_spdb/uniprot_sprot.fasta"
# snakemake.input[2] = "utils/uniprot_spdb/uniprotkb_AND_reviewed_true_2025_02_22.tsv"
# snakemake.output[0] = "data/{sample}-homologous_proteins.faa"
if __name__ == "__main__":
    # load sequences from input
    sequences = SeqIO.parse(snakemake.input[0], "fasta")

    # load uniprot-rhea dicts
    uniprot2rhea, rhea2uniprot = parse_uniprot_mapper(snakemake.input[2])

    homologous_rhea = set()
    for seq in sequences:
        u_id = seq.id
        r_ids = uniprot2rhea.get(u_id)
        for r_id in r_ids:
            homologous_rhea.add(r_id)
        print(u_id)
        print(r_ids)
        
    print(len(homologous_rhea))
    homologous_ids = set()

    for r_id in homologous_rhea:
        u_ids = rhea2uniprot.get(r_id)
        for u_id in u_ids:
            homologous_ids.add(u_id)

    print(len(homologous_ids))

    # load all sp sequences from sp database
    # sp_sequencest = SeqIO.parse(snakemake.input[1], "fasta")
    # write them to a file
