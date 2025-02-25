from Bio import SeqIO
from snakemake.script import snakemake
import re


def get_uniprot_from_ncbi(ncbi_id: str) -> str | None:
    pass  # TODO

def get_uniprot_from_genbank(genbank_id: str) -> str | None:
    pass  # TODO

def get_uniprot_from_pdb(pdb_id: str) -> str | None:
    pass  # TODO

def get_uniprot_from_tremble(tremble_id: str) -> str | None:
    pass  # TODO

uniprot_id_pattern = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

id_patterns = {
    # "Uniprot/Swiss-Prot": re.compile(r"sp|(.*)|"),
    "NBCI/RefSeq": re.compile(r"ref|(.*)|"),
    "GenBank": re.compile(r"gb|(.*)|"),
    "PDB": re.compile(r"pdb|(.*)|"),
    "TrEMBLE": re.compile(r"tr|(.*)|")
}

conversion_functions = {
    "NBCI/RefSeq": get_uniprot_from_ncbi,
    "GenBank": get_uniprot_from_genbank,
    "PDB": get_uniprot_from_pdb,
    "TrEMBLE": get_uniprot_from_tremble
}

def classify_id(seq_id: str) -> tuple[str, str] | None:
    for format_name, pattern in id_patterns.items():
        id_match = pattern.match(seq_id)
        if id_match is not None:
            return format_name, id_match.group()
        
def extract_uniprot_id(seq_id: str) -> str | None:
    id_match = uniprot_id_pattern.search(seq_id)    
    if id_match is not None:
        id_match = id_match.group()
    return id_match


# snakemake.input[0] = "data/{sample}.faa"
# snakemake.output[0] = "data/{sample}-uniprot_mapped.faa"
if __name__ == "__main__":
    raise NotImplementedError("ID mapping has not been implemented yet. If creating a reference package from proteins, please provided a data/{sample}-uniprot_mapped.faa with proteins that are annotated with UniProt Accession IDs. If you are creating a reference package from a Reah function, please provide a data/{sample}-reahid.txt file.")

    # load sequences
    sequences = list(SeqIO.parse(snakemake.input[0], "fasta"))
    for seq in sequences:
        # try to extract a uniprot id from the sequence id
        uniprot_id = extract_uniprot_id(seq.id)

        if uniprot_id is None:
            # get id type
            id_info = classify_id(seq.id)
            if id_info is None:
                # did not recognize the id type
                raise RuntimeError(f"Could not determine id type for: {seq.id}")
            # convert id
            uniprot_id = conversion_functions[id_info[0]](id_info[1])
            if uniprot_id is None:
                # conversion failed
                raise RuntimeError(f"Could not get Uniprot Accession ID for {id_info[0]} ID: {id_info[1]}")
        
        # set id to uniprot id
        seq.id = uniprot_id
        # clear out the description (TODO is this necessary?)
        seq.description = ""
    SeqIO.write(sequences, snakemake.output[0], "fasta")
