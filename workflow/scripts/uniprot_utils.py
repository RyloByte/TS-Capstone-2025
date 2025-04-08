import gzip
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def extract_uniprot_accession(record) -> str:
    if isinstance(record, str):
        return record.split("|")[1]
    if isinstance(record, SeqRecord):
        return record.id.split("|")[1]
    raise TypeError("Expected string or SeqRecord")


class UniprotFastaParser:
    def __init__(self, fasta_file):
        fasta_file = Path(fasta_file)
        if fasta_file.suffix == ".gz":
            self.handle = gzip.open(fasta_file, "rt")
            self.parser = SeqIO.parse(self.handle, "fasta")
            result = subprocess.run(
                [f"zcat {fasta_file} | grep -c '^>'"],
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0:
                self.n_records = int(result.stdout.strip())
            else:
                self.n_records = None
        else:
            self.handle = None
            self.parser = SeqIO.parse(fasta_file, "fasta")
            result = subprocess.run(
                [f"grep -c '^>' {fasta_file}"],
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0:
                self.n_records = int(result.stdout.strip())
            else:
                self.n_records = None

    def __del__(self):
        if self.handle is not None:
            self.handle.close()

    def __iter__(self):
        return self

    def __next__(self):
        record = next(self.parser)
        accession = extract_uniprot_accession(record)
        return accession, record

    def __len__(self):
        return self.n_records


# def parse_uniprot_fasta(fasta_path):
#     if fasta_path.endswith(".gz"):
#         f = gzip.open(fasta_path, "rt")
#     else:
#         f = fasta_path
#     for record in SeqIO.parse(f, "fasta"):
#         accession = extract_uniprot_accession(record)
#         yield accession, record
#
#     if fasta_path.endswith(".gz"):
#         f.close()
