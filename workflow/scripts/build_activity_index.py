from snakemake.script import snakemake
from Bio import SwissProt
import gzip
import json
import re
from collections import defaultdict
from tqdm import tqdm


ec_num_pattern = re.compile(r"EC[ |:=]*((\d+\.)*\d+)", flags=re.IGNORECASE)
rhea_pattern = re.compile(r"RHEA[ |:=]*(\d+)", flags=re.IGNORECASE)


if __name__ == "__main__":
    # inputs
    annotated_sequences = snakemake.input[0]

    # outputs
    index_file = snakemake.output[0]

    # iterate sequences
    ec_index, rhea_index = defaultdict(list), defaultdict(list)
    with gzip.open(annotated_sequences, "rt") as handle:
        for record in tqdm(SwissProt.parse(handle), desc="Building activity number index", unit="entries", total=572970):
            entry_name = record.entry_name

            # search comments
            for comment in record.comments:
                rhea_matches = re.findall(rhea_pattern, comment)
                for rhea_match in rhea_matches:
                    rhea_index[rhea_match].append(entry_name)
                ec_num_matches = re.findall(ec_num_pattern, comment)
                for ec_num_match in ec_num_matches:
                    ec_index[ec_num_match[0]].append(entry_name)

            # search description
            rhea_matches = re.findall(rhea_pattern, record.description)
            for rhea_match in rhea_matches:
                rhea_index[rhea_match].append(entry_name)
            ec_num_matches = re.findall(ec_num_pattern, record.description)
            for ec_num_match in ec_num_matches:
                ec_index[ec_num_match[0]].append(entry_name)

    # write index to file
    combined_index = {
        "EC": ec_index,
        "RHEA": rhea_index,
    }
    with open(index_file, "w", newline="\n") as f:
        json.dump(combined_index, f)
