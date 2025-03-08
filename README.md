# Installation:

All the commands in this readme assume working directory of `hyperpackage_create/`.

1. Clone repository/download workflow code
2. Create conda environment with `conda env create -f environment.yaml`
3. Activate environment with `conda activate snakemake_env`

Update conda environment with `conda env update --file environment.yaml`

# Usage:

1. Place a `.fasta` or `.faa` file in `data/`, ex. `data/DsrAB.faa`
2. Request a finished file with `snakemake --use-conda result/DsrAB-hyperpackage.tar.gz`

There are other inputs that you can give, and you can also request intermediate files.

# Flowchart (slightly out of date):

```mermaid
flowchart TD
    A{{download_swissprot}}
    B[data/swissprot.db]
    A-->B
    C{{rename_fasta}}
    D(["data/{sample}.fasta"])
    E(["data/{sample}.faa"])
    D-->C
    C-->E
    F{{id_mapping}}
    G(["data/{sample}-uniprot_mapped.faa"])
    E-->F
    F-->G
    H{{homologous_proteins}}
    I(["data/{sample}-rheaid.txt"])
    B-->H
    I-->H
    H-->G
    J{{structure_clustering}}
    K["data/{sample}-structure_clusters.tar.gz"]
    G-->J
    J-->K
    L{{sequence_clustering}}
    M["data/{sample}-sequence_clusters.tar.gz"]
    K-->L
    L-->M
    N{{treesapp_create}}
    O["result/{sample}-hyperpackage.tar.gz"]
    B-->N
    M-->N
    N-->O
    P{{rename_mapped_fasta}}
    Q(["data/{sample}-uniprot_mapped.fasta"])
    Q-->P
    P-->G
```

# Tests

## Creating a test

The tests for this project are configured to run dynamically based on things in the directories

- `tests/successful_scenarios/`
- `tests/failure_scenarios/`
- `tests/expected_outputs/`

To create a test create a folder in `tests/successful_scenarios/` (or `tests/failure_scenarios/`) with the name of your
test, and a folder with the same name in `tests/expected_outputs/` with the expected output files (if you are making a
failing test you can use empty files since they will just be used to request outputs).

Let's make a test for the following rule:

```snakemake
rule rhea_homologous_proteins:
    input:
        "input/{sample}-rhea_id.txt",
        "utils/swissprot_data.tsv.gz",
        "utils/swissprot_sequences.fasta.gz"
    output:
        "data/{sample}-uniprot_mapped.faa"
    script:
        "scripts/rhea_homologous_proteins.py"
```

We can make `tests/successful_scenarios/rhea_test/` and fill it with the necessary inputs:

```
tests/successful_scenarios/rhea_test/
|-- input/
|   |-- whatever-rhea_id.txt
```

and `tests/expected_outputs/rhea_test/` with the files we want to request and check for accuracy:

```
tests/expected_outputs/rhea_test/
|-- data/
|   |-- whatever-uniprot_mapped.faa
```

If `config` and `util` directories are not provided, they will be linked from the main
directory. If your rule depends on a config value, it is best to include a config in your scenario so that it runs as
expected every time.

When a test runs it will:

1. Call snakemake requesting all the files in your expected output
2. Make sure snakemake exited normally
3. Check that each output file is identical to the expected one

If you are making a scenario for a workflow that will fail, it will:

1. Call snakemake requesting all the files in your expected output (these can be empty files)
2. Make sure snakemake *did not* exit normally
3. Make sure none of the output files were actually produced

## Running tests

When you have the conda environment activated and you are in the correct working directory you can run all tests by
running `pytest`.

You can run tests in parallel with `pytest -n auto`. You could go to funky town if you have a situation like you have
not downloaded the SwissProt DB, and then you have multiple tests wanting to download it at once. 
