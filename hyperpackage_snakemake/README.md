Create environment: `conda env create -f environment.yaml`

Activate environment: `conda activate snakemake_env`

Run workflow for input `data/{sample}.faa` or `data/{sample}.fasta`: `snakemake --use-conda --cores 4 result/{sample}-hyperpackage.tar.gz`

Delete intermediate files: `snakemake --delete-all-output` (might delete your input `data/{sample}.faa` file?)

Installation:

1. Clone repository/download workflow code
2. Create conda environment with `conda env create -f environment.yaml`
3. Activate environment with `conda activate snakemake_env`

Usage:

1. Place a `.fasta` or `.faa` file in `data/`, ex. `data/DsrAB.faa`
2. Request a finished file with `snakemake --use-conda result/DsrAB-hyperpackage.tar.gz`

Flowchart:

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
