Create environment: `conda env create -f environment.yaml`

Activate environment: `conda activate snakemake_env`

Run workflow for input `data/{sample}.faa` or `data/{sample}.fasta`: `snakemake --use-conda --cores 4 result/{sample}-hyperpackage.tar.gz`

Delete intermediate files: `snakemake --delete-all-output` (might delete your input `data/{sample}.faa` file?)
