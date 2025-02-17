Create environment: `conda env create -f config/environment.yaml`

Activate environment: `conda activate snakemake_env`

Run workflow: `snakemake --use-conda --cores 4`

Delete intermediate files: `snakemake --delete-all-output`
