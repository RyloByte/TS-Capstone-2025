from snakemake.script import snakemake
import os
import subprocess


def run_treesapp(input_faa: str):
    result = subprocess.run(["treesapp", input_faa])
    if result.returncode != 0:
        raise Exception(":(")


# snakemake.input[0] = "data/{sample}-sequence_clusters", fill with {0..n}.faa
# snakemake.output[0] = "data/{sample}-referencce_packages", fill with directories {0..n} for ref pkgs
if __name__ == "__main__":
    for file in os.listdir(snakemake.input[0]):
        run_treesapp(os.path.join(snakemake.input[0], file))
