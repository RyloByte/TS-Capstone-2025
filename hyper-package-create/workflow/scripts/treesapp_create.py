from snakemake.script import snakemake
import subprocess


def run_treesapp(input_faa: str):
    result = subprocess.run(["treesapp", input_faa])
    if result.returncode != 0:
        raise Exception(":(")


# snakemake.input = ["data/squence_clusters/0.faa", ... "data/sequence_clusters/{n}.faa"]
# output = "data/reference_packages/{0..n}/"
if __name__ == "__main__":
    for file in snakemake.input:
        pass
