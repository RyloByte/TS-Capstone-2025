from snakemake.script import snakemake
import subprocess
import tarfile
import tempfile
import os
from pathlib import Path
import logging
import string
import random


logger = logging.getLogger()


def generate_refpkg_name(n: int) -> str:
    chars = string.ascii_letters
    return "".join(random.choice(chars, k=n))


def run_treesapp(input_faa: str, refpkg_name: str, output_dir: str):
    logger.debug(f"Running TreeSAPP, input: {input_faa} output: {output_dir}")

    result = subprocess.run(
        [
            "treesapp",
            "create",
            "-i",
            input_faa,
            "-c",
            refpkg_name,
            "-o",
            output_dir,
            "--headless",
        ],
        capture_output=True,
        text=True,
    )

    if result.stdout:
        logger.debug(result.stdout.strip())
    if result.stderr:
        logger.warning(result.stderr.strip())

    if result.returncode != 0:
        raise RuntimeError(f"Got non-zero return code from TreeSAPP: {result.returncode}")


# snakemake.input[0] = "data/{sample}-sequence_clusters.tar.gz", fill with {0..n}.faa
# snakemake.output[0] = "result/{sample}-hyperpackage.tar.gz", fill with directories {0..n} for ref pkgs
if __name__ == "__main__":
    # create temp output
    with tempfile.TemporaryDirectory(delete=True) as temp_output:
        # open input tarfile
        with tarfile.open(snakemake.input[0], "r:gz") as tar:
            # go through each .faa file in tar
            for member in tar.getmembers():
                # write that tar to a tempfile
                with tempfile.NamedTemporaryFile(delete=True) as temp_fasta:
                    temp_fasta.write(tar.extractfile(member).read())
                    temp_fasta.flush()

                    # create output dir
                    refpkg_output = os.path.join(temp_output, Path(member).stem)
                    os.makedirs(refpkg_output, exist_ok=True)

                    # run treesapp
                    refpkg_name = generate_refpkg_name(6)
                    run_treesapp(temp_fasta.name, refpkg_name, refpkg_output)

        # gather the outputs into a tarfile
        with tarfile.open(snakemake.output[0], "w:gz") as tar:
            tar.add(temp_output, arcname=".")
