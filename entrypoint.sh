#!/bin/bash

# Exit on error, treat unset vars as errors, and fail on pipeline errors
set -euo pipefail

echo "Sourcing Conda and activating environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate snakemake_env

echo "Running Snakemake pipeline with x86-64 enforced..."
CONDA_SUBDIR=linux-64 snakemake --use-conda --conda-frontend conda --cores 8 data/hyperpackages/rhea_10596.refpkg.tar.gz