#!/bin/bash
set -e

source activate $ENV_NAME

exec snakemake --use-conda "$@"