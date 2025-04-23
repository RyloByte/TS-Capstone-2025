#!/bin/bash
set -e

mkdir -p data utils results

# check for/create user config file
if [ ! -f "config.yaml" ]; then
  cp "config.yaml.example" "config.yaml"
fi

source activate $ENV_NAME

exec snakemake --use-conda "$@"