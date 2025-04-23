#!/bin/bash

set -e

IMAGE_NAME="ghcr.io/RyloByte/TS-Capstone-2025:latest"

# get number of processors
if command -v nproc &> /dev/null; then
  NUM_PROCS=$(nproc)
else
  NUM_PROCS=$(sysctl -n hw.ncpu)
fi

# run deleting container with mounts
docker run \
  --rm \
  -v "$(pwd)/data:/workflow-dir/data" \
  -v "$(pwd)/utils:/workflow-dir/utils" \
  -v "$(pwd)/results:/workflow-dir/results" \
  -v "$(pwd)/config.yaml:/workflow-dir/config.yaml" \
  -it $IMAGE_NAME \
  -j $NUM_PROCS "$@"