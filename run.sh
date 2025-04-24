#!/bin/bash

set -e

IMAGE_NAME="ghcr.io/RyloByte/TS-Capstone-2025:latest"
LOCAL_IMAGE_NAME="treesapp-hyperpackage-workflow"

if docker pull $IMAGE_NAME >/dev/null 2>&1; then
  echo "Using image from GitHub Container Registry"
  RUN_IMAGE=$IMAGE_NAME
else
  echo "Attempting to use locally built image"
  RUN_IMAGE=$LOCAL_IMAGE_NAME
fi

# make necessary directories
mkdir -p data utils results

# get number of processors
if command -v nproc &> /dev/null; then
  NUM_PROCS=$(nproc)
else
  NUM_PROCS=$(sysctl -n hw.ncpu)
fi

# run deleting container with mounts
if [ -f "config.yaml" ]; then
  docker run \
    --rm \
    -v "$(pwd)/data:/workflow-dir/data" \
    -v "$(pwd)/utils:/workflow-dir/utils" \
    -v "$(pwd)/results:/workflow-dir/results" \
    -v "$(pwd)/config.yaml:/workflow-dir/config.yaml" \
    -it $RUN_IMAGE \
    -j $NUM_PROCS "$@"
else
  docker run \
    --rm \
    -v "$(pwd)/data:/workflow-dir/data" \
    -v "$(pwd)/utils:/workflow-dir/utils" \
    -v "$(pwd)/results:/workflow-dir/results" \
    -it $RUN_IMAGE \
    -j $NUM_PROCS "$@"
fi