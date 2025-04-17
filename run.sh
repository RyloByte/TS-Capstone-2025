#!/bin/bash

IMAGE_NAME="treesapp-hyperpackage-workflow"

mkdir -p data utils results
if [ ! -f "config.yaml" ]; then
  echo "Creating config.yaml from example file"
  cp "config.yaml.example" "config.yaml"
fi

docker run \
  --platform linux/amd64 \
  -v "$(pwd)/data:/workflow-dir/data" \
  -v "$(pwd)/utils:/workflow-dir/utils" \
  -v "$(pwd)/results:/workflow-dir/results" \
  -v "$(pwd)/config.yaml:/workflow-dir/config.yaml" \
  -it $IMAGE_NAME \
  -j $(nproc) "$@"