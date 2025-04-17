#!/bin/bash

IMAGE_NAME="treesapp-hyperpackage-workflow"

# allow compiling for x64 from arm
docker run --privileged --rm tonistiigi/binfmt --install all

# check for build x
docker buildx create --name multiarch --use 2>/dev/null || true

docker buildx build \
  --platform linux/amd64 \
  -t $IMAGE_NAME \
  --load .

echo "Built image: $IMAGE_NAME"