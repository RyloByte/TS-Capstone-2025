#!/bin/bash

set -e

IMAGE_NAME="treesapp-hyperpackage-workflow"

docker build -t $IMAGE_NAME .

echo "Built image: $IMAGE_NAME"