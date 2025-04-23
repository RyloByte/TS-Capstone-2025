#!/bin/bash

IMAGE_NAME="treesapp-hyperpackage-workflow"

docker built -t $IMAGE_NAME .

echo "Built image: $IMAGE_NAME"