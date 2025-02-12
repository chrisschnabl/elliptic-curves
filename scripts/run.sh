#!/bin/bash
set -e;

# Navigate to the directory of the script
cd "$(dirname "$0")"

# Build the base docker
docker build -f Dockerfile.base -t base .

# Build the example docker
docker build -f Dockerfile.example -t example .

# Run the container
docker run --rm example
