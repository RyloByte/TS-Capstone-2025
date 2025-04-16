FROM ubuntu:22.04

# Non-interactive shell
ENV DEBIAN_FRONTEND=noninteractive

# Set up Conda
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
ENV CONDA_SUBDIR=linux-64

# Install required system packages
RUN apt-get update && apt-get install -y \
    wget \
    git \
    vim \
    bzip2 \
    ca-certificates \
    curl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Miniconda (x86_64)
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# Double-lock conda architecture (optional safety)
RUN conda config --system --set subdir linux-64

# Create workspace
WORKDIR /workspace

# Copy environment definition and create env
COPY environment.yaml .
RUN conda env create -f environment.yaml

# Copy entrypoint script
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'set -euo pipefail' >> /entrypoint.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /entrypoint.sh && \
    echo 'conda activate snakemake_env' >> /entrypoint.sh && \
    echo 'echo "Running Snakemake..."' >> /entrypoint.sh && \
    echo 'CONDA_SUBDIR=linux-64 snakemake --use-conda --conda-frontend conda --cores 8 data/hyperpackages/rhea_10596.refpkg.tar.gz' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

# Ensure bash is used in RUN and ENTRYPOINT
SHELL ["/bin/bash", "-c"]

# Set default entrypoint
ENTRYPOINT ["/entrypoint.sh"]