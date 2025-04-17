FROM condaforge/mambaforge:latest

ENV ENV_NAME="snakemake_env"
ENV TRIGGER_ALL="results/hyperpackages/ec_1.refpkg.tar.gz"

WORKDIR "/workflow-dir"

COPY . .

# create root environment
RUN mamba env create -f environment.yaml --quiet

# create other environments
RUN cp config.yaml.example config.yaml && \
    /bin/bash -c "source activate $ENV_NAME && snakemake --use-conda --conda-create-envs-only -j $(nproc) $TRIGGER_ALL" && \
    rm config.yaml.example config.yaml

RUN mv entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]