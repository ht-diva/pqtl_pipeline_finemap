FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="c491f9e40f666e3716468cb047b869cfd8a80767f513b18b196e0e4a503c280f"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/fine_mapping.yml
#   prefix: /conda-envs/cfb439855512bc80795457ea259e5135
#   name: fine_mapping
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#     - R
#   dependencies:
#     - plink2=2.00a5.10-0
#     - gcta=1.94.1-0
#     - r-base=4.3.3
#     - r-optparse=1.7.4
#     - r-tidyverse=2.0.0
#     - r-data.table=1.15.2
#     - r-R.utils=2.12.3
#     - r-coloc=5.2.3
#     - r-cowplot=1.1.3
#     - r-patchwork=1.2.0
RUN mkdir -p /conda-envs/cfb439855512bc80795457ea259e5135
COPY workflow/envs/fine_mapping.yml /conda-envs/cfb439855512bc80795457ea259e5135/environment.yaml

# Conda environment:
#   source: workflow/envs/locus_breaker.yml
#   prefix: /conda-envs/9e7d231b447296924baea9cb371baed0
#   name: locus_breaker
#   channels:
#     - conda-forge
#     - defaults
#     - R
#   dependencies: #R>=4.3
#     - r-base=4.3.3
#     - r-optparse=1.7.4
#     - r-tidyverse=2.0.0
#     - r-data.table=1.15.2
RUN mkdir -p /conda-envs/9e7d231b447296924baea9cb371baed0
COPY workflow/envs/locus_breaker.yml /conda-envs/9e7d231b447296924baea9cb371baed0/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/cfb439855512bc80795457ea259e5135 --file /conda-envs/cfb439855512bc80795457ea259e5135/environment.yaml && \
    mamba env create --prefix /conda-envs/9e7d231b447296924baea9cb371baed0 --file /conda-envs/9e7d231b447296924baea9cb371baed0/environment.yaml && \
    mamba clean --all -y
