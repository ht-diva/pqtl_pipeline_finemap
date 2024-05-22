FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="fdabd2690178ff7860dfef3b02463c7607a55214fea9c8f71ce79c3cd9abb87b"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/r_environment.yml
#   prefix: /conda-envs/9d0481ef3312907f0c0b2a7736a65351
#   name: r_finemap
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#     - R
#   dependencies: #R>=4.3
#     - plink2=2.00a5.10-0
#     - gcta=1.94.1-0
#     - r-base=4.3.3
#     - r-optparse=1.7.4
#     - r-tidyverse=2.0.0
#     - r-R.utils=2.12.3
#     - r-coloc=5.2.3
#     - r-cowplot=1.1.3
#     - r-corrplot=0.92
#     - r-bigsnpr=1.12.2
#     - r-stringi=1.8.3
#     - r-patchwork=1.2.0
#     - r-plyr=1.8.9
#     - r-reshape2=1.4.4
#     - r-RColorBrewer=1.1-3
#     - r-igraph=2.0.2
#     - r-Matrix=1.6-5
#     - bioconductor-Gviz=1.46.1
#     - bioconductor-EnsDb.Hsapiens.v75=2.99.0
RUN mkdir -p /conda-envs/9d0481ef3312907f0c0b2a7736a65351
COPY workflow/envs/r_environment.yml /conda-envs/9d0481ef3312907f0c0b2a7736a65351/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/9d0481ef3312907f0c0b2a7736a65351 --file /conda-envs/9d0481ef3312907f0c0b2a7736a65351/environment.yaml && \
    mamba clean --all -y
