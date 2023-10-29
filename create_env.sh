#!/bin/bash

# Define the name of the conda environment
ENV_NAME="Salish_sea_env"

## option 1
# Create the conda environment
# conda create --name $ENV_NAME python=3.8  --file packages.txt

## option 2
# Create the conda environment with a fresh installation of Python
conda create -y -n $ENV_NAME python=3.8

# Activate the conda environment
conda activate $ENV_NAME

# Install the required packages from packages.txt
while read link; do
    conda install -y -c conda-forge $package
done < packages.txt

# Deactivate the conda environment
conda deactivate
