#!/bin/bash

# stop on errors
set -e

echo "Running snakemake..."

# Remove tmp/ dir and slurm files
snakemake clean --cores 1

# Make fresh tmp/ dir 
mkdir -p tmp

# Run the main analysis on slurm cluster
snakemake \
    --use-conda \
    --conda-prefix env \
    -j 999 \
    --cluster-config cluster.yml \
    --cluster "sbatch -p {cluster.partition} -c {cluster.cpus} -t {cluster.time} -J {cluster.name} -o ./tmp/slurm-%x.%j.out" \
    --latency-wait 60

echo "Run of snakemake complete."
