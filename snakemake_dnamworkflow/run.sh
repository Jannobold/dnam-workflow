#!/usr/bin/env bash
#SBATCH --cluster=omics
#SBATCH --partition=shortterm
#SBATCH --nodes=1
#SBATCH --requeue

export PATH="./software:$PATH"

echo "RUNNING SNAKEMAKE WORKFLOW..."
source ~/.bashrc
# for larger datasets, consider to adjust --mem and --time accordingly
snakemake -j 22 all --cluster "sbatch --mem 10G --time 0-01:00 --partition=shortterm" --rerun-incomplete

#EXAMPLES TO RUN WORKFLOWS BELOW
# snakemake -j 22 all
# snakemake -j 22 --rerun-incomplete

#WHEN JOBS WERE CANCELLED UNEXPECTEDLY (unlock directories & files)
# snakemake --unlock