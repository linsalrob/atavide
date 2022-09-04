#!/bin/bash
#SBATCH --job-name=Atavide
#SBATCH --time=10-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -o atavide-%j.out
#SBATCH -e atavide-%j.err

snakemake --configfile atavide.yaml -s  ~/GitHubs/atavide/workflow/atavide.snakefile --profile slurm --restart-times 5 --rerun-incomplete --reason

