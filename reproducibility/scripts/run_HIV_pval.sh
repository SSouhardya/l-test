#!/bin/bash
#SBATCH -J HIV_ltest  # A single job name for the array
#SBATCH -p # Partition
#SBATCH -c 1 # number of cores
#SBATCH -t 2-00:00  # Running time in the format - D-HH:MM
#SBATCH --mem 4000 # Memory request - 1000 corresponds to 1GB
#SBATCH -o out_%a.out # Standard output
#SBATCH -e err_%a.err # Standard error
Rscript ~/reproducibility/HIV_data_pval.R

