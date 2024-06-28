#!/bin/bash
#SBATCH -J  l_test # Add a job name for the array
#SBATCH -p # Add a partition name
#SBATCH -c 1 # number of cores
#SBATCH -t 1-00:00  # Running time in the format - D-HH:MM
#SBATCH --mem 4000 # Memory request - 1000 corresponds to 1GB
#SBATCH -o out_%a.out # Standard output
#SBATCH -e err_%a.err # Standard error
Rscript ~/reproducibility/single_test.R # Add the file you want to run here
