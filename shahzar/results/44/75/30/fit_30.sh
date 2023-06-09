#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

#SBATCH --mail-user=shahzar@berkeley.edu
#SBATCH --mail-type=ALL

R CMD BATCH --no-save fit_30.R fit_30.out