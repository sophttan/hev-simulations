#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1

#SBATCH --mail-user=shahzar@berkeley.edu
#SBATCH --mail-type=ALL

R CMD BATCH --no-save mix_50_sims.R mix_50_sims.out