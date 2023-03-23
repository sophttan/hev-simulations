#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --nodes=1

#SBATCH --mail-user=shahzar@berkeley.edu
#SBATCH --mail-type=ALL

R CMD BATCH --no-save hev.R hev.out