#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1

#SBATCH --mail-user=shahzar@berkeley.edu
#SBATCH --mail-type=ALL

R CMD BATCH --no-save p2p_05_sims.R p2p_05_sims.out