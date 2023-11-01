#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1

#SBATCH --mail-user=shahzar@berkeley.edu
#SBATCH --mail-type=ALL

R CMD BATCH --no-save p2p_sims_30.R p2p_sims_30.out