rm(list = ls())
gc()
library(dplyr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
# Use detectCores() to find out how many cores are available.
num_cores <- 4
registerDoParallel(num_cores)

source(here::here("analysis/functions/seir_functions.R"))
source(here::here("analysis/functions/mcmc_functions.R"))

#########################
#### SEIR Simulation ####
#########################
time <- 365 # Number of days.
inc <- 28 # Average incubation period length.
inf <- 7 # Average infectious period length.
N <- 1000 # Population size.

target <- c(0.1, 0.25)
start <- c(53.6136261376403, 0.086607502588279)
results <- metropolis(start, target, num_sim = 5, num_iter = 10)
path <- results[[1]]
liks <- results[[2]]
best <- results[[3]]
write.table(path, file = 'path.txt', row.names = F, col.names = F)
write.table(liks, file = 'liks.txt', row.names = F, col.names = F)
write.table(best, file = 'best.txt', row.names = F, col.names = F)
