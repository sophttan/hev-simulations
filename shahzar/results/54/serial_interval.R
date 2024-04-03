rm(list = ls())
gc()
library(dplyr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- 24
registerDoParallel(num_cores)

# Number of simulations
n <- 1000000000

# Compute latent (incubation) period, infectious period, and serial interval.
l <- rlnorm(n, meanlog = log(29.8), sdlog = 0.45) %>% round()
i <- rnorm(n, mean = 7, sd = 1) %>% round()
si <- (l + i/2) %>% round()

# Compute empirical distribution of serial interval
si_dist <- table(factor(si, levels = 0:max(si))) / n

saveRDS(si, file = 'si.rds')
saveRDS(si_dist, file = 'si_dist.rds')
write.table(si, file = 'si.txt', row.names = F, col.names = F)
write.table(as.numeric(si_dist), file = 'si_dist.txt', row.names = F, col.names = F)
