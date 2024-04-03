rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)
library(EpiEstim)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- 4
registerDoParallel(num_cores)

file_dir <- '../50/pan/'
results <- read.csv(paste0(file_dir, '30/1.csv'))


## EpiEstim will define imported cases as those without an obvious infector
## Serial interval: latent period + 1/2 infectious
# incubation (exposure to symptom onset)
# latent period (exposure to infectious)
# infectious period (infectious to cessation of infectiousness)

# R_HH and R_HH/R plots over time
# EpiEstim calculate R with imported cases
# By Friday night
###############
## Panmictic ##
###############

save_dir <- 'pan/'
file_dir <- '../50/pan/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_30.txt'), row.names = F, col.names = F)




##############
### Pure E ###
##############

save_dir <- '0/'
file_dir <- '../50/0/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_30.txt'), row.names = F, col.names = F)




###############
##  25% P2P  ##
###############

save_dir <- '25/'
file_dir <- '../50/25/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_30.txt'), row.names = F, col.names = F)



###############
##  50% P2P  ##
###############

save_dir <- '50/'
file_dir <- '../50/50/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_30.txt'), row.names = F, col.names = F)




###############
##  75% P2P  ##
###############

save_dir <- '75/'
file_dir <- '../50/75/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_30.txt'), row.names = F, col.names = F)




##############
## 100% P2P ##
##############

save_dir <- '100/'
file_dir <- '../50/100/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_30.txt'), row.names = F, col.names = F)

# Examine this and the shrinkage effect over time
# Look at earlier time periods.
#     Surmountable if it works for earlier times but not later ones
# Look at other formulations which can distinguish between E vs HH (R_HH alone?)
# Normal reproduction number method and look at cases which are "imports" (no obvious infectors)
#    Some cases won't have a plausible infector. Label those as environmental.
#    R package which computes this and has some accounting for imports
#    How would environmental look relative to person to person for pure reproduction number methods?
# With detailed spatial data, could make likelihood network based on contacts (closer vs farther)
#    Meta-population of 
# x-axis with proportion of household cases instead of person-to-person transmission
# Big issue: relation between household and person-to-person is not stable over different incidences.
#    # As long as it's household, there'll be that problem.

# Three Routes
#   1. Not a problem early in outbreaks. Later on, community cases overrun household cases. [Check R_HH/R early on]
#   2. Imported case approach works. [Quick check]
#   3. New formulation or adjustment to address this. E.g., just R_HH, or other formulation. [Check other formulations]