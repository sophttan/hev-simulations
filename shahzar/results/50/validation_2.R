rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- 32
registerDoParallel(num_cores)

get_idc <- function(results) {
  # Incidence is the proportion of the population that became infected.
  idc <- mean(!is.na(results$TIME))
  return(idc)
}

get_sar <- function(results) {
  # The SAR is the average SAR for each individual that was infectious.
  sar <- mean(results$I_num / results$S_num, na.rm = T)
  return(sar)
}

get_prp <- function(results) {
  # The proportion of cases caused by household infections.
  cases <- results[!is.na(results$TIME), ]
  prp <- mean((cases$TYPE != 'E') & (cases$TYPE != '0'))
  return(prp)
}


  
##########
## 100% ##
##########

dir <- '100/'

# 5% Incidence
n_sims <- 100000
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  if (file.exists(f)) {
     results <- read.csv(f)
     c(get_idc(results), get_sar(results)) 
  } else {
      -1
  }
}
obs <- matrix(obs[obs != -1], ncol = dim(obs)[2])
cat(length(obs), '\n')
saveRDS(obs, file = paste0(dir, 'obs_05.rds'))
write.table(obs, file = paste0(dir, 'obs_05.txt'), row.names = F, col.names = F)


# 10% Incidence
n_sims <- 100000
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  if (file.exists(f)) {
     results <- read.csv(f)
     c(get_idc(results), get_sar(results)) 
  } else {
      -1
  }
}
obs <- matrix(obs[obs != -1], ncol = dim(obs)[2])
cat(length(obs), '\n')
saveRDS(obs, file = paste0(dir, 'obs_10.rds'))
write.table(obs, file = paste0(dir, 'obs_10.txt'), row.names = F, col.names = F)


# 30% Incidence
n_sims <- 100000
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  if (file.exists(f)) {
     results <- read.csv(f)
     c(get_idc(results), get_sar(results))
  } else {
      -1
  }
}
obs <- matrix(obs[obs != -1], ncol = dim(obs)[2])
cat(length(obs), '\n')
saveRDS(obs, file = paste0(dir, 'obs_30.rds'))
write.table(obs, file = paste0(dir, 'obs_30.txt'), row.names = F, col.names = F)