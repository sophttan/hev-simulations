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

###############
## Panmictic ##
###############

dir <- 'pan/'

# 5% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 10% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}

            
# 30% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}



########
## 0% ##
########

dir <- '0/'

# 5% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 10% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 30% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}



#########
## 25% ##
#########

dir <- '25/'

# 5% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 10% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 30% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}



#########
## 50% ##
#########

dir <- '50/'

# 5% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 10% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 30% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


   
#########
## 75% ##
#########

dir <- '75/'

# 5% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 10% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 30% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


  
##########
## 100% ##
##########

dir <- '100/'

# 5% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 10% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}


# 30% Incidence
n_sims <- 2500
obs <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  results <- read.csv(f)
  
  detected <- rbinom(1000, 1, 0.5)
  detect_time <- round(runif(1000, -0.5, 7.5))
  observation_time <- ifelse(detected, detect_time, NA)
  results$TIME <- results$TIME + observation_time
  
  write.csv(results, file = f, row.names = F)
}