set.seed(1234)

library(foreach)
library(doParallel)
library(tidyverse)
library(here)
library(latex2exp)

# choosing the number of cores to do parallelisation on
numCores = detectCores()-2
registerDoParallel(numCores)

source(here::here("analysis/functions/seir_functions.R"))

#######################
# Simulating the data #
#######################

time <- 365 # Number of days.
inf <- 7 # Average infectious period length.
N <- 1000 # Population size.
sims <- 500 # Number of simulations

# start = vector of values with (beta_h, beta_c, beta_e)
# output: sims many simulations of 1 parameter set
sim_blend = function(start, sims, inf) {
  incidence<-NULL
  sar<-NULL
  results_data<-NULL
  
  for (i in 1:sims) {
    results <- SEIR_blend(start, inf)
    res <- metrics_blend(results)
    incidence<-c(incidence,res[1])
    sar<-c(sar,res[2])
    
    results <- results %>% filter(!is.na(TIME)) %>% mutate(i=i)
    results_data <- rbind(results_data, results)
  }
  return(results_data)
}

# inputting the calibrated values for each incidence and then concatenating into one big matrix
# cols: beta_h, beta_c, beta_e and rows: 25%, 50%, 75%
param_5 = c(10, 0.005, 0.0001, 28, 0.005, 0.0001, 42, 0.03, 0.00005)
param_5 = matrix(param_5, nrow = 3, ncol = 3, byrow = TRUE)

param_10 = c(12, 0, 0.00025, 27, 0.005, 0.00015, 43, 0.035, 0.00005)
param_10 = matrix(param_10, nrow = 3, ncol = 3, byrow = TRUE)

param_30 = c(12, 0.005, 0.0008, 29, 0.01, 0.00055, 41, 0.04, 0.00025)
param_30 = matrix(param_30, nrow = 3, ncol = 3, byrow = TRUE)

param = rbind(param_5, param_10, param_30)

# getting the simulations
blended = foreach(j=1:9, .combine = rbind, .export = ls(envir = globalenv()), .packages = c("here", "tidyverse")) %dopar% {
  sim_blend(param[j,], sims, inf) %>% mutate(param = j)
}

write.csv(blended, file = here("nila/hev/estimate-prop/blend-sims/blended.csv"))