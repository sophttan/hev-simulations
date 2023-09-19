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

#########################
#### SEIR Simulation ####
#########################
time <- 365 # Number of days.
inf <- 7 # Average infectious period length.
N <- 1000 # Population size.

create_hh <- function() {
  # Randomly sample household sizes such that total population is 1000 
  # individuals.
  hh_size <- sample(x = c(3, 4, 5, 6), size = 340, replace = T)
  
  # Keep households such that total population is < 1000.
  hh_size <- hh_size[which(cumsum(hh_size) < N)]
  
  leftover <- N - sum(hh_size)
  if (leftover < 3) {
    hh <- 1:length(hh_size)
    sampled <- sample(hh[hh_size < 6], leftover)
    hh_size[sampled] <- hh_size[sampled] + 1
  } else {
    hh_size <- c(hh_size, leftover)
  }
  return(hh_size)
}

SEIR_env <- function(params, inf = 7, verbose = F) {
  hh_size <- create_hh()
  
  # Create frame for running the simulation.
  # ID: ID of individual.
  # SIZE: size of individual's household.
  # HH: ID of individual's household.
  # S: susceptibility status.
  # E: exposed status.
  # E_count: number of days since exposed.
  # I: infectious status.
  # I_count: number of days since infectious.
  # R: recovered status.
  # INC: incubation period.
  # INF: infectious period.
  data <- data.frame(ID = 1:N,
                     SIZE = rep(hh_size, times = hh_size),
                     HH = rep(1:length(hh_size), times = hh_size), 
                     S = rep(1, N), 
                     E = rep(0, N),
                     E_count = rep(0, N), 
                     I = 0,
                     I_count = 0, 
                     R = 0, 
                     INC = rep(0, N),
                     INF = 0)
  
  # Create frame for storing results.
  # ID: ID of individual.
  # SIZE: size of individual's household.
  # HH: ID of individual's household.
  # TYPE: the kind of infection: household (H), community (C), or both (B).
  # TIME: when the individual became infectious.
  # S_NUM: number of susceptible people in individual's household when their 
  #        infectious period begins.
  # I_NUM: number of people in household that this individual infected over 
  #        their infectious period.
  results <- data[, 1:3] %>% mutate(TYPE = NA, TIME = NA, S_NUM = NA, I_NUM = 0)
    
  for(t in 1:time) {
    if (verbose) {
      if (t %% 10 == 0) {
        cat(paste0(t, ' '))
      }
    }
    
    # Anyone who has been infectious for as many days as their infectious period
    # is now recovered.
    recovered <- (data$INF > 0) & (data$I_count == data$INF)
    if(sum(recovered, na.rm = T) > 0) {
      data$R[recovered] <- 1
      data$I[recovered] <- 0
      data$I_count[recovered] <- 0 
    }
    
    # Anyone who has been incubating for as many days as their incubation period
    # is now infectious.
    new_inf <- (data$INC > 0) & (data$E_count == data$INC)
    num_new_inf <- sum(new_inf, na.rm = T)
    if(num_new_inf > 0) {
      # Change status to newly infectious and add infectious period.
      data$I[new_inf] <- 1
      random_inf <- rnorm(num_new_inf, mean = inf, sd = 1) %>% round()
      data$INF[new_inf] <- random_inf
      
      # Remove exposure status and exposure count.
      data$E[new_inf] <- 0
      data$E_count[new_inf] <- 0 
      
      results$TIME[new_inf] <- t
    }
    
    beta_E <- params[1]
    risk_E <- pmin(beta_E * data$S, 1)
    
    new_exposed <- rbinom(N, 1, risk_E)
    num_new_exposed <- sum(new_exposed, na.rm = T)
    if(num_new_exposed > 0) {
      # Change status to newly exposed and add incubation period.
      data$E[new_exposed == 1] <- 1
      random_inc <- rlnorm(num_new_exposed, meanlog = log(29.8), sdlog = 0.45) %>% round()
      data$INC[new_exposed == 1] <- random_inc
      
      # Remove susceptible status.
      data$S[new_exposed == 1] <- 0
        
      # Label infection types.
      results$TYPE[new_exposed == 1] <- 'E'
    }
    
    # Increment exposure and infectious counters.
    data$E_count[data$E == 1] <- data$E_count[data$E == 1] + 1
    data$I_count[data$I == 1] <- data$I_count[data$I == 1] + 1
  }
  return(results)
}

metrics <- function(results) {
  # Incidence is the proportion of the population that became infected.
  idc <- mean(!is.na(results$TIME))
  return(c(idc))
}

# 5% Cumulative Incidence
beta_E <- 0.00015
params <- c(beta_E)

n_sims <- 2500
i <- 1001
j <- 0
while (i <= n_sims) {
    results <- SEIR_env(params)
    j <- j + 1
    cat(paste0(j, '\t', i, '\n'))
    if (metrics(results) > 0.04 & metrics(results) < 0.06) {
        write.csv(results, file = paste0('5/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}
cat(paste0('05%:\t', j, '\n'))

# 10% Cumulative Incidence
beta_E <- 0.00031
params <- c(beta_E)

n_sims <- 2500
i <- 1001
j <- 0
while (i <= n_sims) {
    results <- SEIR_env(params)
    j <- j + 1
    cat(paste0(j, '\t', i, '\n'))
    if (metrics(results) > 0.09 & metrics(results) < 0.11) {
        write.csv(results, file = paste0('10/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}
cat(paste0('10%:\t', j, '\n'))

# 30% Cumulative Incidence
beta_E <- 0.00107
params <- c(beta_E)

n_sims <- 2500
i <- 1001
j <- 0
while (i <= n_sims) {
    results <- SEIR_env(params)
    j <- j + 1
    cat(paste0(j, '\t', i, '\n'))
    if (metrics(results) > 0.29 & metrics(results) < 0.31) {
        write.csv(results, file = paste0('30/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}
cat(paste0('30%:\t', j, '\n'))
