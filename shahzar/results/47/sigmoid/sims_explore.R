rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- detectCores()
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

SEIR <- function(params, inf = 7, verbose = F) {
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
                     S = c(0, rep(1, N - 1)), 
                     E = c(1, rep(0, N - 1)),
                     E_count = c(1, rep(0, N - 1)), 
                     I = 0,
                     I_count = 0, 
                     R = 0, 
                     INC = c(round(rlnorm(1, meanlog = log(29.8), sdlog = 0.45)), rep(0, N - 1)),
                     INF = 0)
  
  # Create frame for storing results.
  # ID: ID of individual.
  # SIZE: size of individual's household.
  # HH: ID of individual's household.
  # TYPE: the kind of infection: household (H), community (C), or both (B).
  # TIME: when the individual became infectious.
  # S_num: number of susceptible people in individual's household when their 
  #        infectious period begins.
  # I_num: number of people in household that this individual infected over 
  #        their infectious period.
  results <- data[, 1:3] %>% mutate(TYPE = NA, TIME = NA, S_num = NA, I_num = 0)
  results$TYPE[1] <- '0'
    
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
    
    # Risk is growing logistically, with rate determined by params[2]
    beta_E <- params[1]
    gamma <- params[2]
    risk_E <- pmin(beta_E * data$S / (1 + exp(-gamma * (t - 100))), 1)
    
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


### GAMMA = 0.1 ### 
params <- c(0.0001, 0.1)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    if (metrics(results) > 10/1000) {
        write.csv(results, file = paste0('explore/1/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}

params <- c(0.0010, 0.1)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    if (metrics(results) > 10/1000) {
        write.csv(results, file = paste0('explore/2/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}

params <- c(0.0100, 0.1)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    if (metrics(results) > 10/1000) {
        write.csv(results, file = paste0('explore/3/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}

### GAMMA = 0.01

params <- c(0.0001, 0.01)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    if (metrics(results) > 10/1000) {
        write.csv(results, file = paste0('explore/4/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}

params <- c(0.0010, 0.01)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    if (metrics(results) > 10/1000) {
        write.csv(results, file = paste0('explore/5/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}

params <- c(0.0100, 0.01)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    if (metrics(results) > 10/1000) {
        write.csv(results, file = paste0('explore/6/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}

### GAMMA = 0.001

params <- c(0.0001, 0.001)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    write.csv(results, file = paste0('explore/7/', i, '.csv'), row.names = F)
    i <- i + 1
}

params <- c(0.0010, 0.001)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    write.csv(results, file = paste0('explore/8/', i, '.csv'), row.names = F)
    i <- i + 1
}

params <- c(0.0100, 0.001)
n_sims <- 1000
i <- 0
while (i < n_sims) {
    results <- SEIR(params, inf)
    write.csv(results, file = paste0('explore/9/', i, '.csv'), row.names = F)
    i <- i + 1
}