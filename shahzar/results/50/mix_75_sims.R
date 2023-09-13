rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- 24
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

SEIR_mix <- function(params, inf = 7, verbose = F) {
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
                     INC = c(rlnorm(1, meanlog = log(29.8), sdlog = 0.45), rep(0, N - 1)),
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
      
      # Record time at which infectious period starts.
      results$TIME[new_inf] <- t
      
      # Save the number of susceptible people in each infectious individual's 
      # household.
      S_data <- data %>% group_by(HH) %>% 
        mutate(S_TOT = sum(S)) %>% 
        select(HH, S_TOT)
      results$S_NUM[new_inf == 1] <- S_data$S_TOT[new_inf == 1]
    }
    
    # I_H is the number of infections inside each household.
    # I_C is the number of infections outside each household.
    I_data <- data %>% group_by(HH) %>% 
      mutate(I_H = sum(I)) %>% 
      ungroup() %>% 
      mutate(I_C = sum(I) - I_H)
    
    beta_H <- params[1]
    beta_C <- params[2]
    beta_E <- params[3]
    risk_H <- pmin(beta_H * data$S * I_data$I_H / N, 1)
    risk_C <- pmin(beta_C * data$S * I_data$I_C / N, 1)
    risk_E <- pmin(beta_E * data$S, 1)

    new_inf_H <- rbinom(N, 1, risk_H)
    new_inf_C <- rbinom(N, 1, risk_C)
    new_inf_E <- rbinom(N, 1, risk_E)
    new_exposed <- (new_inf_H == 1) | (new_inf_C == 1) | (new_inf_E == 1)
    num_new_exposed <- sum(new_exposed, na.rm = T)
    if (num_new_exposed > 0) {
      # Change status to newly exposed and add incubation period.
      data$E[new_exposed] <- 1
      random_inc <- rlnorm(num_new_exposed, meanlog = log(29.8), sdlog = 0.45) %>% round()
      data$INC[new_exposed] <- random_inc
      
      # Remove susceptible status.
      data$S[new_exposed] <- 0
      
      # Get number of new infections in each household.
      I_data <- I_data %>%
        select(ID, HH, I, I_H) %>%
        mutate(new_I_H = new_inf_H) %>%
        group_by(HH) %>%
        # Find households with at least 1 currently infectious individual. If 
        # exactly 1 infectious individual in household, assign all new H 
        # exposures to that individual. If there are multiple infectious 
        # individuals, assign all infections to the infectious individual with 
        # the first ID.
        mutate(new_I_H = ifelse(I == 1 & ID == first(ID[I == 1]), 
                                sum(new_I_H), 0))
      
      results$I_NUM <- results$I_NUM + I_data$new_I_H
      
      # Label infections types.
      results$TYPE[new_inf_H == 1] <- 'H'
      results$TYPE[new_inf_C == 1] <- 'C'
      results$TYPE[new_inf_E == 1] <- 'E'
      results$TYPE[(new_inf_H == 1) & (new_inf_C == 1)] <- 'HC'
      results$TYPE[(new_inf_H == 1) & (new_inf_E == 1)] <- 'HE'
      results$TYPE[(new_inf_C == 1) & (new_inf_E == 1)] <- 'CE'
      results$TYPE[(new_inf_H == 1) & (new_inf_C == 1) & (new_inf_E == 1)] <- 'HCE'
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
  
  cases <- results[!is.na(results$TIME), ]
  prp <- mean((cases$TYPE == 'H') | (cases$TYPE == 'HC') | (cases$TYPE == 'HCE'))
    
  return(c(idc, prp))
}

#########
## 75% ##
#########

# 5% Cumulative Incidence
beta_H <- 42
beta_C <- 0.030
beta_E <- 0.00005
params <- c(beta_H, beta_C, beta_E)

n_sims <- 1000
i <- 1
j <- 0
while (i <= n_sims) {
    results <- SEIR_mix(params)
    j <- j + 1
    if (all(metrics(results) > c(0.04, 0.70)) & all(metrics(results) < c(0.06, 0.80))) {
        write.csv(results, file = paste0('75/5/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}
cat(paste0('05%:\t', j, '\n'))

# 10% Cumulative Incidence
beta_H <- 43
beta_C <- 0.035
beta_E <- 0.00005
params <- c(beta_H, beta_C, beta_E)

n_sims <- 1000
i <- 1
j <- 0
while (i <= n_sims) {
    results <- SEIR_mix(params)
    j <- j + 1
    if (all(metrics(results) > c(0.09, 0.70)) & all(metrics(results) < c(0.11, 0.80))) {
        write.csv(results, file = paste0('75/10/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}
cat(paste0('10%:\t', j, '\n'))

# 30% Cumulative Incidence
beta_H <- 41
beta_C <- 0.040
beta_E <- 0.00025
params <- c(beta_H, beta_C, beta_E)

n_sims <- 1000
i <- 1
j <- 0
while (i <= n_sims) {
    results <- SEIR_mix(params)
    j <- j + 1
    if (all(metrics(results) > c(0.29, 0.70)) & all(metrics(results) < c(0.31, 0.80))) {
        write.csv(results, file = paste0('75/30/', i, '.csv'), row.names = F)
        i <- i + 1
    }
}
cat(paste0('30%:\t', j, '\n'))