rm(list = ls())
gc()
library(dplyr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
# Use detectCores() to find out how many cores are available.
message(detectCores())
num_cores <- 24
registerDoParallel(num_cores)

time <- 365 # Number of days.
inc <- 28 # Average incubation period length.
inf <- 7 # Average infectious period length.
pop <- 1000 # Population size.

create_hh <- function() {
  # Randomly sample household sizes such that total population is 1000 
  # individuals.
  hh_size <- sample(x = c(3, 4, 5, 6), size = 340, replace = T)
  
  # Keep households such that total population is < 1000.
  hh_size <- hh_size[which(cumsum(hh_size) < pop)]
  
  leftover <- pop - sum(hh_size)
  if (leftover < 3) {
    hh <- 1:length(hh_size)
    sampled <- sample(hh[hh_size < 6], leftover)
    hh_size[sampled] <- hh_size[sampled] + 1
  } else {
    hh_size <- c(hh_size, leftover)
  }
  return(hh_size)
}

SEIR <- function(beta_H, beta_C, inc, inf, verbose = 0) {
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
  data <- data.frame(ID = 1:pop,
                    SIZE = rep(hh_size, times = hh_size),
                    HH = rep(1:length(hh_size), times = hh_size), 
                    S = c(0, rep(1, pop - 1)), 
                    E = c(1, rep(0, pop - 1)),
                    E_count = c(1, rep(0, pop - 1)), 
                    I = 0,
                    I_count = 0, 
                    R = 0, 
                    INC = c(round(rnorm(1, inc, 2)), rep(0, pop - 1)),
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
        message(t)
      }
    }
    
    # Anyone who has been infectious for as many days as their infectious
    # period is now recovered.
    recovered <- (data$INF > 0) & (data$I_count == data$INF)
    if(sum(recovered, na.rm = T) > 0) {
      data$R[recovered] <- 1
      data$I[recovered] <- 0
      data$I_count[recovered] <- 0 
    }
    
    # Anyone who has been incubating for as many days as their incubation
    # period is now infectious.
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
    }
    
    # I_H is the number of infections inside each household.
    # I_C is the number of infections outside each household.
    I_data <- data %>% group_by(HH) %>% 
      mutate(I_H = sum(I)) %>% 
      ungroup() %>% 
      mutate(I_C = sum(I) - I_H)
    
    # Calculate household risk and community risk.
    risk_H <- pmin(beta_H * data$S * I_data$I_H / pop, 1)
    risk_C <- pmin(beta_C * data$S * I_data$I_C / pop, 1)
    
    # Each individual is infected from their household or 
    # community independently with probabilities risk_H
    # and risk_C.
    new_inf_H <- rbinom(pop, 1, risk_H)
    new_inf_C <- rbinom(pop, 1, risk_C)
    
    new_exposed <- (new_inf_H == 1) | (new_inf_C == 1)
    num_new_exposed <- sum(new_exposed, na.rm = T)
    if (num_new_exposed > 0) {
      # Change status to newly exposed and add incubation period.
      data$E[new_exposed] <- 1
      random_inc <- rnorm(num_new_exposed, mean = inc, sd = 2) %>% round()
      data$INC[new_exposed] <- random_inc
      
      # Remove susceptible status.
      data$S[new_exposed] <- 0
      
      # Label infections.
      results$TYPE[new_inf_C == 1] <- 'C'
      results$TYPE[new_inf_H == 1] <- 'H'
      results$TYPE[(new_inf_H == 1) & (new_inf_C == 1)] <- 'B'
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
  
  # If incidence is 0, the proportion of household infections is undefined.
  prp <- NA
  if (idc != 0) {
    # The proportion of household infections is the proportion of infections that
    # were from the household.
    prp <- mean(results[!is.na(results$TIME), ]$TYPE == 'H')
  }
  return(c(idc, prp))
}

beta_Hs <- seq(501, 630, 1)
beta_Cs <- seq(0, 0.2, 0.01)

a <- length(beta_Hs)
b <- length(beta_Cs)

reps <- 100
idcs <- array(rep(NA, a * b * reps), dim = c(a, b, reps))
prps <- array(rep(NA, a * b * reps), dim = c(a, b, reps))
t_tot <- 0
for (i in 1:a) {
  for (j in 1:b) {
    beta_H <- beta_Hs[i]
    beta_C <- beta_Cs[j]
    
    cat(paste0(beta_H, '/630\t', 
               format(beta_C, nsmall = 2), '/0.20\t'))
      
    t_0 <- Sys.time()
    vals <- foreach (k = 1:reps, .combine = 'c') %dopar% {
      results <- SEIR(beta_H, beta_C, inc, inf, verbose = F) 
      metrics(results)
    }
      
    t_1 <- Sys.time()
    t_tot <- t_tot + (t_1 - t_0)
    cat(paste0(format(t_tot, nsmall = 2), 
               '\t(', format(t_1 - t_0, nsmall = 2), ')\t'))
    
    vals <- matrix(vals, reps, byrow = T)
    idcs[i, j, ] <- vals[, 1]
    prps[i, j, ] <- vals[, 2]
      
    cat(paste0(format(mean(vals[, 1]), nsmall = 3), '\t', 
               format(mean(vals[, 2]), nsmall = 3)))
      
    write.table(idcs, file = 'idcs_7.txt', row.names = F, col.names = F)
    write.table(prps, file = 'prps_7.txt', row.names = F, col.names = F)
    cat('\n')
  }
  message('\n')
}