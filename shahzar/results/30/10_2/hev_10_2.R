rm(list = ls())
gc()
library(dplyr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
num_cores <- 24
registerDoParallel(num_cores)

#########################
#### SEIR Simulation ####
#########################
time <- 365 # Number of days.
inc <- 28 # Average incubation period length.
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

SEIR <- function(params, inc, inf, verbose = F) {
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
                    INC = c(round(rnorm(1, inc, 2)), rep(0, N - 1)),
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
      
      # Record time at which infectious period starts.
      results$TIME[new_inf] <- t
      
      # Save the number of susceptible people in each infectious individual's 
      # household.
      S_data <- data %>% group_by(HH) %>% 
        mutate(S_tot = sum(S)) %>% 
        select(HH, S_tot)
      results$S_num[new_inf == 1] <- S_data$S_tot[new_inf == 1]
    }
    
    # I_H is the number of infections inside each household.
    # I_C is the number of infections outside each household.
    I_data <- data %>% group_by(HH) %>% 
      mutate(I_H = sum(I)) %>% 
      ungroup() %>% 
      mutate(I_C = sum(I) - I_H)
    
    # Calculate household risk and community risk.
    beta_H <- params[1]
    beta_C <- params[2]
    risk_H <- beta_H * data$S * I_data$I_H / N
    risk_C <- beta_C * data$S * I_data$I_C / N
    
    # Each individual is infected from their household or community 
    # independently with probabilities risk_H and risk_C.
    new_inf_H <- rbinom(nrow(data), 1, risk_H)
    new_inf_C <- rbinom(nrow(data), 1, risk_C)
    
    new_exposed <- (new_inf_H == 1) | (new_inf_C == 1)
    num_new_exposed <- sum(new_exposed, na.rm = T)
    if (num_new_exposed > 0) {
      # Change status to newly exposed and add incubation period.
      data$E[new_exposed] <- 1
      random_inc <- rnorm(num_new_exposed, mean = inc, sd = 2) %>% round()
      data$INC[new_exposed] <- random_inc
      
      # Remove susceptible status.
      data$S[new_exposed] <- 0
      
      # Label community infections with C and household infections with H.
      results$TYPE[new_inf_C == 1] <- 'C'
      results$TYPE[new_inf_H == 1] <- 'H'
      
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
      
      results$I_num <- results$I_num + I_data$new_I_H
        
      # Label individuals with both a household and community infection with B.
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
  
  # If incidence is 0, the SAR is undefined.
  sar <- NA
  if (idc != 0) {
    # The SAR is the average SAR for each individual that was infectious.
    sar <- mean(results$I_num / results$S_num, na.rm = T)
  }
  return(c(idc, sar))
}

beta_Hs <- c(55.450)
beta_Cs <- seq(0.0850, 0.0860, 0.0001)

a <- length(beta_Hs)
b <- length(beta_Cs)

reps <- 1000
idcs <- array(rep(NA, a * b * reps), dim = c(a, b, reps))
sars <- array(rep(NA, a * b * reps), dim = c(a, b, reps))
t_tot <- 0
for (i in 1:a) {
  for (j in 1:b) {
    beta_H <- beta_Hs[i]
    beta_C <- beta_Cs[j]
    params <- c(beta_H, beta_C)
      
    t_0 <- Sys.time()
    vals <- foreach (k = 1:reps, .combine = 'c') %dopar% {
      results <- SEIR(params, inc, inf, verbose = F) 
      metrics(results)
    }
    t_1 <- Sys.time()
    t_tot <- t_tot + (t_1 - t_0)
    vals <- matrix(vals, reps, byrow = T)
    idcs[i, j, ] <- vals[, 1]
    sars[i, j, ] <- vals[, 2]
    message(paste0(format(beta_H, digits = 5, nsmall = 3), '/55.450\t', 
                   format(beta_C, digits = 5, nsmall = 4), '/0.0860\t',  
                   format(t_tot, nsmall = 2), '\t(', format(t_1 - t_0, nsmall = 2), ')\t', 
                   format(mean(vals[, 1]), nsmall = 3), '\t', 
                   format(mean(vals[, 2]), nsmall = 3)))
    write.table(idcs, file = 'idcs_10_2.txt', row.names = F, col.names = F)
    write.table(sars, file = 'sars_10_2.txt', row.names = F, col.names = F)
  }
  message('\n')
}
