rm(list = ls())
gc()
library(dplyr)
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
  
  # The SAR is the average SAR for each individual that was infectious.
  sar <- mean(results$I_NUM / results$S_NUM, na.rm = T)
  
  # The proportion of cases caused by household infections.
  cases <- results[!is.na(results$TIME), ]
  prp <- mean((cases$TYPE == 'H') | (cases$TYPE == 'HC') | (cases$TYPE == 'HCE'))

  return(c(idc, sar, prp))
}

beta_Hs <- seq(10, 13, 1) # 4
beta_Cs <- seq(0, 0.005, 0.001) # 6
beta_Es <- seq(0.00015, 0.00035, 0.00001) # 21

d_H <- length(beta_Hs)
d_C <- length(beta_Cs)
d_E <- length(beta_Es)

reps <- 200
idcs <- array(rep(NA, d_H * d_C * d_E * reps), dim = c(d_H, d_C, d_E, reps))
sars <- array(rep(NA, d_H * d_C * d_E * reps), dim = c(d_H, d_C, d_E, reps))
prps <- array(rep(NA, d_H * d_C * d_E * reps), dim = c(d_H, d_C, d_E, reps))
t_tot <- 0
for (i in 1:d_H) {
  for (j in 1:d_C) {
    for (k in 1:d_E) {
      beta_H <- beta_Hs[i]
      beta_C <- beta_Cs[j]
      beta_E <- beta_Es[k]
      
      params <- c(beta_H, beta_C, beta_E)
      
      cat(paste0(format(beta_H, nsmall = 0), '/13\t', 
                 format(beta_C, nsmall = 3), '/0.005\t',
                 format(beta_E, nsmall = 5, scientific = F), '/0.00035\t'))
        
      t_0 <- Sys.time()
      vals <- foreach (l = 1:reps, .combine = 'c') %dopar% {
        results <- SEIR_mix(params) 
        metrics(results)
      }
      t_1 <- Sys.time()
      t_tot <- t_tot + (t_1 - t_0)
        
      cat(paste0(format(t_tot, nsmall = 2, digits = 4), '\t(', 
                 format(t_1 - t_0, nsmall = 2, digits = 4), ')\t'))
      
      vals <- matrix(vals, reps, byrow = T)
      idcs[i, j, k, ] <- vals[, 1]
      sars[i, j, k, ] <- vals[, 2]
      prps[i, j, k, ] <- vals[, 3]
        
      cat(paste0(format(round(mean(vals[, 1]), 3), nsmall = 3), '\t', 
                 format(round(mean(vals[, 2]), 3), nsmall = 3), '\t', 
                 format(round(mean(vals[, 3]), 3), nsmall = 3), '\n'))
    }
    cat('\n')
    saveRDS(idcs, file = 'idcs.rds')
    saveRDS(sars, file = 'sars.rds')
    saveRDS(prps, file = 'prps.rds')
    
    write.table(idcs, file = 'idcs.txt', row.names = F, col.names = F)
    write.table(sars, file = 'sars.txt', row.names = F, col.names = F)
    write.table(prps, file = 'prps.txt', row.names = F, col.names = F)
  }
  message('\n')
}
