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

SEIR <- function(params, inf, verbose = F) {
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
      
      results$I_num <- results$I_num + I_data$new_I_H
      
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

prp_true <- function(results) {
  cases <- results[!is.na(results$TIME), ]
  cases_HH = (cases$TYPE == 'H') | (cases$TYPE == 'HC') | (cases$TYPE == 'HCE')
  prp <- mean(cases_HH)
  return(prp)
}

method_HH <- function(results) {
  f <- results %>% filter(!is.na(TIME))
  f <- f %>% group_by(HH) %>% 
    mutate(day_limits = list(TIME)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(IS_H = any((TIME - unlist(day_limits)) < 45 & (TIME - unlist(day_limits)) > 7))
  
  return(mean(f$IS_H))
}

probability <- function(cases, index, rel_p_hh=1) {
  inc_prim <- cases$TIME
  inc_sec <- cases$TIME[index]
  hh_prim <- cases$HH
  hh_sec <- cases$HH[index]
  # calculate relative likelihood that primary case caused secondary case based on their incidences
  # serial interval is approximate
  results <- (rel_p_hh*(hh_prim==hh_sec)+(hh_prim!=hh_sec))*dnorm(inc_sec-inc_prim, mean=31.5, sd=4)
  # if primary case happens after secondary case, set probability to 0
  results[inc_prim >= inc_sec] <- 0
  if(sum(results) == 0){
    return(results)
  }
  
  return(results/sum(results))
}

method_R <- function(results) {
  cases <- results %>% 
    filter(!is.na(TIME)) %>%
    mutate(ID = 1:n())
  
  R <- rep(0, nrow(cases))
  R_HH <- rep(0, nrow(cases))
  
  for (k in 1:nrow(cases)) {
    probs <- probability(cases, k, 1)
    HH_probs <- probs
    HH_probs[cases$HH != cases$HH[k]] <- 0
    R <- R + probs
    R_HH <- R_HH + HH_probs
  }
  return(mean(R_HH) / mean(R))
}

values <- function(results) {
  return(c(prp_true(results), method_HH(results), method_R(results)))
}

# 5% Cumulative Incidence
beta_H <- 42
beta_C <- 0.030
beta_E <- 0.00005
params <- c(beta_H, beta_C, beta_E)

n_sims <- 1000
vals <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- SEIR(params, inf)
  write.csv(results, file = paste0('5/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims, byrow = T)
saveRDS(vals, file = '5/vals.rds')
write.table(vals, file = '5/vals.txt', row.names = F, col.names = F)


# 10% Cumulative Incidence
beta_H <- 43
beta_C <- 0.035
beta_E <- 0.00005
params <- c(beta_H, beta_C, beta_E)

n_sims <- 1000
vals <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- SEIR(params, inf)
  write.csv(results, file = paste0('10/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims, byrow = T)
saveRDS(vals, file = '10/vals.rds')
write.table(vals, file = '10/vals.txt', row.names = F, col.names = F)


# 30% Cumulative Incidence
beta_H <- 41
beta_C <- 0.040
beta_E <- 0.00025
params <- c(beta_H, beta_C, beta_E)

n_sims <- 1000
vals <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- SEIR(params, inf)
  write.csv(results, file = paste0('30/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims, byrow = T)
saveRDS(vals, file = '30/vals.rds')
write.table(vals, file = '30/vals.txt', row.names = F, col.names = F)
