rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- 64
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
  # E_COUNT: number of days since exposed.
  # I: infectious status.
  # I_COUNT: number of days since infectious.
  # R: recovered status.
  # INC: incubation period.
  # INF: infectious period.
  data <- data.frame(ID = 1:N,
                     SIZE = rep(hh_size, times = hh_size),
                     HH = rep(1:length(hh_size), times = hh_size), 
                     S = c(0, rep(1, N - 1)), 
                     E = c(1, rep(0, N - 1)),
                     E_COUNT = c(1, rep(0, N - 1)), 
                     I = 0,
                     I_COUNT = 0, 
                     R = 0, 
                     INC = c(round(rlnorm(1, meanlog = log(29.8), sdlog = 0.45)), rep(0, N - 1)),
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
    recovered <- (data$INF > 0) & (data$I_COUNT == data$INF)
    if(sum(recovered, na.rm = T) > 0) {
      data$R[recovered] <- 1
      data$I[recovered] <- 0
      data$I_COUNT[recovered] <- 0 
    }
    
    # Anyone who has been incubating for as many days as their incubation period
    # is now infectious.
    new_inf <- (data$INC > 0) & (data$E_COUNT == data$INC)
    num_new_inf <- sum(new_inf, na.rm = T)
    if(num_new_inf > 0) {
      # Change status to newly infectious and add infectious period.
      data$I[new_inf] <- 1
      random_inf <- rnorm(num_new_inf, mean = inf, sd = 1) %>% round()
      data$INF[new_inf] <- random_inf
      
      # Remove exposure status and exposure count.
      data$E[new_inf] <- 0
      data$E_COUNT[new_inf] <- 0 
      
      # Record time at which infectious period starts.
      results$TIME[new_inf] <- t
      
      # Save the number of susceptible people in each infectious individual's 
      # household.
      S_data <- data %>% group_by(HH) %>% 
        mutate(S_tot = sum(S)) %>% 
        select(HH, S_tot)
      results$S_NUM[new_inf == 1] <- S_data$S_tot[new_inf == 1]
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
    risk_H <- pmin(beta_H * data$S * I_data$I_H / N, 1)
    risk_C <- pmin(beta_C * data$S * I_data$I_C / N, 1)
    
    # Each individual is infected from their household or community 
    # independently with probabilities risk_H and risk_C.
    new_inf_H <- rbinom(N, 1, risk_H)
    new_inf_C <- rbinom(N, 1, risk_C)
    
    new_exposed <- (new_inf_H == 1) | (new_inf_C == 1)
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
      
      # Label infection types.
      results$TYPE[new_inf_C == 1] <- 'C'
      results$TYPE[new_inf_H == 1] <- 'H'
      results$TYPE[(new_inf_H == 1) & (new_inf_C == 1)] <- 'HC'
    }
    
    # Increment exposure and infectious counters.
    data$E_COUNT[data$E == 1] <- data$E_COUNT[data$E == 1] + 1
    data$I_COUNT[data$I == 1] <- data$I_COUNT[data$I == 1] + 1
  }
  return(results)
}

metrics <- function(results) {
  # Incidence is the proportion of the population that became infected.
  idc <- mean(!is.na(results$TIME))
  
  # If incidence is 0, the SAR is undefined.
  sar <- NA
  prp <- NA
  if (idc != 0) {
    # The SAR is the average SAR for each individual that was infectious.
    sar <- mean(results$I_NUM / results$S_NUM, na.rm = T)
    cases <- results[!is.na(results$TIME), ]
    prp <- mean((cases$TYPE == 'H') | (cases$TYPE == 'HC') | (cases$TYPE == 'HCE'))
  }
  return(c(idc, sar, prp))
}

##############################
#### Metropolis Algorithm ####
##############################
score <- function(fit, tgt) {
  # The score is the L² distance of the observed values from the target.
  return(sum((fit - tgt)^2))
}

# The likelihood is calculated by first averaging the incidence and SAR over n
# simulations with the state parameters. The likelihood is the negative log
# score of the average incidence and SAR.
likelihood <- function(state, tgt, n = 300) {
  # If any parameter is negative, do not transition to that state.
  if (any(state < 0)) {
    return(-Inf)
  }
  # Otherwise, find the average incidence and SAR and compute likelihood.
  vals <- foreach (i = 1:n, .combine = c) %dopar% {
    results <- SEIR(state)
    metrics(results)
  }
  vals <- matrix(vals, n, byrow = T)
  fit <- colMeans(vals)
  return(-log(score(fit, tgt)))
}

# Proposal function
q <- function(state, sds = c(0.01, 0.001)) {
  # Sample from a multivariate normal distributions centered at the current 
  # state. The SDs roughly correspond to the step-size of the chain for each 
  # parameter.
  return(rnorm(n = 2, mean = state, sd = sds))
}

# MCMC
metropolis <- function(start, tgt, num_sim, num_iter) {
  path <- matrix(NA, num_iter + 1, 2)
  liks <- rep(NA, num_iter + 1)
  
  # Initialize current state.
  curr <- start
  curr_lik <- likelihood(curr, tgt, num_sim)
  
  # Initialize best state.
  best <- curr
  best_lik <- curr_lik
  for (i in 1:num_iter) {
    # Save the current state and its likelihood.
    path[i, ] <- curr
    liks[i] <- curr_lik
      
    cat('[', 
        format(round(curr[1], 3), nsmall = 3), '\t', 
        format(round(curr[2], 5), nsmall = 5), ']\t',
        format(round(curr_lik, 3), nsmall = 3), '\t', sep = '')
    
    # Get a proposed state and calculate its likelihood.
    prop <- q(curr)
    prop_lik <- likelihood(prop, tgt, num_sim)
    
    cat('[', 
        format(round(prop[1], 3), nsmall = 3), '\t', 
        format(round(prop[2], 5), nsmall = 5), ']\t',
        format(round(prop_lik, 3), nsmall = 3), '\t', sep = '')
      
    # Compute the ratio of the scores of the two states and generate a uniform 
    # bit.
    r <- exp(prop_lik - curr_lik)
    p <- runif(1)
    
    # Print the current progress.
    cat(format(round(r, 3), nsmall = 3), '\t', format(round(p, 3), nsmall = 3), '\n', sep = '')
    
    # Transition if the proposed state is better or if the coin flip succeeds.
    if (p < r) { 
      curr <- prop
      curr_lik <- prop_lik
      
      # If the new likelihood is better than the best we've seen so far, replace 
      # the best.
      if (curr_lik > best_lik) {
        best <- curr
        best_lik <- curr_lik
      }
    }
    
    # Save the path, best state, and likelihoods so far.
    write.table(path, file = 'path.txt', row.names = F, col.names = F)
    write.table(liks, file = 'liks.txt', row.names = F, col.names = F)
    write.table(best, file = 'best.txt', row.names = F, col.names = F)
  }
  path[num_iter + 1, ] <- curr
  liks[num_iter + 1] <- curr_lik
  return(list(path, liks, best))
}

# Solve for optimal values via MCMC.
tgt <- c(0.10, 0.25)
start <- c(51, 0.10)
results <- metropolis(start, tgt, num_sim = 300, num_iter = 1000)
path <- results[[1]]
liks <- results[[2]]
best <- results[[3]]
write.table(path, file = 'path.txt', row.names = F, col.names = F)
write.table(liks, file = 'liks.txt', row.names = F, col.names = F)
write.table(best, file = 'best.txt', row.names = F, col.names = F)

