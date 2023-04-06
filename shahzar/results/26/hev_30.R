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
likelihood <- function(state, tgt, n = 300, verbose = T) {
  # If either parameter is nonpositive, do not transition to that state.
  if (any(state <= 0)) {
    return(-Inf)
  }
  # Otherwise, find the average incidence and SAR and compute likelihood.
  vals <- foreach (i = 1:n, .combine = c) %dopar% {
    results <- SEIR(state, inc, inf)
    metrics(results)
  }
  vals <- matrix(vals, n, byrow = T)
  fit <- colMeans(vals)
  
  if (verbose) {
    cat('(', 
        format(fit[1], digits = 4, nsmall = 4), ' ',
        format(fit[2], digits = 4, nsmall = 4), ')\t', sep = '')
  }
      
  return(-log(score(fit, tgt)))
}

# Proposal function
q <- function(state, sds = c(0.05, 0.0005)) {
  # Sample from a multivariate normal distributions centered at the current 
  # state. The SDs roughly correspond to the step-size of the chain for each 
  # parameter.
  return(rnorm(n = 2, mean = state, sd = sds))
}

# MCMC
metropolis <- function(start, tgt, num_sim, num_iter) {
  path <- matrix(NA, num_iter + 1, 2)
  liks <- rep(NA, num_iter + 1)
  
  cat('START:', '\t', sep = '')
  
  # Initialize current state.
  curr <- start
  curr_lik <- likelihood(curr, tgt, num_sim)
  
  cat('\n')
  
  # Initialize best state.
  best <- curr
  best_lik <- curr_lik
  for (i in 1:num_iter) {
    # Save the current state and its likelihood.
    path[i, ] <- curr
    liks[i] <- curr_lik
      
    cat(i, '\t[', 
        format(curr[1], digits = 5, nsmall = 3), ' ', 
        format(curr[2], digits = 3, nsmall = 4), ']\t', 
        format(curr_lik, digits = 3, nsmall = 3), '\t', sep = '')
    
    # Get a proposed state and calculate its likelihood.
    prop <- q(curr)
    
    cat('[', 
        format(prop[1], digits = 5, nsmall = 3), ' ', 
        format(prop[2], digits = 3, nsmall = 4), '] ', sep = '')
    
    prop_lik <- likelihood(prop, tgt, num_sim)
    
    cat(format(prop_lik, digits = 3, nsmall = 3), '\t', sep = '')
      
    # Compute the ratio of the scores of the two states and generate a uniform 
    # bit.
    r <- exp(prop_lik - curr_lik)
    p <- runif(1)
    
    # Print the current progress.
    cat(format(r, digits = 3, nsmall = 3), '\t', 
        format(p, digits = 3, nsmall = 3), '\n', sep = '')
    
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
    write.table(path, file = '30/path.txt', row.names = F, col.names = F)
    write.table(liks, file = '30/liks.txt', row.names = F, col.names = F)
    write.table(best, file = '30/best.txt', row.names = F, col.names = F)
  }
  path[num_iter + 1, ] <- curr
  liks[num_iter + 1] <- curr_lik
  return(list(path, liks, best))
}

# Solve for optimal values via MCMC.
target <- c(0.30, 0.25)
start <- c(52.6740892462374, 0.121665484115338)
results <- metropolis(start, target, num_sim = 2000, num_iter = 200)
path <- results[[1]]
liks <- results[[2]]
best <- results[[3]]
write.table(path, file = '30/path.txt', row.names = F, col.names = F)
write.table(liks, file = '30/liks.txt', row.names = F, col.names = F)
write.table(best, file = '30/best.txt', row.names = F, col.names = F)
