library(tidyverse)

time = 365
inc = 28
inf = 7
pop = 1000 # Population size

create_hh <- function() {
  # Randomly sample household sizes such that total population is 1000 individuals.
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
  hh_size = create_hh()
  
  # Create frame for running the simulation
  # ID: ID of individual
  # SIZE: size of individual's household
  # HH: ID of individual's household
  # S: susceptibility status
  # E: exposed status
  # E_count: number of days since exposed
  # I: infectious status
  # I_count: number of days since infectious
  # R: recovered status
  # INC: incubation period
  # INF: infectious period
  data = data.frame(ID = 1:pop,
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
  
  # Create frame for storing results
  # ID: ID of individual
  # SIZE: size of individual's household
  # HH: ID of individual's household
  # TYPE: the kind of infection: household (H), community (C), or both (B)
  # TIME: when the individual became infectious
  # S_num: number of susceptible people in individual's household when their infectious period begins
  # I_num: number of people in household that this individual infected over their infectious period
  results <- data[, 1:3] %>% mutate(TYPE = NA, TIME = NA, S_num = NA, I_num = 0)
  results$TYPE[1] = '0'
  
  for(t in 1:time) {
    if (verbose) {
      if (t %% 10 == 0) {
        print(t)
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
    num_new_inf = sum(new_inf, na.rm = T)
    if(num_new_inf > 0) {
      # Change status to newly infectious and add infectious period.
      data$I[new_inf] <- 1
      random_inf <- rnorm(num_new_inf, mean = inf, sd = 1) %>% round()
      data$INF[new_inf] <- random_inf
      
      # Remove exposure status and exposure count.
      data$E[new_inf] <- 0
      data$E_count[new_inf] <- 0 
      
      # Record time at which infectious period starts.
      results$TIME[new_inf == 1] <- t
      
      # Save the number of susceptible people in each infectious 
      # individual's household.
      S_data = data %>% group_by(HH) %>% mutate(S_tot = sum(S)) %>% select(HH, S_tot)
      results$S_num[new_inf == 1] = S_data$S_tot[new_inf == 1]
    }
    
    # I_H is the number of infections inside each household.
    # I_C is the number of infections outside each household.
    I_data <- data %>% group_by(HH) %>% 
      mutate(I_H = sum(I)) %>% 
      ungroup() %>% 
      mutate(I_C = sum(I) - I_H)
    
    # Calculate household risk and community risk.
    risk_H <- beta_H * data$S * I_data$I_H / pop
    risk_C <- beta_C * data$S * I_data$I_C / pop
    
    # Each individual is infected from their household or 
    # community independently with probabilities risk_H
    # and risk_C.
    new_inf_H <- rbinom(nrow(data), 1, risk_H)
    new_inf_C <- rbinom(nrow(data), 1, risk_C)
    
    new_exposed = (new_inf_H == 1) | (new_inf_C == 1)
    num_new_exposed = sum(new_exposed, na.rm = T)
    if (sum(new_exposed) > 0) {
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
      rr <- results %>% filter(new_inf_H == 1) %>%
        group_by(HH) %>%
        mutate(I_tot = sum(TYPE == 'H')) %>%
        ungroup()
      rr <- unique(rr[, c('HH', 'I_tot')])
      
      # Get people with the smallest I_counts and the households with the
      # new infections.
      dd = data %>% filter(data$I == 1) %>% group_by(HH) %>% slice(which.min(I_count))
      min_IDs <- dd$ID
      new_HHs <- results[new_inf_H == 1, ]$HH
      
      # Every individual with the smallest I_count in each household with 
      # at least one new H infection gets the number of new H infections
      # added to their I_num.
      idx <- (data$ID %in% min_IDs) & (data$HH %in% new_HHs)
      
      # For debugging purposes.
      if (length(results[idx, ]$I_num) != length(rr$I_tot)) {
        return(list(data, results, new_inf_H))
      }
      
      results[idx, ]$I_num <- results[idx, ]$I_num + rr$I_tot
      
      # Label individuals with both a household and community infection with B.
      results$TYPE[(new_inf_H == 1) & (new_inf_C == 1)] <- 'B'
    }
    # Increment exposure and infectious counters.
    data$E_count[data$E == 1] <- data$E_count[data$E == 1] + 1
    data$I_count[data$I == 1] <- data$I_count[data$I == 1] + 1
  }
  return(results)
}

metrics = function(results) {
  # Incidence is the proportion of the population that became infected.
  idc = mean(!is.na(results$TYPE))
  
  # If incidence is 0, the SAR is undefined.
  sar = NA
  if (idc != 0) {
    # The SAR is the average SAR for each individual that was infectious.
    sar = mean(results$I_num / results$S_num, na.rm = T)
  }
  return(c(idc, sar))
}

score = function(obs, target) {
  # The score is the L2 distance of the observed values from the target.
  return(sum((obs - target)^2))
}

# The likelihood is calculated by first averaging the incidence and SAR over N
# simulations with the state parameters. The likelihood is the negative log
# score of the average incidence and SAR.
likelihood = function(state) {
  beta_H = state[1]
  beta_C = state[2]
  
  vals = matrix(0, N, 2)
  for (i in 1:N) {
    results = SEIR(beta_H, beta_C, inc, inf)
    vals[i, ] = metrics(results)
  }
  
  avg_vals = colMeans(vals)
  print(avg_vals)
  return(-log(score(avg_vals, target)))
}

#### Metropolis algorithm ####

# Proposal function
q = function(state) {
  beta_H = state[1]
  beta_C = state[2]
  r = c(beta_H^2, beta_C^2 / 0.0001)
  v = c(beta_H^2, beta_C^2 / 0.0001)
  
  return(pmax(rgamma(n = 2, shape = r, rate = v), 1e-3))
}

# MCMC
metropolis = function(start, num_iter) {
  chain = matrix(0, num_iter + 1, 2)
  liks = matrix(0, num_iter + 1, 2)
  
  # Initialize current state.
  curr = start
  curr_lik = likelihood(curr)
  
  # Initialize best state.
  best = curr
  best_lik = curr_lik
  for (i in 1:num_iter) {
    # Save the current state and its likelihood.
    chain[i, ] = curr
    liks[i, ] = curr_lik
    
    # Print current state and likelihood.
    paste0(i, '\t', curr, " ", signif(curr_lik, 3), '\t')
    
    # Get a proposed state and calculate its likelihood.
    prop = q(curr)
    prop_lik = likelihood(prop)
    paste0(prop, " ", signif(prop_lik, 3), '\t')
    
    # Compute the ratio of the scores of the two states
    # and flip a coin.
    r = exp(prop_lik - curr_lik)
    p = runif(1)
    paste(signif(r, 3), signif(p, 3))
    
    # Transition if the proposed state is better or
    # if the coin flip succeeds.
    if (p < r) { 
      curr = prop
      curr_lik = prop_lik
      
      # If the new likelihood is better than the
      # best we've seen so far, replace the best.
      if (curr_lik > best_lik) {
        best = curr
        best_lik = curr_lik
      }
    }
    
    # Save the chain, best state, and likelihoods
    # so far.
    save(chain, file = 'chain.Rdata')
    save(liks, file = 'liks.Rdata')
    save(best, file = 'best.Rdata')
  }
  return(list(chain, liks, best))
}

# Solve for optimal values via MCMC.
target = c(0.3, 0.25) # Target values.
N = 300 # Number of times over which to average likelihood.

#metropolis_results = metropolis(c(30, 0.12), 10)
#chain = metropolis_results[[1]]
#liks = metropolis_results[[2]]
#best = metropolis_results[[3]]
#save(chain, file = "chain.Rdata")
#save(liks, file = "liks.Rdata")
#save(best, file = "best.Rdata")

t_0 = Sys.time()
beta_H = 30
beta_C = 0.15
N = 1000
vals = matrix(0, N, 2)
for (i in 1:N) {
  results = SEIR(beta_H, beta_C, inc, inf)
  vals[i, ] = metrics(results)
  save(vals, file = 'vals.Rdata')
}
t_1 = Sys.time()
print(t_1 - t_0)
