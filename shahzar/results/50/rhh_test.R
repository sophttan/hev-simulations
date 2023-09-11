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

probability <- function(cases, index, rel_p_hh = 1) {
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

proportion_HH <- function(results, relative = T){
  sim_data <- results %>% filter(!is.na(TIME))
  sim_data_HH <- results %>% select(HH, SIZE) %>% distinct(HH, .keep_all = T)
  
  # household sizes in the data
  k <- sort(unique(sim_data_HH$SIZE))
  # how many households of each size
  n <- table(sim_data_HH$SIZE)
  # matrix with rows=number of infections, cols=household size
  a <- matrix(0,nrow=max(k)+1,ncol=length(k))
  if (length(unique(sim_data$HH))>1){
    cases_per_HH <- table(sim_data$HH,sim_data$SIZE)
  } else {
    cases_per_HH <- table(sim_data$HH)
    cases_per_HH <- t(cases_per_HH)
    colnames(cases_per_HH) <- unique(sim_data$SIZE)[[1]]
  }
  for (idx in 1:length(k)){
    size <- k[idx]
    if (as.character(size) %in% colnames(cases_per_HH)){
      case_counts <- table(cases_per_HH[,as.character(size)])
      for (cases in labels(case_counts)[[1]]){
        a[as.integer(cases)+1,idx] <- case_counts[cases]
      }
    }
  }
  # above loop gives incorrect zero-cases, instead find zero-case
  # households by finding households not included in case data
  a[1,] <- 0
  for (idx in 1:length(k)){
    size <- k[idx]
    for (HH in sim_data_HH[sim_data_HH$SIZE==size,]$HH){
      if (!(HH %in% sim_data$HH)){
        a[1,idx] <- a[1,idx] + 1
      }
    }
  }
  # B_est is probability of escaping infection from community transmission
  B_est <- sum(n*(a[1,]/n)**(1/k))/sum(n)
  
  # phi is avg number infected per household
  phi <- sum((0:max(k))*rowSums(a))/sum(n)
  # theta is the household attack rate
  theta <- sum((0:max(k))*rowSums(t(t(a)/k)))/sum(n)
  # Q_est is probability of escaping infection from household transmission
  # estimator can give probabilities greater than 1 so set max of 1
  Q_est <- min(1,((1-theta)/B_est)**(1/phi))
  
  # give relative probability of household infectiohn
  if (relative){
    return((1-Q_est)/(1-B_est))
  }
  
  # proportion of household infection expressed as conditional probability
  # of household infection given infection at all
  p_HH <- (1-Q_est)/(1-B_est*Q_est)
  return(p_HH)
}

method_R <- function(results) {
  cases <- results %>% 
    filter(!is.na(TIME)) %>%
    mutate(ID = 1:n())
  
  R <- rep(0, nrow(cases))
  R_HH <- rep(0, nrow(cases))
  
  for (k in 1:nrow(cases)) {
    rel_p_hh <- proportion_HH(results)
    #probs <- probability(cases, k, rel_p_hh)
    probs <- probability(cases, k)
    HH_probs <- probs
    HH_probs[cases$HH != cases$HH[k]] <- 0
    R <- R + probs
    R_HH <- R_HH + HH_probs
  }
  ratio <- mean(R_HH / R, na.rm = T)
  cat(paste0(format(round(mean(ratio), 2), nsmall = 2), '\t'))
  return(ratio)
}

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

SEIR_p2p <- function(params, inf = 7, verbose = F) {
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
      
      results$I_num <- results$I_num + I_data$new_I_H
      
      # Label infection types.
      results$TYPE[new_inf_C == 1] <- 'C'
      results$TYPE[new_inf_H == 1] <- 'H'
      results$TYPE[(new_inf_H == 1) & (new_inf_C == 1)] <- 'HC'
    }
    
    # Increment exposure and infectious counters.
    data$E_count[data$E == 1] <- data$E_count[data$E == 1] + 1
    data$I_count[data$I == 1] <- data$I_count[data$I == 1] + 1
  }
  return(results)
}

R_seq <- function(results) {
  Rs <- rep(NA, 365)
  min_t <- min(results$TIME, na.rm = T) + 1
  for (t in min_t:365) {
    Rs[t] <- method_R(results %>% filter(TIME < t))
  }
  cat('\n')
  return(Rs)
}

# 30% Cumulative Incidence
n_sims <- 10
params <- c(51, 0.140)
vals <- array(rep(NA, n_sims*365), dim = c(n_sims, 365))
for (i in 1:n_sims) {
  results <- SEIR_p2p(params)
  vals[i, ] <- R_seq(results)
}
saveRDS(vals, file = 'vals_30.rds')
write.table(vals, file = 'vals_30.txt', row.names = F, col.names = F)

# 10% Cumulative Incidence
n_sims <- 10
params <- c(51, 0.100)
vals <- array(rep(NA, n_sims*365), dim = c(n_sims, 365))
for (i in 1:n_sims) {
  results <- SEIR_p2p(params)
  vals[i, ] <- R_seq(results)
}
saveRDS(vals, file = 'vals_10.rds')
write.table(vals, file = 'vals_10.txt', row.names = F, col.names = F)
