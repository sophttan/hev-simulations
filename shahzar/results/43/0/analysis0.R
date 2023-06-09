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
    
    beta_E <- params[1]
    risk_E <- pmin(beta_E * data$S, 1)
    
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

proportion_HH <- function(sim_data_HH,sim_data,relative=TRUE){
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
beta_E <- 0.00015
params <- c(beta_E)

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
beta_E <- 0.00031
params <- c(beta_E)

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
beta_E <- 0.00107
params <- c(beta_E)

n_sims <- 1000
vals <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- SEIR(params, inf)
  write.csv(results, file = paste0('30/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims, byrow = T)
saveRDS(vals, file = '30/vals.rds')
write.table(vals, file = '30/vals.txt', row.names = F, col.names = F)
