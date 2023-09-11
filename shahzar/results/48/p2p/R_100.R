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
    probs <- probability(cases, k, rel_p_hh)
    HH_probs <- probs
    HH_probs[cases$HH != cases$HH[k]] <- 0
    R <- R + probs
    R_HH <- R_HH + HH_probs
  }
  return(c(mean(R_HH), mean(R), mean(R_HH / R, na.rm = T)))
    # originally return mean(R_HH) / mean(R)
}

# 5% Cumulative Incidence
n_sims <- 1000
Rs <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- read.csv(paste0('5/', i, '.csv'))
  method_R(results)
}
Rs <- matrix(Rs, n_sims, byrow = T)
saveRDS(Rs, file = '5/Rs_all_2.rds')
write.table(Rs, file = '5/Rs_all_2.txt', row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 1000
Rs <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- read.csv(paste0('10/', i, '.csv'))
  method_R(results)
}
Rs <- matrix(Rs, n_sims, byrow = T)
saveRDS(Rs, file = '10/Rs_all_2.rds')
write.table(Rs, file = '10/Rs_all_2.txt', row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 1000
Rs <- foreach (i = 1:n_sims, .combine = 'c') %dopar% {
  results <- read.csv(paste0('30/', i, '.csv'))
  method_R(results)
}
Rs <- matrix(Rs, n_sims, byrow = T)
saveRDS(Rs, file = '30/Rs_all_2.rds')
write.table(Rs, file = '30/Rs_all_2.txt', row.names = F, col.names = F)
