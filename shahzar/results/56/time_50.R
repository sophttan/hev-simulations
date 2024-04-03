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

probability <- function(cases, index, rel_p_hh = 1) {
  inc_prim <- cases$TIME
  inc_sec <- cases$TIME[index]
  hh_prim <- cases$HH
  hh_sec <- cases$HH[index]
  # calculate relative likelihood that primary case caused secondary case based on their incidences
  # serial interval is approximate
  # try scaling rel_p_hh - 1/incidence or 2, 5, 10
  #results <- (rel_p_hh * (hh_prim == hh_sec) + (hh_prim != hh_sec)) * dnorm(inc_sec - inc_prim, mean = 31.5, sd = 4)
  if (any(hh_prim == hh_sec)) {
    # If there are household cases, use only those
    results <- (hh_prim == hh_sec)
  } else {
    # Otherwise, assign weights normally.
    results <- dnorm(inc_sec - inc_prim, mean = 31.5, sd = 4)
  }
  # if primary case happens after secondary case, set probability to 0
  results[inc_prim >= inc_sec] <- 0
  if(sum(results) == 0){
    return(results)
  }
  
  return(results/sum(results))
}

method_R <- function(results, rel_p_hh = 1) {
  cases <- results %>% filter(!is.na(TIME))
  if (nrow(cases) < 1) { 
    return(NaN)  
  }
  
  cases <- cases %>%
    mutate(ID = 1:n())
  
  R <- rep(0, nrow(cases))
  R_HH <- rep(0, nrow(cases))
  
  for (k in 1:nrow(cases)) {
    probs <- probability(cases, k, rel_p_hh)
    HH_probs <- probs
    HH_probs[cases$HH != cases$HH[k]] <- 0
    R <- R + probs
    R_HH <- R_HH + HH_probs
  }
  # R_HH = c * R_HH where c depends on the incidence  
  #return(c(mean(R_HH), mean(R), mean(R_HH / R, na.rm = T)))
  #return(mean(R_HH / R, na.rm = T))
  return(c(mean(R_HH), mean(R)))
}

#values <- function(results) {
#  return(c(method_R(results, rel_p_hh = ), method_R(results, rel_p_hh = 2), method_R(results, rel_p_hh = 5)))
#}

values <- function(results) {
  return(method_R(results))
}

#############
## 50% P2P ##
##############

save_dir <- '50/'
file_dir <- '../50/50/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  time_vals <- array(NA, dim = c(365, 2))
  for (t in 1:365) {
    time_vals[t, ] <- values(results %>% filter(TIME <= t))
  }
  time_vals
}
vals <- aperm(array(vals, dim = c(365, 2500, 2)), c(2, 1, 3))
saveRDS(vals, file = paste0(save_dir, 'R_HH_t_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_t_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  time_vals <- array(NA, dim = c(365, 2))
  for (t in 1:365) {
    time_vals[t, ] <- values(results %>% filter(TIME <= t))
  }
  time_vals
}
vals <- aperm(array(vals, dim = c(365, 2500, 2)), c(2, 1, 3))
saveRDS(vals, file = paste0(save_dir, 'R_HH_t_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_t_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  time_vals <- array(NA, dim = c(365, 2))
  for (t in 1:365) {
    time_vals[t, ] <- values(results %>% filter(TIME <= t))
  }
  time_vals
}
vals <- aperm(array(vals, dim = c(365, 2500, 2)), c(2, 1, 3))
saveRDS(vals, file = paste0(save_dir, 'R_HH_t_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_t_30.txt'), row.names = F, col.names = F)