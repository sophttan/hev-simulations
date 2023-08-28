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
  return(method_R(results))
}

# 5% Cumulative Incidence
n_sims <- 999
vals <- foreach (i = 0:n_sims, .combine = 'c') %dopar% {
  results <- read.csv(paste0('5/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims + 1, byrow = T)
saveRDS(vals, file = '5/no_longini.rds')
write.table(vals, file = '5/no_longini.txt', row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 999
vals <- foreach (i = 0:n_sims, .combine = 'c') %dopar% {
  results <- read.csv(paste0('10/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims + 1, byrow = T)
saveRDS(vals, file = '10/no_longini.rds')
write.table(vals, file = '10/no_longini.txt', row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 999
vals <- foreach (i = 0:n_sims, .combine = 'c') %dopar% {
  results <- read.csv(paste0('30/', i, '.csv'))
  values(results)
}
vals <- matrix(vals, n_sims + 1, byrow = T)
saveRDS(vals, file = '30/no_longini.rds')
write.table(vals, file = '30/no_longini.txt', row.names = F, col.names = F)


