rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

# Set up the number of cores used for parallelization.
message(detectCores())
num_cores <- 24
registerDoParallel(num_cores)

probability <- function(cases, index, rel_p_hh = 1) {
  inc_prim <- cases$TIME
  inc_sec <- cases$TIME[index]
  hh_prim <- cases$HH
  hh_sec <- cases$HH[index]
  # calculate relative likelihood that primary case caused secondary case based on their incidences
  # serial interval is approximate
  # try scaling rel_p_hh - 1/incidence or 2, 5, 10
  results <- (rel_p_hh * (hh_prim == hh_sec) + (hh_prim != hh_sec)) * dnorm(inc_sec - inc_prim, mean = 31.5, sd = 4)
  # if primary case happens after secondary case, set probability to 0
  results[inc_prim >= inc_sec] <- 0
  if(sum(results) == 0){
    return(results)
  }
  
  return(results/sum(results))
}

method_R <- function(results, rel_p_hh = 1) {
  cases <- results %>% 
    filter(!is.na(TIME)) %>%
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
  return(mean(R_HH))
}

#values <- function(results) {
#  return(c(method_R(results, rel_p_hh = ), method_R(results, rel_p_hh = 2), method_R(results, rel_p_hh = 5)))
#}

values <- function(results) {
  return(method_R(results))
}

###############
## Panmictic ##
###############

save_dir <- 'pan/'
file_dir <- '../50/pan/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_30.txt'), row.names = F, col.names = F)




##############
### Pure E ###
##############

save_dir <- '0/'
file_dir <- '../50/0/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_30.txt'), row.names = F, col.names = F)




###############
##  25% P2P  ##
###############

save_dir <- '25/'
file_dir <- '../50/25/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_30.txt'), row.names = F, col.names = F)



###############
##  50% P2P  ##
###############

save_dir <- '50/'
file_dir <- '../50/50/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_30.txt'), row.names = F, col.names = F)




###############
##  75% P2P  ##
###############

save_dir <- '75/'
file_dir <- '../50/75/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_30.txt'), row.names = F, col.names = F)




##############
## 100% P2P ##
##############

save_dir <- '100/'
file_dir <- '../50/100/'

# 5% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '5/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_05.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_05.txt'), row.names = F, col.names = F)


# 10% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '10/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_10.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_10.txt'), row.names = F, col.names = F)


# 30% Cumulative Incidence
n_sims <- 2500
vals <- foreach (i = 1:n_sims, .combine = 'rbind') %dopar% {
  f <- paste0(file_dir, '30/', i, '.csv')
  results <- read.csv(f)
  values(results)
}
vals <- matrix(vals, ncol = dim(vals)[2])
saveRDS(vals, file = paste0(save_dir, 'R_HH_30.rds'))
write.table(vals, file = paste0(save_dir, 'R_HH_30.txt'), row.names = F, col.names = F)

# Examine this and the shrinkage effect over time
# Look at earlier time periods.
#     Surmountable if it works for earlier times but not later ones
# Look at other formulations which can distinguish between E vs HH (R_HH alone?)
# Normal reproduction number method and look at cases which are "imports" (no obvious infectors)
#    Some cases won't have a plausible infector. Label those as environmental.
#    R package which computes this and has some accounting for imports
#    How would environmental look relative to person to person for pure reproduction number methods?
# With detailed spatial data, could make likelihood network based on contacts (closer vs farther)
#    Meta-population of 
# x-axis with proportion of household cases instead of person-to-person transmission
# Big issue: relation between household and person-to-person is not stable over different incidences.
#    # As long as it's household, there'll be that problem.

# Three Routes
#   1. Not a problem early in outbreaks. Later on, community cases overrun household cases. [Check R_HH/R early on]
#   2. Imported case approach works. [Quick check]
#   3. New formulation or adjustment to address this. E.g., just R_HH, or other formulation. [Check other formulations]