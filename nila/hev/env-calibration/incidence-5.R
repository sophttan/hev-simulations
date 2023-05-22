rm(list = ls())
gc()
library(dplyr)
library(here)

source(here::here("nila/hev/env-calibration/seir_functions.R"))

#########################
#### SEIR Simulation ####
#########################
time <- 365 # Number of days.
inf <- 7 # Average infectious period length.
N <- 1000 # Population size.

start <- 0.1399777/N
incidence<-NULL
results_data<-NULL
demg_data<-NULL

sims <- 1000
for (i in 1:sims) {
  seir <- SEIR_environment(start, inf)
  results <- seir[[1]]
  colnames(results)[colnames(results) == 'time'] <- 'TIME'
  demg <- seir[[2]]
  res <- metrics(results)[1]
  incidence<-c(incidence,res)
  
  results <- results %>% filter(!is.na(TIME)) %>% mutate(i=i)
  results_data <- rbind(results_data, results)
  
  demg <- demg %>% mutate(i=i)
  demg_data <- rbind(demg_data, demg)
}

# to store only iteration, HH_size, and HH ID
demg_data = demg_data[!duplicated(demg[3:5]),] %>% select(3:5)

write.csv(results_data, here::here("nila/hev/env-calibration/data/env-results-5.csv"))
write.csv(demg_data, here::here("nila/hev/env-calibration/data/env-demg-5.csv"))
mean(incidence)
