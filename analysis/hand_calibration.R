rm(list = ls())
gc()
library(dplyr)

source(here::here("analysis/functions/seir_functions.R"))

#########################
#### SEIR Simulation ####
#########################
time <- 365 # Number of days.
inc <- 28 # Average incubation period length.
inf <- 7 # Average infectious period length.
N <- 1000 # Population size.

start <- c(53.6136261376403, 0.086607502588279)
incidence<-NULL
sar<-NULL
results_data<-NULL

sims <- 100
for (i in 1:sims) {
  results <- SEIR(start, inc, inf)
  res <- metrics(results)
  incidence<-c(incidence,res[1])
  sar<-c(sar,res[2])
  
  results <- results %>% filter(!is.na(TIME)) %>% mutate(i=i)
  results_data <- rbind(results_data, results)
}

results_data
mean(incidence)
mean(sar)
