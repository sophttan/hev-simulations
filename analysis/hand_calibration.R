rm(list = ls())
gc()
library(dplyr)
library(here)
set.seed(1234)
source(here::here("analysis/functions/seir_functions.R"))

#########################
#### SEIR Simulation ####
#########################
time <- 365 # Number of days.
inf <- 7 # Average infectious period length.
N <- 1000 # Population size.

start <- c(32.5,	0.075,	0.0005) # for 10% idc and 67.2% prp
incidence<-NULL
sar<-NULL
results_data<-NULL

sims <- 500
for (i in 1:sims) {
  results <- SEIR_blend(start, inf)
  res <- metrics_blend(results)
  incidence<-c(incidence,res[1])
  sar<-c(sar,res[2])
  
  results <- results %>% filter(!is.na(TIME)) %>% mutate(i=i)
  results_data <- rbind(results_data, results)
}

write.csv(results_data, file = here("nila", "hev", "blend-sims", "idc-10-prp-67.csv"))
mean(incidence)
mean(sar)
