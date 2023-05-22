rm(list = ls())
gc()
library(dplyr)
library(foreach)
library(doParallel)

inf <- 7 # Average infectious period length.

source(here::here("analysis/functions/seir_functions.R"))
source(here::here("analysis/functions/mcmc_functions.R"))

# Set up the number of cores used for parallelization.
# Use detectCores() to find out how many cores are available.
num_cores <- 1
registerDoParallel(num_cores)

args <- commandArgs(trailingOnly = T) %>% as.numeric()
incidence <- args[1]/100
sar <- args[2]/100
target <- c(incidence, sar)

if(incidence==0.05) {
  start<-c(56, 0.07)
}else if(incidence==0.1){
  start<-c(60.3, 0.08)
}else if(incidence==0.3){
  start<-c(53.7, 0.12)
}

results <- metropolis(start, target, num_sim = 300, num_iter = 1000)
path <- results[[1]]
liks <- results[[2]]
best <- results[[3]]
write.table(path, file = paste0('path_inc_', incidence, '.txt'), row.names = F, col.names = F)
write.table(liks, file = paste0('liks_inc_', incidence, '.txt'), row.names = F, col.names = F)
write.table(best, file = paste0('best_inc_', incidence, '.txt'), row.names = F, col.names = F)