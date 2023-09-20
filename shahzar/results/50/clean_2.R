rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)



##########
## 100% ##
##########

dir <- '100/'

n_sims <- 100000
obs <- foreach (i = 1:n_sims) %dopar% {
  f <- paste0(dir, '5/', i, '.csv')
  if (file.exists(f)) {
    results <- read.csv(f)
    colnames(results) <- c('ID', 'SIZE', 'HH', 'TYPE', 'TIME', 'S_NUM', 'I_NUM')
    write.csv(results, file = f, row.names = F)
  }
}

n_sims <- 100000
obs <- foreach (i = 1:n_sims) %dopar% {
  f <- paste0(dir, '10/', i, '.csv')
  if (file.exists(f)) {
    results <- read.csv(f)
    colnames(results) <- c('ID', 'SIZE', 'HH', 'TYPE', 'TIME', 'S_NUM', 'I_NUM')
    write.csv(results, file = f, row.names = F)
  }
}

n_sims <- 100000
obs <- foreach (i = 1:n_sims) %dopar% {
  f <- paste0(dir, '30/', i, '.csv')
  if (file.exists(f)) {
    results <- read.csv(f)
    colnames(results) <- c('ID', 'SIZE', 'HH', 'TYPE', 'TIME', 'S_NUM', 'I_NUM')
    write.csv(results, file = f, row.names = F)
  }
}