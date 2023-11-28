rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

message(detectCores())
num_cores <- 24
registerDoParallel(num_cores)

##########
## 100% ##
##########

j <- 1
for (i in 1:100000) {
    ld <- paste0('5/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('5_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('10/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('10_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('30/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('30_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}