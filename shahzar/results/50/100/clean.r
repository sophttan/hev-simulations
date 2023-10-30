rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

message(detectCores())
num_cores <- 32
registerDoParallel(num_cores)

j <- 1
for (i in 1:100000) {
    ld <- paste0('100/5/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('100/5_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('100/10/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('100/10_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('100/30/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('100/30_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}