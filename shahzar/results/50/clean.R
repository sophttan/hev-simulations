rm(list = ls())
gc()
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

###############
## Panmictic ##
###############

# 5% Incidence
# Already cleaned

# 10% Incidence
# Already cleaned

# 30% Incidence
# Already cleaned


########
## 0% ##
########

# 5% Incidence
# Already cleaned

# 10% Incidence
# Already cleaned

# 30% Incidence
# Already cleaned

#########
## 25% ##
#########

j <- 1
for (i in 1:100000) {
    ld <- paste0('25/5/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('25/5_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('25/10/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('25/10_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('25/30/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('25/30_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

#########
## 50% ##
#########

j <- 1
for (i in 1:100000) {
    ld <- paste0('50/5/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('50/5_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('50/10/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('50/10_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('50/30/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('50/30_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

#########
## 75% ##
#########

j <- 1
for (i in 1:100000) {
    ld <- paste0('75/5/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('75/5_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('75/10/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('75/10_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

j <- 1
for (i in 1:100000) {
    ld <- paste0('75/30/', i, '.csv')
    if (file.exists(ld)) {
        results <- read.csv(ld)
        sv <- paste0('75/30_clean/', j, '.csv')
        write.csv(results, file = sv, row.names = F)
        j <- j + 1
    }
}

##########
## 100% ##
##########

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