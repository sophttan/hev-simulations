rm(list = ls())
gc()
library(dplyr)
library(foreach)
library(doParallel)

message(detectCores())