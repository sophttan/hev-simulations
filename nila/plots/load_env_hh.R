setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

# load_env_inc - loading person-to-person & environmental data sets
# ~inputs~
# hh_relpath: vector of 3 different strings in the incidence order 5, 10, 30
#              for person-to-person transmission only
# env_relpath: vector of 3 different strings in the incidence order 5, 10, 30
#              for environmental transmission only
#~ouput~
# list containing two datasets: inc, env

load_env_hh <- function(hh_relpath, env_relpath){
  # loading person-to-person dataset
  # i_percent identifies their cumulative incidence
  inc_5 = read.csv(hh_relpath[1]) %>% mutate(i_percent = 5)
  inc_10 = read.csv(hh_relpath[2]) %>% mutate(i_percent = 10)
  inc_30 = read.csv(hh_relpath[3]) %>% mutate(i_percent = 30)
  
  # binding together all three person-to-person files
  inc = rbind(inc_5, inc_10, inc_30)
  inc$i_percent <- as.factor(inc$i_percent)
  # to distinguish environment & person-to-person data
  inc$enviro <- rep("Person-to-person", nrow(inc))
  
  # loading environmental datasets
  env_5 = read.csv(env_relpath[1]) %>% mutate(i_percent = 5)
  env_10 = read.csv(env_relpath[2]) %>% mutate(i_percent = 10)
  env_30 = read.csv(env_relpath[3]) %>% mutate(i_percent = 30)
  
  # binding together all three environment files
  env = rbind(env_5, env_10, env_30)
  env$i_percent <- as.factor(env$i_percent)
  env$enviro <- rep("Environmental", nrow(env))
  list(inc, env)
}

# # uncomment to run example code
# hh_relpath = c("data_inc_5.csv", "data_inc_10.csv", "data_inc_30.csv")
# env_relpath =
#   c("../results/separate_models/cumulative_inc_5/simulated_data/environmental.csv",
#     "../results/separate_models/cumulative_inc_10/simulated_data/environmental.csv",
#     "../results/separate_models/cumulative_inc_30/simulated_data/environmental.csv")
# 
# hh_env = load_env_hh(hh_relpath, env_relpath)
# hh = hh_env[[1]]
# env = hh_env[[2]]
