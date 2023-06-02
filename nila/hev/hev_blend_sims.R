set.seed(1234)

library(foreach)
library(doParallel)
library(tidyverse)
library(here)
library(latex2exp)

# choosing the number of cores to do parallelisation on
numCores = detectCores()-2
registerDoParallel(numCores)

source(here::here("analysis/functions/seir_functions.R"))
source(here::here("nila/hev/r0_estimation.R"))

#######################
# Simulating the data #
#######################

time <- 365 # Number of days.
inf <- 7 # Average infectious period length.
N <- 100 # Population size.
sims <- 5 # Number of simulations

# start = vector of values with beta_h, beta_c, beta_e
sim_blend = function(start, sims, inf) {
  incidence<-NULL
  sar<-NULL
  results_data<-NULL
  
  for (i in 1:sims) {
    results <- SEIR_blend(start, inf)
    res <- metrics_blend(results)
    incidence<-c(incidence,res[1])
    sar<-c(sar,res[2])
    
    results <- results %>% filter(!is.na(TIME)) %>% mutate(i=i)
    results_data <- rbind(results_data, results)
  }
  return(results_data)
}

# inputting the calibrated values for each incidence and then concatenating into one big matrix
# cols: beta_h, beta_c, beta_e and rows: 25%, 50%, 75%
param_5 = c(10, 0.005, 0.0001, 28, 0.005, 0.0001, 42, 0.03, 0.00005)
param_5 = matrix(param_5, nrow = 3, ncol = 3, byrow = TRUE)

param_10 = c(12, 0, 0.00025, 27, 0.005, 0.00015, 43, 0.035, 0.00005)
param_10 = matrix(param_10, nrow = 3, ncol = 3, byrow = TRUE)

param_30 = c(12, 0.005, 0.0008, 29, 0.01, 0.00055, 41, 0.04, 0.00025)
param_30 = matrix(param_30, nrow = 3, ncol = 3, byrow = TRUE)

param = rbind(param_5, param_10, param_30)

# getting the simulations
blended = foreach(j=1:9, .combine = rbind, .export = ls(envir = globalenv()), .packages = c("here", "tidyverse")) %dopar% {
  sim_blend(param[j,], sims, inf) %>% mutate(param = j)
}

######
# R0 #
######

# calculating the true r0 value
# inputs:
# data is the big dataset containing all parameter set simulations
# N is population size
true_r0 = function(data, N) {
  # r0 values for each simulation of each parameter set
  r0_i = data %>% group_by(param, i) %>% summarise(r0 = mean(TYPE != "E"))
  # averaging the r0 values for each parameter set
  r0_avg = r0_i %>% group_by(param) %>% summarise(avg_r0 = mean(r0))
  return(r0_avg)
}

r0_true = true_r0(blended, N)$avg_r0

# data is the big dataset containing all parameter set simulations (in this case, blended)
estimate_r0 = function(data) {
  r0_est = foreach(j=1:9, .combine = rbind, .export = ls(envir = globalenv()), .packages = c("here", "tidyverse")) %dopar% {
    data %>% filter(param == j) %>% res_ben() %>% mutate(param = j)
  }
  return(r0_est)
}

r = estimate_r0(blended)
r0_est = r$R_hh / r$R

parameters = c(1:3,1:3,1:3)*0.25

# 3 tables with 25%, 50%, and 70%
indices = c(5,5,5,10,10,10,30,30,30)
r0 = data.frame(indices,parameters, r0_true, r0_est)
r0_5 = r0 %>% slice(1:3)
r0_10 = r0 %>% slice(4:6)
r0_30 = r0 %>% slice(7:9)

# we need to get the environmental and only person-to-person models as well
# note that the true value for environmental and person-to-person are 0 and 1 respectively

# ~work in progress~

# # this is the part takes a really long time to load
# env_inc_5 <- read.csv(here("nila/hev/env-calibration/data/env-results-5.csv")) %>% mutate(inc = 5)
# env_inc_10 <- read.csv(here("nila/hev/env-calibration/data/env-results-10.csv")) %>% mutate(inc = 10)
# env_inc_30 <- read.csv(here("nila/hev/env-calibration/data/env-results-30.csv")) %>% mutate(inc = 30)
# env = rbind(env_inc_5, env_inc_10, env_inc_30)
# 
# env_est_r0 = c()
# for(j in c(5,10, 30)) {
#   estimate = env %>% filter(inc == j) %>% res_ben()
#   env_est_r0[j] = estimate$R
# }

############
# Plotting #
############

cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")

r0_pretty = r0 %>% pivot_longer(r0, cols = 3:4, names_to = "type")
r0_pretty$type = str_replace(r0_pretty$type, "r0_est", "estimate")
r0_pretty$type = str_replace(r0_pretty$type, "r0_true", "true")

r0_pretty %>%
  ggplot(aes(x = parameters, y = value, colour = type)) +
  geom_line() +
  facet_grid(cols = vars(indices), labeller = as_labeller(cum_labels)) +
  ylab("proportion of household cases") +
  xlab("proportion of person-to-person") +
  scale_color_discrete(name = "type") +
  ggtitle(TeX("estimating via $R_{hh}/R$ with a normally distrubted infection interval"))+
  theme_minimal()
