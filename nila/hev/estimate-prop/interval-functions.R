gc()
library(foreach)
library(doParallel)
library(tidyverse)
library(here)
library(latex2exp)

# choosing the number of cores to do parallelisation on
numCores = detectCores()-2
registerDoParallel(numCores)

# data should have a column labelled inc indicating the cumulative incidences: 5, 10, 30%
# parameter is the percentage of person-to-person ie: 0% for env, 100% for ptp
estimate_int_rhh <- function(data, sims, parameters, prop_true){
  est = foreach(j=c(5,10,30), .combine = rbind, .export = ls(envir = globalenv()), .packages = c("here", "tidyverse")) %dopar% {
    incidence = j
    
    estimate = data %>% filter(inc == j) %>% res_ben()
    r_hh = c(estimate$R_hh / estimate$R)[1]
    
    estimate = data %>%
      filter(inc == j) %>%
      interval_est(sims = sims) %>%
      group_by(i) %>%
      summarise(prop = mean(has_hh)) %>%
      summarise(avg_prop = mean(prop))
    interval = estimate$avg_prop
    
    est = data.frame(incidence, parameters, prop_true, r_hh, interval)
    est
  }
  return(est)
}

interval_est <- function(data, sims){
  prop_hh <- rep(0, sims)
  inc <- rep(0, sims)
  inf_type <- NULL
  for (j in 1:sims) {
    # restructure table
    f <- data %>% filter(i==j) %>% group_by(HH) %>% mutate(day_limits = list(TIME)) %>%
      ungroup() %>% rowwise() %>%
      mutate(has_hh = any((TIME - unlist(day_limits)) < 45 & (TIME - unlist(day_limits)) > 7))
    
    # checking if it's the big dataset or not
    if("param" %in% colnames(f))
    {
      f = f %>% select(c(TIME, SIZE, HH, TYPE, has_hh, param))
    } else {
      f = f %>% select(c(TIME, SIZE, HH, TYPE, has_hh))
    }
    
    if(nrow(f)==1) {
      prop_hh[j] <- NA
    }else{
      prop_hh[j] <- sum(f$TYPE=="H"|f$TYPE=="B",na.rm=T)/nrow(f)
    }
    
    inf_type <- rbind(inf_type, cbind(i=j, f))
  }
  return(inf_type)
}

# for data with multiple parameter sets, organised with column 'param'
interval_blend <- function(data, sims) {
  inf_data = NULL
  for(j in 1:9){
    # setting variable values
    blend = data %>% filter(param == j)
    inf_type <- NULL
    inf_type <- interval_est(blend, sims)
    inf_data <- rbind(inf_data, inf_type)
  }
  
  interval_prop = inf_data %>%
    group_by(param, i) %>%
    summarise(prop = mean(has_hh)) %>%
    summarise(avg_prop = mean(prop))
  return(interval_prop)
}