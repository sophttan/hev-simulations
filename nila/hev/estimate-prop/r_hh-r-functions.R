# estimate within hh r estimation
gc()
library(tidyverse)
library(purrr)
library(here)

# updated Wallinga-Teunis method includes prior for relative probability of 
# infection within a household to non-household infections
probability <- function(cases, index, rel_p_hh=1) {
  inc_prim <- cases$TIME
  inc_sec <- cases$TIME[index]
  hh_prim <- cases$HH
  hh_sec <- cases$HH[index]
  # calculate relative likelihood that primary case caused secondary case based on their incidences
  # serial interval is approximate
  results <- (rel_p_hh*(hh_prim==hh_sec)+(hh_prim!=hh_sec))*dnorm(inc_sec-inc_prim, mean=31.5, sd=4)
  # if primary case happens after secondary case, set probability to 0
  results[inc_prim>=inc_sec] <- 0
  if(sum(results)==0){return(results)}
  
  return(results/sum(results))
}

res_ben <- function(data){
  t <- data %>% group_by(i) %>% 
    mutate(id=1:n()) 
  
  # res stores results
  res<-NULL
  for (j in unique(t$i)) {
    # for each simulation in 1000 simulations, estimate total R and total R from households
    data <- t %>% filter(i==j) %>% mutate(id=1:n())
    R<-rep(0, nrow(data))
    R_hh<-rep(0, nrow(data))
    for (k in 1:nrow(data)) {
      # for each case in the simulated outbreak
      # find relative likelihoods of infection from all other cases in the population
      probs <- probability(data, k, 1)
      hh_probs <- probs
      # set likelihoods to 0 if not within same household
      hh_probs[data$HH!=data$HH[k]] <- 0
      # R is the sum of the likelihoods for each case
      R <- R + probs
      R_hh <- R_hh + hh_probs
    }
    res<-res %>% rbind(data%>%mutate(R=R, R_hh=R_hh))
  }
  
  res_summary = (res %>% group_by(i) %>% summarise(R=mean(R), R_hh=mean(R_hh), mean(R_hh/R))) %>% summarise_all(mean)
  return(res_summary)
}

# data is the big dataset containing all parameter set simulations (in this case, blended)
estimate_prop = function(data) {
  prop_est = foreach(j=1:9, .combine = rbind, .export = ls(envir = globalenv()), .packages = c("here", "tidyverse")) %dopar% {
    data %>% filter(param == j) %>% res_ben() %>% mutate(param = j)
  }
  return(prop_est)
}

