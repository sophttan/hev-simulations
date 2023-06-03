set.seed(1234)

library(tidyverse)
library(here)

# calculating the true prop value
# inputs:
# data is the big dataset containing all parameter set simulations
true_prop = function(data) {
  # prop values for each simulation of each parameter set
  prop_i = data %>% group_by(param, i) %>% summarise(prop = mean(str_detect(TYPE, "H")))
  # averaging the prop values for each parameter set
  prop_avg = prop_i %>% group_by(param) %>% summarise(avg_prop = mean(prop))
  return(prop_avg)
}

true_hh = function(data){
  prop = data %>%
    group_by(i) %>%
    summarise(prop = mean(str_detect(TYPE, "H"))) %>%
    summarise(avg_prop = mean(prop))
  return(prop$avg_prop)
}
