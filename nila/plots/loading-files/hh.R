library(tidyverse)
library(here)

# load_env_inc - loading person-to-person & environmental data sets
# ~inputs~
# hh_relpath: vector of 3 different strings in the incidence order 5, 10, 30
#              for person-to-person transmission only
#~ouput~
# combined dataset with column i_percent for final cumulative incidence

load_hh <- function(hh_relpath){
  inc_5 = read.csv(hh_relpath[1]) %>% mutate(i_percent = 5)
  inc_10 = read.csv(hh_relpath[2]) %>% mutate(i_percent = 10)
  inc_30 = read.csv(hh_relpath[3]) %>% mutate(i_percent = 30)
  inc = rbind(inc_5, inc_10, inc_30)
  inc$i_percent <- as.factor(inc$i_percent)
  inc
}

# loading the actual person-to-person example data
hh_relpath = here("nila", "plots", "data", c("inc_5.csv", "inc_10.csv", "inc_30.csv"))
inc = load_hh(hh_relpath)