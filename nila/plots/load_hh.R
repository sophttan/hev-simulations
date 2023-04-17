setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

load_hh <- function(hh_relpath){
  inc_5 = read.csv(hh_relpath[1]) %>% mutate(i_percent = 5)
  inc_10 = read.csv(hh_relpath[2]) %>% mutate(i_percent = 10)
  inc_30 = read.csv(hh_relpath[3]) %>% mutate(i_percent = 30)
  inc = rbind(inc_5, inc_10, inc_30)
  inc$i_percent <- as.factor(inc$i_percent)
  inc
}