# estimate within hh r estimation

rm(list=ls())
gc()
library(tidyverse)
library(purrr)
library(EpiEstim)

setwd(here::here("results"))

#data <- read_csv("hev-comparison/data_inc_30.csv") %>% rename("time"="TIME")
data <- read_csv("separate_models/cumulative_inc_5/simulated_data/environmental.csv")

data
?wallinga_teunis

probability <- function(inc_prim, inc_sec) {
  # calculate relative likelihood that primary case caused secondary case based on their incidences
  # serial interval is approximate
  results <- dnorm(inc_sec-inc_prim, mean=31.5, sd=4)
  # if primary case happens after secondary case, set probability to 0
  results[inc_prim>=inc_sec] <- 0
  if(sum(results)==0){return(results)}

  return(results/sum(results))
}

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
    probs <- probability(data$time, data$time[k])
    hh_probs <- probs
    # set likelihoods to 0 if not within same household - but should this be reweighted instead?
    hh_probs[data$HH!=data$HH[k]] <- 0
    # R is the sum of the likelihoods for each case
    R <- R + probs
    R_hh <- R_hh + hh_probs
  }
  res<-res %>% rbind(data%>%mutate(R=R, R_hh=R_hh))
}

res

(res %>% group_by(i) %>% summarise(R=mean(R), R_hh=mean(R_hh), R_hh/R)) %>% summarise_all(mean)


# t <- inf2 %>% group_by(i, HH, time) %>% 
#   arrange(HH, time) %>% summarise(incid=n()) %>% ungroup()
# 
# R0 <- NULL
# for (j in unique(t$i)) {
#   data <- t %>% filter(i==j) %>% select(!HH) %>% 
#     full_join(expand.grid(i=j, time=1:364)) %>% replace_na(list(incid=0)) %>%
#     arrange(time)
#   res <- estimate_R(data$incid,
#                     method="parametric_si", 
#                     config = list(mean_si=28+4, std_si=4,
#                                   t_start=2:(364-45+1),
#                                   t_end=46:364))
#   R0 <- rbind(R0, res$R %>% mutate(i=j))
# }
# 
# (R0 %>% group_by(i) %>% summarise(meanR0=mean(`Mean(R)`,na.rm=T)))$meanR0 %>% mean()
