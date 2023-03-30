# estimate within hh r estimation

rm(list=ls())
gc()
library(tidyverse)
library(purrr)
library(EpiEstim)

setwd(here::here("results/separate_models/cumulative_inc_10/"))

inf1 <- read_csv("simulated_data/environmental.csv")
inf2 <- read_csv("simulated_data/hh_75h25c.csv")
# inf3 <- read_csv("simulated_data/hh_50h50c.csv")
inf4 <- read_csv("simulated_data/hh_25h75c.csv")

inf1
?wallinga_teunis


t <- inf2 %>% group_by(i, HH, time) %>% 
  arrange(HH, time) %>% summarise(incid=n()) %>% ungroup()

R0 <- NULL
for (j in unique(t$i)) {
  data <- t %>% filter(i==j) %>% select(!HH) %>% 
    full_join(expand.grid(i=j, time=1:364)) %>% replace_na(list(incid=0)) %>%
    arrange(time)
  res <- estimate_R(data$incid,
                    method="parametric_si", 
                    config = list(mean_si=28+4, std_si=4,
                                  t_start=2:(364-45+1),
                                  t_end=46:364))
  R0 <- rbind(R0, res$R %>% mutate(i=j))
}

(R0 %>% group_by(i) %>% summarise(meanR0=mean(`Mean(R)`,na.rm=T)))$meanR0 %>% mean()


probability <- function(inc_prim, inc_sec) {
  # calculate probability that primary case caused secondary case based on their incidences
  results <- dnorm(inc_sec-inc_prim, mean=31.5, sd=4)
  results[inc_prim>=inc_sec] <- 0
  return(results/sum(results))
}

t <- inf4 %>% group_by(i) %>% 
  mutate(id=1:n()) 

res<-NULL
for (j in unique(t$i)) {
  data <- t %>% filter(i==j) %>% mutate(id=1:n())
  primary<-NULL
  hh_primary<-NULL
  for (inc in data$time) {
    probs <- probability(data$time, inc)
    most_likely <- data$id[probs==max(probs)]
    most_likely <- ifelse(all(most_likely%>%is.na()), NA, most_likely)
    hh_most_likely <- data$HH[probs==max(probs)]
    hh_most_likely <- ifelse(all(hh_most_likely%>%is.na()), NA, hh_most_likely)
    primary <- c(primary, most_likely)
    hh_primary <- c(hh_primary, hh_most_likely)
  }
  res<-res %>% rbind(data%>%mutate(primary=primary, hh_primary=hh_primary))
}

(res %>% group_by(i, primary) %>% summarise(R0_HH=n()) %>% 
    group_by(i) %>% summarise(R0_HH=mean(R0_HH)))$R0_HH %>% mean()
(res %>% group_by(i, primary) %>% summarise(R0_HH=sum(hh_primary==HH,na.rm=T)) %>% 
    group_by(i) %>% summarise(R0_HH=mean(R0_HH)))$R0_HH %>% mean()

(res %>% mutate(same_hh=HH==hh_primary) %>% group_by(i) %>% summarise(same_hh=mean(same_hh,na.rm=T)))$same_hh%>%mean(na.rm=T)
res
