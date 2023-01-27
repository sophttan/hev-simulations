# estimate within hh attack rates

rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/separate_models/cumulative_inc_10/"))

inf1 <- read_csv("simulated_data/environmental.csv")
inf2 <- read_csv("simulated_data/hh_75h25c.csv")
inf3 <- read_csv("simulated_data/hh_50h50c.csv")
inf4 <- read_csv("simulated_data/hh_25h75c.csv")

find_secondary_and_susceptible <- function(inf) {
  inf <- inf %>% arrange(i, HH, time) %>% group_by(i, HH) %>% 
    mutate(otherinfs=list(time), inf=1:n()) %>% 
    rowwise() %>% 
    mutate(secondaryinfs = list(otherinfs[unlist(otherinfs)-time < 45 & unlist(otherinfs)-time > 7])) %>%
    select(!otherinfs)
  
  inf_unnest <- inf %>% unnest(secondaryinfs, keep_empty = T) %>% ungroup()
  
  # identify secondary infections that can be attributed to many primary infections
  check <- inf_unnest %>% group_by(i, HH, secondaryinfs) %>% filter(!is.na(secondaryinfs)) %>% filter(n()>1)
  nomult <- inf_unnest %>% group_by(i, HH, secondaryinfs) %>% filter(is.na(secondaryinfs)|n()==1)
  
  # identify most likely primary case and assign - if equally likely, randomly assign
  assign_prob <- check %>% mutate(difftime=secondaryinfs-time, prob=dnorm(difftime, mean=31, sd=sqrt(5)))
  assign_single_case <- assign_prob %>% group_by(i, HH, secondaryinfs) %>% 
    arrange(i, HH, time, desc(prob)) %>% 
    filter(prob==max(prob)) %>% sample_n(size = 1, replace=F)
  
  # combine datasets
  nomult <- nomult %>% group_by(i, HH, inf) %>% 
    mutate(secondaryinfs=ifelse(n()==1 & is.na(secondaryinfs), 0, n())) %>%
    distinct(.keep_all = T)
  mult <- check %>% ungroup() %>% 
    distinct(i, HH, inf, .keep_all = T) %>% 
    full_join(assign_single_case %>% select(i, HH, inf, secondaryinfs), by=c("i", "HH", "inf"))
  mult <- mult %>% group_by(i, HH, inf) %>% 
    mutate(secondaryinfs=ifelse(n()==1 & is.na(secondaryinfs.y), 0, n())) %>% 
    select(!c(secondaryinfs.x, secondaryinfs.y)) %>% 
    distinct(.keep_all=T)
  
  final <- nomult %>% 
    full_join(mult, by=c("i", "time", "HHsize", "HH", "Type", "has_hh", "month", "inf")) %>% 
    arrange(i, HH, inf)                  
  final <- final %>% rowwise() %>% 
    mutate(secondaryinfs = sum(secondaryinfs.x,secondaryinfs.y,na.rm=T)) %>% 
    select(!c(secondaryinfs.x, secondaryinfs.y))  
  
  
  # estimate susceptible population in household
  nosecondary <- final %>% group_by(i, HH) %>% filter(n()==1) %>% mutate(susceptible=HHsize-inf)
  hassecondary <- final %>% group_by(i, HH) %>% filter(n()>1) %>% 
    mutate(secondaryrep=ifelse(secondaryinfs<=1, 1, secondaryinfs), 
           secondaryrep=c(1, secondaryrep[1:(n()-1)])) %>% 
    mutate(groupsecondary=rep(1:n(), times=secondaryrep)[1:n()]) %>% 
    select(!secondaryrep)
  
  susceptibles <- hassecondary %>% group_by(i, HH, groupsecondary) %>% summarise(HHsize=first(HHsize), infs=n()) %>% 
    mutate(susceptible=HHsize-cumsum(infs))
  
  hassecondary <- hassecondary %>% left_join(susceptibles) %>% select(!c(groupsecondary, infs))
  
  final <- rbind(nosecondary, hassecondary) %>% arrange(i, HH, time)
  final
}

d1 <- find_secondary_and_susceptible(inf1) %>% mutate(attackrate=secondaryinfs/susceptible) 
d2 <- find_secondary_and_susceptible(inf2) %>% mutate(attackrate=secondaryinfs/susceptible) 
d3 <- find_secondary_and_susceptible(inf3) %>% mutate(attackrate=secondaryinfs/susceptible) 
d4 <- find_secondary_and_susceptible(inf4) %>% mutate(attackrate=secondaryinfs/susceptible) 

attack_rate <- function(d) {
  (d %>% filter(susceptible>0) %>% group_by(i) %>% 
    summarise(attackrate=mean(attackrate,na.rm=T)))$attackrate %>% mean()
}

d1 %>% attack_rate()
d2 %>% attack_rate()
d3 %>% attack_rate()
d4 %>% attack_rate()

attack_rate_month <- function(d) {
  d %>% filter(susceptible>0) %>% group_by(i, month) %>% 
    summarise(attackrate=mean(attackrate,na.rm=T)) %>% 
    group_by(month) %>% summarise(attackrate=mean(attackrate))  
}

ggplot() + 
  geom_line(data = attack_rate_month(d1), aes(month, attackrate, color="Environmental transmission only")) + 
  geom_line(data = attack_rate_month(d2), aes(month, attackrate, color="75% household transmission")) + 
  geom_line(data = attack_rate_month(d3), aes(month, attackrate, color="50% household transmission")) + 
  geom_line(data = attack_rate_month(d4), aes(month, attackrate, color="25% household transmission")) + 
  scale_y_continuous("Attack rate", expand=c(0.01,0)) + 
  scale_x_continuous("Time (months)", breaks=1:13) + 
  scale_color_brewer(palette="Dark2", direction=-1, name=element_blank()) + 
  theme(panel.grid.minor = element_blank())


1-(1-0.0011)**((mean(3:6)-1)*7)  

