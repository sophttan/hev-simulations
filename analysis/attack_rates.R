# estimate within hh attack rates

rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/blended_model/cumulative_inc_30/"))

inf <- read_csv("simulated_data/25p75e_hhrisk5.csv")

find_secondary_and_susceptible <- function(inf) {
  inf <- inf %>% arrange(i, HH, time) %>% group_by(i, HH) %>% 
    mutate(otherinfs=list(time), inf=1:n()) %>% 
    rowwise() %>% 
    mutate(secondaryinfs = list(otherinfs[unlist(otherinfs)-time<=45 & unlist(otherinfs)-time > 7])) %>%
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

d <- find_secondary_and_susceptible(inf)
d <- d %>% mutate(attackrate=secondaryinfs/susceptible)
d %>% filter(susceptible>0) %>% group_by(i, HH) %>% summarise(attackrate=mean(attackrate,na.rm=T)) %>% 
  group_by(i) %>% summarise(attackrate=mean(attackrate,na.rm=T)) %>% summary()

