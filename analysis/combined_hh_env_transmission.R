rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/cumulative_inc_30/"))

initial_pop <- data.frame(No=1:1000, 
                          HH=rep(1:200, each=5), 
                          S=1, 
                          E=0, Ecounter=0, 
                          I=0, Icounter=0, 
                          R=0) # 200 households, 4 compartments (SEIR)
results <- initial_pop[,1:2] %>% mutate(Itype=NA) # number of infections per household over time

time <- 52*7 # num days to simulate (12 months total)
months <- rep(1:13, each=30)[1:time]
days_months <- data.frame(day=1:time, month=months)

# incubation period and infectious period fixed
num_weeks_inc = 4*7
num_weeks_inf = 1*7

SEIR_blend <- function(d, res, rr, b_hh, b_e, inc, inf) {
  
  total_inf <- data.frame(HH=1:200, I=0)
  
  for (i in 1:time) {
    d$R[d$Icounter==inf] <- 1
    d$I[d$Icounter==inf] <- 0
    d$Icounter[d$Icounter==inf] <- 0
    
    d$I[d$Ecounter==inc] <- 1
    d$E[d$Ecounter==inc] <- 0
    d$Ecounter[d$Ecounter==inc] <- 0
    
    summary_data <- d %>% group_by(HH) %>% 
      mutate(Ih=sum(I)) %>% ungroup() %>% mutate(Ic=sum(I)-Ih) 
    res <- res %>% cbind(I=summary_data$I)
    
    risk_hh <- d$S*rr*b_hh*summary_data$Ih/4
    risk_c <- d$S*b_hh*summary_data$Ic/995
    risk_e <- b_e*d$S
    
    new_inf <- rbinom(1000, 1, risk_e)
    new_inf_hh <- rbinom(1000, 1, risk_hh)
    new_inf_c <- rbinom(1000, 1, risk_c)
    
    d$E[new_inf_hh==1|new_inf_c==1|new_inf==1] <- 1
    case_type <- case_when(new_inf==1~"E",
                           new_inf==1&(new_inf_hh==1|new_inf_c==1)~"B",
                           new_inf_hh==1~"H", 
                           new_inf_c==1~"C")
    
    res$Itype <- ifelse(!is.na(res$Itype), res$Itype, case_type)
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
  }
  return(res)
}

rr <- 1
beta_hh <- 0.08
beta_env <- 0.0003

for (i in 1:100) {
  final <- SEIR_blend(initial_pop, results, rr, beta_hh, beta_env, num_weeks_inc, num_weeks_inf)
  
  # restructure table
  names(final) <- c("No", "HH", "Type", 1:time)
  f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
  f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) 
  f <- f %>% group_by(HH) %>% mutate(day_limits = list(time)) %>% 
    ungroup() %>% rowwise() %>% 
    mutate(has_hh = any((time - unlist(day_limits)) < 45 & (time - unlist(day_limits)) > 7)) %>% 
    select(c(time, HH, Type, has_hh))
  
  if(i==1){
    inf_type <- cbind(i=i, f)
  }else{
    inf_type <- rbind(inf_type, cbind(i=i, f))
  }
}

inf_type <- inf_type %>% left_join(days_months, by=c("time"="day"))

(inf_type %>% nrow())/100
inf_type %>% group_by(Type) %>% summarise(n=n()/nrow(.)) 

p <- inf_type %>% group_by(month, Type) %>% summarise(count=n()/100) %>%
  ggplot() + geom_line(aes(month, count, group=Type, color=Type)) + 
  scale_color_discrete("Source of infection", 
                       breaks=c("H", "C", "E"),
                       labels=c("Household", "Community", "Environment")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with known source of infection", 
       subtitle="Household relative risk = 1, Cumulative incidence ~ 30%")
p

