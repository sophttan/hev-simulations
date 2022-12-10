rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/blended_model/cumulative_inc_10/"))

time <- 52*7 # num days to simulate (12 months total)
months <- rep(1:13, each=30)[1:time]
days_months <- data.frame(day=1:time, month=months)

# incubation period and infectious period fixed
num_weeks_inc = 4*7
num_weeks_inf = 1*7

SEIR_blend <- function(d, res, rr, b_hh, b_e, inc, inf) {
  
  # vary household size
  # randomly sample household sizes such that total population is 1000 individuals
  hh_size <- sample(x = c(3, 4, 5, 6), size=340, replace=T)
  hh_size <- hh_size[which(cumsum(hh_size) < 1000)]
  leftover <- 1000-sum(hh_size)
  if(leftover < 3) {
    hh <- 1:length(hh_size)
    sampled <- sample(hh[hh_size<6], leftover)
    hh_size[sampled] <- hh_size[sampled] + 1
  } else {
    hh_size <- c(hh_size, leftover)
  }
  
  d <- data.frame(No=1:1000, 
                  HHsize=rep(hh_size, times=hh_size),
                  HH=rep(1:length(hh_size), times=hh_size), 
                  S=1, 
                  E=0, Ecounter=0, 
                  I=0, Icounter=0, 
                  R=0, 
                  inc=0, inf=0) # variable number of households (size 3-6), total pop of 1000
  res <- d[,1:3] %>% mutate(Itype=NA)
  
  for (i in 1:time) {
    recovered <- d$inf>0 & d$Icounter==d$inf
    if(sum(recovered,na.rm=T)>0) {
      d$R[recovered] <- 1
      d$I[recovered] <- 0
      d$Icounter[recovered] <- 0 
    }
    
    new_inf <- d$inc>0 &d$Ecounter==d$inc
    if(sum(new_inf,na.rm=T)>0) {
      random_inf <- rnorm(sum(new_inf, na.rm=T), mean=inf, sd=1) %>% round()
      d$I[new_inf] <- 1
      d$inf[new_inf] <- random_inf
      d$E[new_inf] <- 0
      d$Ecounter[new_inf] <- 0 
    }
    
    summary_data <- d %>% group_by(HH) %>% 
      mutate(Ih=sum(I)) %>% ungroup() %>% mutate(Ic=sum(I)-Ih) 
    res <- res %>% cbind(I=summary_data$I)
    
    risk_hh <- d$S*rr*b_hh*summary_data$Ih/4
    risk_c <- d$S*b_hh*summary_data$Ic/995
    risk_e <- b_e*d$S
    
    new_inf <- rbinom(1000, 1, risk_e)
    new_inf_hh <- rbinom(1000, 1, risk_hh)
    new_inf_c <- rbinom(1000, 1, risk_c)
    new_exposed <- new_inf==1|new_inf_hh==1|new_inf_c==1
    
    if(sum(new_exposed)>0) {
      d$E[new_exposed] <- 1
      d$inc[new_exposed] <- rnorm(sum(new_exposed, na.rm=T), mean=inc, sd=2) %>% round()
      case_type <- case_when(new_inf==1~"E",
                             new_inf==1&(new_inf_hh==1|new_inf_c==1)~"B",
                             new_inf_hh==1~"H", 
                             new_inf_c==1~"C")
      res$Itype <- ifelse(!is.na(res$Itype), res$Itype, case_type) 
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
  }
  return(res)
}

rr <- 5
beta_hh <- 0.031
beta_env <- 0.00008
sims <- 1000
num_hh <- rep(0, sims)

for (i in 1:sims) {
  final <- SEIR_blend(initial_pop, results, rr, beta_hh, beta_env, num_weeks_inc, num_weeks_inf)
  num_hh[i] <- max(final$HH)
  
  # restructure table
  names(final) <- c("No", "HHsize", "HH", "Type", 1:time)
  f <- final %>% pivot_longer(cols = 5:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
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

(inf_type %>% nrow())/sims
inf_type %>% group_by(Type) %>% summarise(n=n()/nrow(.)) 

write_csv(inf_type, "simulated_data/75p25e_hhrisk5.csv")

# compare observed and predicted fraction of cases with prior household infection
(inf_type %>% 
    mutate(has_hh_obs=ifelse(Type=="H", T, F)) %>% 
    group_by(i) %>% summarise(prop_obs=mean(has_hh_obs), prop_pred=mean(has_hh))) %>% summary()

mean((inf_type %>% group_by(i) %>% summarise(HH=length(unique(HH)))/num_hh)$HH)

# calculate average hh attack rate across hh and simulations
inf_type %>% group_by(i, HH) %>% summarise(attack_rate_hh=n()/5) %>% 
  group_by(i) %>% summarise(attack_rate=mean(attack_rate_hh)) %>% 
  summary()

p <- inf_type %>% group_by(month, Type) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=Type, color=Type)) + 
  scale_color_discrete("Source of infection", 
                       breaks=c("H", "C", "E"),
                       labels=c("Household", "Community", "Environment")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with known source of infection", 
       subtitle="Household relative risk = 5\n75:25 ratio of transmission from person-person contact:environment\nCumulative incidence ~ 10%")
p
p %>% ggsave(filename = "figures/75p25e_hhrisk5_obs.jpg")

p <- inf_type %>% group_by(month, has_hh) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=has_hh, color=has_hh)) + 
  scale_color_discrete("Predicted source of infection", 
                       breaks=c(TRUE, FALSE),
                       labels=c("Household infection", "Not household infection")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with predicted source of infection", 
       subtitle="Household relative risk = 5\n75:25 ratio of transmission from person-person contact:environment\nCumulative incidence ~ 10%")
p
p %>% ggsave(filename = "figures/75p25e_hhrisk5_pred.jpg")
