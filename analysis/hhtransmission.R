rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/separate_models/cumulative_inc_30/"))

time <- 52*7 # num days to simulate (12 months total)
months <- rep(1:13, each=30)[1:time]
days_months <- data.frame(day=1:time, month=months)

# incubation period and infectious period fixed
num_weeks_inc = 4*7
num_weeks_inf = 1*7

SEIR <- function(rr, b, inc, inf) {
  
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
                  S=c(0,rep(1,999)), 
                  E=c(1,rep(0,999)), Ecounter=c(28,rep(0,999)), 
                  I=0, Icounter=0, 
                  R=0, 
                  inc=c(28,rep(0,999)), inf=0) # variable number of households (size 3-6), total pop of 1000
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
    
    risk_hh <- (summary_data$S)*rr*b*(summary_data$Ih)/(summary_data$HHsize)
    risk_c <- (summary_data$S)*b*(summary_data$Ic)/(1000-summary_data$HHsize)
    
    new_inf_hh <- rbinom(nrow(d), 1, risk_hh)
    new_inf_c <- rbinom(nrow(d), 1, risk_c)
    new_exposed <- new_inf_hh==1|new_inf_c==1
    
    if(sum(new_exposed)>0) {
      d$E[new_exposed] <- 1
      d$inc[new_exposed] <- rnorm(sum(new_exposed, na.rm=T), mean=inc, sd=2) %>% round()
      case_type <- case_when(new_inf_hh==1&new_inf_c==1~"B",
                             new_inf_hh==1~"H", 
                             new_inf_c==1~"C")
      res$Itype <- ifelse(!is.na(res$Itype), res$Itype, case_type) 
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
    # if(i==210){
    #   view(d)
    # }
  }
#  print(d)
  return(res)
}


# final <- SEIR(initial_pop, results, 0.03, num_weeks_inc, num_weeks_inf)
# names(final) <- c("No", "HH", "Type", 1:time)
# f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
# f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) %>% arrange(time)
# f
# f %>%
#   group_by(time) %>% summarise(inf=sum(value)) %>% ggplot(aes(time, inf)) + geom_line()

rr <- 2
beta <- 0.12
sims <- 1000
num_hh <- rep(0, sims)

for (i in 1:sims) {
  final <- SEIR(rr, beta, num_weeks_inc, num_weeks_inf)
  num_hh[i] <- max(final$HH)
  
  # restructure table
  names(final) <- c("No", "HHsize", "HH", "Type", 1:time)
  f <- final %>% pivot_longer(cols = 5:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
  f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) 
  f <- f %>% group_by(HH) %>% mutate(day_limits = list(time)) %>% 
    ungroup() %>% rowwise() %>% 
    mutate(has_hh = any((time - unlist(day_limits)) < 45 & (time - unlist(day_limits)) > 7)) %>% 
    select(c(time, HHsize, HH, Type, has_hh))
  
  if(i==1){
    inf_type <- cbind(i=i, f)
  }else{
    inf_type <- rbind(inf_type, cbind(i=i, f))
  }
}

inf_type <- inf_type %>% left_join(days_months, by=c("time"="day"))

# average total infections
(inf_type %>% nrow())/sims

write_csv(inf_type, "simulated_data/hh_risk2.csv")



inf_type_overall <- inf_type %>% group_by(month) %>% summarise(count=n())


p <- inf_type %>% group_by(month, Type) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=Type, color=Type)) + 
  scale_color_discrete("Source of infection", labels=c("Both community and household", "Community", "Household", "Index case")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with known source of infection", 
       subtitle="Household relative risk = 2, Cumulative incidence ~ 30%")
p %>% ggsave(filename = "figures/obs_hh_risk2.jpg")


p <- inf_type %>% group_by(month, has_hh) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=has_hh, color=has_hh)) + 
  scale_color_discrete("Predicted source of infection", 
                       breaks=c(TRUE, FALSE),
                       labels=c("Household infection", "Not household infection")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with predicted source of infection", 
       subtitle="Household relative risk = 2, Cumulative incidence ~ 30%")
p %>% ggsave(filename = "figures/pred_hh_risk2.jpg")


inf_type <- inf_type %>% group_by(i, month)
inf_type_check <- inf_type %>% mutate(Observed=ifelse(!is.na(Type)&(Type=="H"|Type=="B"), "Household infection", "Not household"),
                                      Predicted=ifelse(has_hh, "Household infection", "Not household")) %>% 
  pivot_longer(cols = c("Observed", "Predicted")) %>% 
  group_by(month, name, value) %>% summarise(count=n()/sims)
p <- inf_type_check %>% ggplot(aes(month, count, color=name, shape=value, group=interaction(value, name))) + 
  geom_line() + 
  geom_point(size=2) + 
  scale_color_brewer(palette="Dark2", direction=-1, name=element_blank()) + 
  scale_shape(name="Source of infection") + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Number of household infections") 
p %>% ggsave(filename = "figures/hh_risk2_classification.jpg")


# compare observed and predicted fraction of cases with prior household infection
(inf_type %>% 
    mutate(has_hh_obs=ifelse(!is.na(Type)&(Type=="H"|Type=="B"), T, F)) %>% 
    group_by(i) %>% summarise(prop_obs=mean(has_hh_obs), prop_pred=mean(has_hh))) %>% summary()

mean((inf_type %>% group_by(i) %>% summarise(HH=length(unique(HH)))/num_hh)$HH)


# calculate average hh attack rate across hh and simulations
inf_type %>% group_by(i, HH) %>% summarise(attack_rate_hh=n()/5) %>% 
  group_by(i) %>% summarise(attack_rate=mean(attack_rate_hh)) %>% 
  summary()

