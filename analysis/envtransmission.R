rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/separate_models/cumulative_inc_5/"))

time <- 52*7 # num days to simulate (12 months total)
months <- rep(1:13, each=30)[1:time]
days_months <- data.frame(day=1:time, month=months)

# incubation period and infectious period fixed
num_weeks_inc = 4*7
num_weeks_inf = 1*7

SEIR_environment <- function(b, inc, inf) {
  
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
    
    res <- res %>% cbind(I=d$I)
    
    risk <- b*d$S
    
    new_exposed <- rbinom(nrow(d), 1, risk)

    if(sum(new_exposed)>0) {
      d$E[new_exposed==1] <- 1
      d$inc[new_exposed==1] <- rnorm(sum(new_exposed, na.rm=T), mean=inc, sd=2) %>% round()
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
  }
  
  return(res)
}


beta <- 0.00016
sims <- 1000
num_hh <- rep(0, sims)

for (i in 1:sims) {
  final <- SEIR_environment(beta, num_weeks_inc, num_weeks_inf)
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
    if(nrow(f)==0) {next}
    inf_type <- cbind(i=i, f)
  }else{
    if(nrow(f)==0) {next}
    inf_type <- rbind(inf_type, cbind(i=i, f))
  }
}

inf_type <- inf_type %>% left_join(days_months, by=c("time"="day"))

# average total infections
(inf_type %>% nrow())/sims

write_csv(inf_type, "simulated_data/environmental.csv")

p <- inf_type %>% group_by(month, has_hh) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=has_hh, color=has_hh)) + 
  scale_color_discrete("Predicted source of infection", 
                       breaks=c(TRUE, FALSE),
                       labels=c("Household infection", "Not household infection")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with predicted source of infection", 
       subtitle="Only environmental source transmission, Cumulative incidence ~ 5%")
p %>% ggsave(filename = "figures/pred_environment.jpg")


mean((inf_type %>% group_by(i) %>% summarise(HH=length(unique(HH)))/num_hh)$HH)

(inf_type %>% group_by(i) %>% 
    summarise(prop_pred=mean(has_hh))) %>% summary()


# inf_type <- inf_type %>% group_by(i, month)
# inf_type_check <- inf_type %>% summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
#   group_by(month) %>% summarise_all(mean)
# p <- inf_type_check %>% ggplot(aes(month)) + 
#   geom_line(aes(y=has_hh)) + 
#   scale_x_continuous("Time (months)") +
#   scale_y_continuous("Proportion of household infections") +
#   theme(legend.title = element_blank())
# p %>% ggsave(filename = "results/environmental_incidence/environmental_classification.jpg")



