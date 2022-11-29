rm(list=ls())
gc()
library(tidyverse)
library(purrr)

setwd(here::here("results/separate_models/cumulative_inc_5/"))

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

SEIR_environment <- function(d, res, b, inc, inf) {
  
  for (i in 1:time) {
    d$R[d$Icounter==inf] <- 1
    d$I[d$Icounter==inf] <- 0
    d$Icounter[d$Icounter==inf] <- 0
    
    d$I[d$Ecounter==inc] <- 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$E[d$Ecounter==inc] <- 0
    d$Ecounter[d$Ecounter==inc] <- 0
    
    summary_data <- d %>% group_by(HH) %>% 
      mutate(Ih=sum(I)) %>% ungroup() %>% mutate(Ic=sum(I)-Ih) 
    res <- res %>% cbind(I=summary_data$I)
    
    risk <- b*d$S
    
    new_inf <- rbinom(1000, 1, risk)
    
    d$E[new_inf==1] <- 1
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$S[d$E==1] <- 0
    
  }
  
  return(res)
}


beta <- 0.00016

for (i in 1:100) {
  final <- SEIR_environment(initial_pop, results, beta, num_weeks_inc, num_weeks_inf)
  
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

# average total infections
(inf_type %>% nrow())/100

write_csv(inf_type, "simulated_data/environmental.csv")

p <- inf_type %>% group_by(month, has_hh) %>% summarise(count=n()/100) %>%
  ggplot() + geom_line(aes(month, count, group=has_hh, color=has_hh)) + 
  scale_color_discrete("Predicted source of infection", 
                       breaks=c(TRUE, FALSE),
                       labels=c("Household infection", "Not household infection")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with predicted source of infection", 
       subtitle="Only environmental source transmission, Cumulative incidence ~ 10%")
p %>% ggsave(filename = "figures/pred_environment.jpg")


(inf_type %>% group_by(i) %>% summarise(prop=mean(has_hh)))$prop %>% mean()

# calculate average hh attack rate across hh and simulations
inf_type %>% group_by(i, HH) %>% summarise(attack_rate_hh=n()/5) %>% 
  group_by(i) %>% summarise(attack_rate=mean(attack_rate_hh)) %>% 
  summary()

# inf_type <- inf_type %>% group_by(i, month)
# inf_type_check <- inf_type %>% summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
#   group_by(month) %>% summarise_all(mean)
# p <- inf_type_check %>% ggplot(aes(month)) + 
#   geom_line(aes(y=has_hh)) + 
#   scale_x_continuous("Time (months)") +
#   scale_y_continuous("Proportion of household infections") +
#   theme(legend.title = element_blank())
# p %>% ggsave(filename = "results/environmental_incidence/environmental_classification.jpg")



