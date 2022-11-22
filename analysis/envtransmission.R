rm(list=ls())
gc()
library(tidyverse)
library(purrr)

initial_pop <- data.frame(No=1:1000, 
                          HH=rep(1:200, each=5), 
                          S=1, 
                          E=0, Ecounter=0, 
                          I=0, Icounter=0, 
                          R=0) # 200 households, 4 compartments (SEIR)
results <- initial_pop[,1:2] %>% mutate(Itype=NA) # number of infections per household over time

time <- 52*7*1.5 # num days to simulate (12 months total)
months <- rep(1:19, each=30)[1:time]
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

# final <- SEIR_environment(initial_pop, results, 0.001, num_weeks_inc, num_weeks_inf)
# names(final) <- c("No", "HH", "Type", 1:time)
# f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
# f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) %>% arrange(time) %>% select(!Type)
# f %>% 
#   group_by(time) %>% summarise(inf=sum(value)) %>% ggplot(aes(time, inf)) + geom_line()

beta <- .00075

for (i in 1:100) {
  final <- SEIR_environment(initial_pop, results, beta, num_weeks_inc, num_weeks_inf)
  
  # restructure table
  names(final) <- c("No", "HH", "Type", 1:time)
  f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
  f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) 
  f <- f %>% group_by(HH) %>% mutate(day_limits = list(time)) %>% 
    ungroup() %>% rowwise() %>% 
    mutate(has_hh = any(unlist(day_limits) >= (time-30) & unlist(day_limits) < time)) %>% 
    select(c(time, Type, has_hh))
  
  if(i==1){
    inf_type <- cbind(i=i, f)
  }else{
    inf_type <- rbind(inf_type, cbind(i=i, f))
  }
}

inf_type <- inf_type %>% left_join(days_months, by=c("time"="day"))

write_csv(inf_type, "results/simulated_data/environmental.csv")

inf_type_overall <- inf_type %>% group_by(month) %>% summarise(count=n()/100)

p <- inf_type_overall %>%
  ggplot() + geom_line(aes(month, count)) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") 
p %>% ggsave(filename = "results/environmental_incidence/environmental.jpg")

inf_type <- inf_type %>% group_by(i, month)
inf_type_check <- inf_type %>% summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
  group_by(month) %>% summarise_all(mean)
p <- inf_type_check %>% ggplot(aes(month)) + 
  geom_line(aes(y=has_hh)) + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Proportion of household infections") +
  theme(legend.title = element_blank())
p %>% ggsave(filename = "results/environmental_incidence/environmental_classification.jpg")



