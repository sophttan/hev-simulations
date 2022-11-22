rm(list=ls())
gc()
library(tidyverse)
library(purrr)

initial_pop <- data.frame(No=1:1000, 
                          HH=rep(1:200, each=5), 
                          S=c(0,rep(1,999)), 
                          E=c(1,rep(0,999)), Ecounter=c(28,rep(0,999)), 
                          I=0, Icounter=0, 
                          R=0) # 200 households, 4 compartments (SEIR)
results <- initial_pop[,1:2] %>% mutate(Itype=NA) # number of infections per household over time

time <- 52*7*1.5 # num days to simulate (12 months total)
months <- rep(1:19, each=30)[1:time]
days_months <- data.frame(day=1:time, month=months)

# incubation period and infectious period fixed
num_weeks_inc = 4*7
num_weeks_inf = 1*7

SEIR <- function(d, res, b, inc, inf) {
  
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
    
    risk_hh <- summary_data$S*2*b*summary_data$Ih/4
    risk_c <- summary_data$S*b*summary_data$Ic/995
    
    new_inf_hh <- rbinom(nrow(d), 1, risk_hh)
    new_inf_c <- rbinom(nrow(d), 1, risk_c)
    
    d$E[new_inf_hh==1|new_inf_c==1] <- 1
    case_type <- case_when(new_inf_hh==1&new_inf_c==1~"B",
                           new_inf_hh==1~"H", 
                           new_inf_c==1~"C")
    res$Itype <- ifelse(!is.na(res$Itype), res$Itype, case_type)
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
    # if(i==210){
    #   view(d)
    # }
  }
  return(res)
}

 
# final <- SEIR(initial_pop, results, 0.03, num_weeks_inc, num_weeks_inf)
# names(final) <- c("No", "HH", "Type", 1:time)
# f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
# f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) %>% arrange(time)
# f
# f %>%
#   group_by(time) %>% summarise(inf=sum(value)) %>% ggplot(aes(time, inf)) + geom_line()

res <- rep(0, time)
beta <- 0.09
for (i in 1:100) {
  final <- SEIR(initial_pop, results, beta, num_weeks_inc, num_weeks_inf)
  
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

# average total infections
(inf_type %>% nrow())/100

write_csv(inf_type, "results/hh_risk2.csv")

inf_type_overall <- inf_type %>% group_by(month) %>% summarise(count=n())

p <- inf_type %>% group_by(month, Type) %>% summarise(count=n()/100) %>%
  ggplot() + geom_line(aes(month, count, group=Type, color=Type)) + 
  scale_color_discrete("Source of infection", labels=c("Both community and household", "Community", "Household", "Index case")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") 
p %>% ggsave(filename = "results/household_incidence/hh_risk2.jpg")


inf_type <- inf_type %>% group_by(i, month)
inf_type_check <- inf_type %>% summarise(hh = sum(Type=="H",na.rm=T)/n(), has_hh = mean(has_hh,na.rm=T)) %>% 
  group_by(month) %>% summarise_all(mean)
p <- inf_type_check %>% ggplot(aes(month)) + 
  geom_line(aes(y=hh, color="blue")) + 
  geom_line(aes(y=has_hh, color="black")) + 
  scale_color_brewer(palette="Dark2", direction=-1, 
                     labels=c("Household infections (observed)", "Household infections (predicted)")) + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Proportion of household infections") +
  theme(legend.title = element_blank())
p %>% ggsave(filename = "results/household_classification/hh_risk2_classification.jpg")


inf_type %>% mutate(correct=(has_hh&(Type=="H"|Type=="B"))|(!has_hh&(Type%>%is.na()|Type=="C")))
