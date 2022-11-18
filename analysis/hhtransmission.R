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
    
    risk_hh <- summary_data$S*4*b*summary_data$Ih/1000
    risk_c <- summary_data$S*b*summary_data$Ic/1000
    
    # print(risk_hh[1:10])
    # print(risk_c[1:10])

    new_inf_hh <- rbinom(nrow(d), 1, risk_hh)
    new_inf_c <- rbinom(nrow(d), 1, risk_c)
    
    # print(new_inf_hh[1:10])
    # print(new_inf_c[1:10])
    
    d$E[new_inf_hh==1|new_inf_c==1] <- 1
    case_type <- case_when(new_inf_hh==1&new_inf_c==1~"B",
                           new_inf_hh==1~"H", 
                           new_inf_c==1~"C")
    res$Itype <- ifelse(!is.na(res$Itype), res$Itype, case_type)
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
  }
  return(res)
}

 
# final <- SEIR(initial_pop, results, 0.35, num_weeks_inc, num_weeks_inf)
# names(final) <- c("No", "HH", "Type", 1:time)
# f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
# f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) %>% arrange(time) 
# f
# f %>% 
#   group_by(time) %>% summarise(inf=sum(value)) %>% ggplot(aes(time, inf)) + geom_line()

hh <- 0
community <- 0
res <- rep(0, time)
beta <- 0.35
for (i in 1:100) {
  final <- SEIR(initial_pop, results, beta, num_weeks_inc, num_weeks_inf)
  names(final) <- c("No", "HH", "Type", 1:time)
  f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
  f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) 
  total_inf <- f %>% group_by(time) %>% summarise(inf=sum(value)) %>% full_join(expand.grid(time=1:time)) %>% replace_na(list(inf=0)) %>% arrange(time)
  res <- res+total_inf$inf
  hh <- hh+sum(f$Type=="H", na.rm=T)
  community <- community+sum(f$Type=="C", na.rm=T)
}
res <- res/100
hh <- hh/100
community <- community/100
res

average <- data.frame(time=1:time, res)
p <- average %>% ggplot(aes(time, res)) + geom_line() + 
  scale_x_continuous("Time (days)") + 
  scale_y_continuous("Incidence (number of new infections")

p
ggsave(p, filename="results/hh.jpg")
write_csv(average, "results/hh.csv")

hh
community

# 800-850 total cases