rm(list=ls())
gc()
library(tidyverse)
library(purrr)

initial_pop <- data.frame(No=1:1000, 
                          HH=rep(1:200, each=5), 
                          S=c(0,rep(1,999)), 
                          E=0, Ecounter=0, 
                          I=0, Icounter=0, 
                          R=0) # 200 households, 4 compartments (SEIR)
results <- initial_pop[,1:2] %>% mutate(Itype=NA) # number of infections per household over time

time <- 52*7*1.5 # num days to simulate (12 months total)

# incubation period and infectious period fixed
num_weeks_inc = 4*7
num_weeks_inf = 1*7

SEIR_environment <- function(d, res, b, inc, inf) {
  
  for (i in 1:time) {
    d$R[d$Icounter==inf] <- 1
    d$I[d$Icounter==inf] <- 0
    d$Icounter[d$Icounter==inf] <- 0
    
    d$I[d$Ecounter==inc] <- 1
    d$E[d$Ecounter==inc] <- 0
    d$Ecounter[d$Ecounter==inc] <- 0
    
    summary_data <- d %>% group_by(HH) %>% summarise(Sh=sum(S), Ih=sum(I)) %>% 
      mutate(Sc=sum(Sh)-Sh, Ic=sum(Ih)-Ih) %>% ungroup()
    res <- res %>% cbind(I=summary_data$Ih)
    
    risk <- b*d$S
    
    new_inf <- rbinom(1000, 1, risk)
    
    d$E[new_inf==1] <- 1
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
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

res <- rep(0, time)
beta <- 0.001
for (i in 1:100) {
  final <- SEIR_environment(initial_pop, results, beta, num_weeks_inc, num_weeks_inf)
  names(final) <- c("No", "HH", "Type", 1:time)
  f <- final %>% pivot_longer(cols = 4:ncol(.), names_to = "time") %>% mutate(time=as.numeric(time))
  f <- f %>% group_by(No, value) %>% summarise_all(first) %>% filter(value==1) 
  total_inf <- f %>% group_by(time) %>% summarise(inf=sum(value)) %>% full_join(expand.grid(time=1:time)) %>% replace_na(list(inf=0)) %>% arrange(time)
  res <- res+total_inf$inf
}
res <- res/100
res

average <- data.frame(time=1:time, res)
p<-average %>% ggplot(aes(time, res)) + geom_line() + 
  scale_x_continuous("Time (days)") + 
  scale_y_continuous("Incidence (number of new infections")
p

ggsave(p, filename="results/environmental.jpg")
write_csv(average, "results/environmental.csv")


# compare <- expand.grid(HH=1:200, time=1:364) 
# names(final) <- c("No", "HH", 1:364)
# names(final2) <- c("HH", 1:364)
# f2 <- final2 %>% pivot_longer(cols = 2:(time+1), names_to = "time") %>% mutate(time=as.numeric(time))
# 
# compare <- compare %>% left_join(f, by=c("HH", "time")) %>% left_join(f2, by=c("HH","time"))
# compare
# 
# compare %>% group_by(time) %>% summarise_all(sum) %>% select(!HH) %>% ggplot(aes(time)) + 
#   geom_line(aes(y=value.x), color="blue") + 
#   geom_line(aes(y=value.y), color="black") + 
#   scale_y_continuous("Number of infections") + 
#   scale_x_continuous("Weeks")
# 

averagehh <- read_csv("results/hh.csv")

p +
  geom_line(aes(y=averagehh$res), color="darkblue") 
