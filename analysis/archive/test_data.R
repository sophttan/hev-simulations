rm(list=ls())
gc()
library(tidyverse)
library(purrr)

time <- 60 # num days to simulate (12 months total)
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
                  I=0, Icounter=0, Itype=NA,
                  R=0, inc=c(28,rep(0,999)), inf=0) # variable number of households (size 3-6), total pop of 1000
  res <- NULL
  
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
    } else {
      case_type <- c("NA")
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
    summary <- d %>% select(S,E,I,R) %>% summarise_all(sum) %>% 
      mutate(exposed = sum(new_exposed),
             incidence = sum(new_inf), 
             hh_inf = sum(case_type=="H",na.rm=T), 
             community_inf = sum(case_type=="C",na.rm=T), 
             both_hh_c = sum(case_type=="B",na.rm=T))
    res <- res %>% rbind(summary)
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

final <- SEIR(rr, beta, num_weeks_inc, num_weeks_inf)

final %>% write_csv("/Users/sophiatan/Dropbox/ID Modeling/Transmission Methods Modeling/sample_data_inc30_rr2.csv")
