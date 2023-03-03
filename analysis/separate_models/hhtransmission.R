rm(list=ls())
gc()
library(tidyverse)
library(purrr)
library(OOR)

setwd(here::here("results/separate_models/cumulative_inc_10/"))

#### Model schematics #### 
time <- 365 # num days to simulate (12 months total)
months <- rep(1:13, each=30)[1:time] # months roughly 30 days long
days_months <- data.frame(day=1:time, month=months) # use months for plotting

# incubation period and infectious period fixed (days)
num_weeks_inc = 28 
num_weeks_inf = 7

create_hh <- function() {
  # vary household size
  # randomly sample household sizes such that total population is 1000 individuals
  hh_size <- sample(x = c(3, 4, 5, 6), size=340, replace=T)
  # keep households such that total population is < 1000
  hh_size <- hh_size[which(cumsum(hh_size) < 1000)]
  
  # assign last household so it has 3-6 people
  # if last household has <3 people, randomly select other households that have <6
  # people and add 1 person so total is 1000
  leftover <- 1000-sum(hh_size)
  if(leftover < 3) {
    hh <- 1:length(hh_size)
    sampled <- sample(hh[hh_size<6], leftover)
    hh_size[sampled] <- hh_size[sampled] + 1
  } else {
    hh_size <- c(hh_size, leftover)
  }
  return(hh_size)
}

#### SEIR model ####
SEIR <- function(bh, bc, inc, inf) {
  #### SET UP POPULATION ####
  # initialize random households
  hh_size <- create_hh()
  
  # initialize population data frame
  # initially all individuals in population are susceptible except for 1 exposed person
  # No = individual id
  # HHsize = household size
  # HH = household id
  # S, E, I, R = binary variables that reflect individual states (susceptible, exposed, infectious, recovered)
  # Ecounter, Icounter = number of days spent in exposed and infectious states
  # inc = incubation period - randomly drawn from normal distribution around 28 days (assigned after exposure)
  # inf = infectious period - randomly drawn from normal distribution around 7 days (assigned after infection)
  d <- data.frame(No=1:1000, 
                  HHsize=rep(hh_size, times=hh_size),
                  HH=rep(1:length(hh_size), times=hh_size), 
                  S=c(0,rep(1,999)), 
                  E=c(1,rep(0,999)), Ecounter=c(1,rep(0,999)), 
                  I=0, Icounter=0, 
                  R=0, 
                  inc=c(round(rnorm(1, mean=inc, sd=2)),rep(0,999)), inf=0) 
  
  # results data frame has No, HHsize, HH, Type of infection (household/community), time is first day of infection
  # assumes perfect ascertainment at first day of infection period
  res <- d[,1:3] %>% mutate(Type=NA, time=NA)
  
  #### SIMULATION ####
  for (i in 1:time) {
    # identify individuals that have recovered
    recovered <- d$inf>0 & d$Icounter==d$inf
    if(sum(recovered,na.rm=T)>0) {
      d$R[recovered] <- 1
      d$I[recovered] <- 0
      d$Icounter[recovered] <- 0 
    }
    
    # identify individuals that are newly infectious
    new_inf <- d$inc>0 &d$Ecounter==d$inc
    if(sum(new_inf,na.rm=T)>0) {
      random_inf <- rnorm(sum(new_inf, na.rm=T), mean=inf, sd=1) %>% round()
      d$I[new_inf] <- 1
      d$inf[new_inf] <- random_inf
      d$E[new_inf] <- 0
      d$Ecounter[new_inf] <- 0 
      
      res$time[new_inf==1] <- i
    }
    
    # find total number of infectious individuals within and outside of household
    summary_data <- d %>% group_by(HH) %>% 
      mutate(Ih=sum(I)) %>% ungroup() %>% mutate(Ic=sum(I)-Ih) 

    # estimate risk of infection
    risk_hh <- (summary_data$S)*bh*(summary_data$Ih)/1000
    risk_c <- (summary_data$S)*bc*(summary_data$Ic)/1000
    
    # infection occurs as a bernoulli event for each individual based on risk of infection
    new_inf_hh <- rbinom(nrow(d), 1, risk_hh)
    new_inf_c <- rbinom(nrow(d), 1, risk_c)
    new_exposed <- new_inf_hh==1|new_inf_c==1
    
    # identify individuals with new exposure
    if(sum(new_exposed)>0) {
      d$E[new_exposed] <- 1
      d$inc[new_exposed] <- rnorm(sum(new_exposed, na.rm=T), mean=inc, sd=2) %>% round()
      case_type <- case_when(new_inf_hh==1&new_inf_c==1~"B",
                             new_inf_hh==1~"H", 
                             new_inf_c==1~"C")
      res <- res %>% mutate(Type = ifelse(!is.na(Type), Type, case_type))
    }
    
    d$Ecounter[d$E==1] <- d$Ecounter[d$E==1] + 1
    d$Icounter[d$I==1] <- d$Icounter[d$I==1] + 1
    d$S[d$E==1] <- 0
    
  }
  return(res)
}


#### Calibrating parameters - optimizing function #### 
# dist_to_target <- function(params,target){
#   # Function to compute sum of square scaled differences between simulation output
#   # and targets
#   # Targets are cumulative incidence (5%, 10%, 30%) and proportion of household infections (25%, 50%, 75%)
#   
#   betah <- params[1]
#   betac <- params[2]
# 
#   state <- SEIR(bh=betah, bc=betac, num_weeks_inc, num_weeks_inf)
#   state <- state %>% filter(!is.na(time))
# 
#   incidence <- nrow(state)/1000*100
#   
#   prop_types <- state %>% filter(!Type %>% is.na()) %>% group_by(Type) %>% summarise(prop=n()/nrow(.))
#   prop_hh <- (prop_types$prop[prop_types$Type=="H"]) * 100
# 
#   outvector <- (c(incidence, prop_hh)-target)/target
#   
#   return(sum((outvector)**2)*100)
# }


#### Calibration ####
# goal <- c(10, 25)
# fit <- StoSOO(c(NA, NA),
#               function(x){return(dist_to_target(x,goal))},
#               lower=c(0,0),upper=c(30,1),nb_iter=1000,
#               control = list(verbose = 0, type = "sto", max = FALSE, light = FALSE))
# cat(c("Params:",fit$par,"\n"))
# cat(c("Value:",fit$value,"\n"))

sims <- 10
num_hh <- rep(0, sims)
inc <- rep(0, sims)
prop_hh <- rep(0, sims)

# test out fitted parameters
betah <- 26.210651842940592
betac <- 0.11921066815772194

inf_type <- NULL

for (i in 1:sims) {
  final <- SEIR(betah, betac, num_weeks_inc, num_weeks_inf)
  num_hh[i] <- max(final$HH)
  
  # restructure table
  f <- final %>% filter(!is.na(time))
  f <- f %>% group_by(HH) %>% mutate(day_limits = list(time)) %>%
    ungroup() %>% rowwise() %>%
    mutate(has_hh = any((time - unlist(day_limits)) < 45 & (time - unlist(day_limits)) > 7)) %>%
    select(c(time, HHsize, HH, Type, has_hh))
  
  if(nrow(f)==1) {
    prop_hh[i] <- NA
  }else{
    prop_hh[i] <- sum(f$Type=="H"|f$Type=="B",na.rm=T)/nrow(f)
  }
  
  inf_type <- rbind(inf_type, cbind(i=i, f))
  inc[i] <- nrow(f)/1000
}

inf_type %>% head()
res <- data.frame(i=1:sims, inc=inc, prop_hh=prop_hh)

# average infections
res$inc %>% summary()

# average fraction household transmission
res$prop_hh %>% mean(na.rm=T)

library(patchwork)
p1 <- res %>% ggplot(aes(inc)) + 
  geom_histogram() + 
  scale_x_continuous(expand=c(0.05,0)) + 
  scale_y_continuous(expand=c(0,0))
p2 <- res %>% ggplot(aes(prop_hh)) + 
  geom_histogram() + 
  scale_x_continuous(expand=c(0.05,0), limits=c(0,1)) + 
  scale_y_continuous(expand=c(0,0))
p1|p2

inf_type <- inf_type %>% left_join(days_months, by=c("time"="day"))

# save dataset
write_csv(inf_type, "simulated_data/hh_25h75c.csv")

#### Plot results ####
inf_type_overall <- inf_type %>% group_by(month) %>% summarise(count=n())

# plot incidence over time stratified by observed source of infection
p <- inf_type %>% group_by(month, Type) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=Type, color=Type)) + 
  scale_color_discrete("Source of infection", labels=c("Both community and household", "Community", "Household", "Index case")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with known source of infection", 
       subtitle="25:75 Household:community transmission, Cumulative incidence ~ 30%")
p %>% ggsave(filename = "figures/obs_hh_25h75c.jpg")


# plot incidence over time stratified by predicted source of infection
p <- inf_type %>% group_by(month, has_hh) %>% summarise(count=n()/sims) %>%
  ggplot() + geom_line(aes(month, count, group=has_hh, color=has_hh)) + 
  scale_color_discrete("Predicted source of infection", 
                       breaks=c(TRUE, FALSE),
                       labels=c("Household infection", "Not household infection")) + 
  scale_x_continuous("Time (months)") + 
  scale_y_continuous("Incidence (number of new infections)") + 
  labs(title="Incidence over time with predicted source of infection", 
       subtitle="25:75 Household:community transmission, Cumulative incidence ~ 30%")
p %>% ggsave(filename = "figures/pred_hh_25h75c.jpg")


# plot comparison between observed and predicted cases
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
p %>% ggsave(filename = "figures/hh_25h75c_classification.jpg")


# compare observed and predicted fraction of cases with prior household infection
(inf_type %>% 
    mutate(has_hh_obs=ifelse(!is.na(Type)&(Type=="H"|Type=="B"), T, F)) %>% 
    group_by(i) %>% summarise(prop_obs=mean(has_hh_obs), prop_pred=mean(has_hh))) %>% summary()

# proportion of households with any infection
mean((inf_type %>% group_by(i) %>% summarise(HH=length(unique(HH)))/num_hh)$HH)



