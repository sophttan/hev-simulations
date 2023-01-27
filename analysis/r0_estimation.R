# estimate within hh r estimation

rm(list=ls())
gc()
library(tidyverse)
library(purrr)
library(EpiEstim)

setwd(here::here("results/separate_models/cumulative_inc_10/"))

inf1 <- read_csv("simulated_data/environmental.csv")
inf2 <- read_csv("simulated_data/hh_75h25c.csv")
inf3 <- read_csv("simulated_data/hh_50h50c.csv")
inf4 <- read_csv("simulated_data/hh_25h75c.csv")

inf1
?wallinga_teunis


t <- inf1 %>% group_by(i, HH, time) %>% 
  arrange(HH, time) %>% summarise(incid=n()) %>% ungroup()

R0 <- NULL
for (j in unique(t$i)) {
  data <- t %>% filter(i==j) 
  
  for (h in unique(data$HH)) {
    data2 <- data %>% filter(HH==h)

    full_data2 <- data2 %>%
      full_join(expand.grid(i=j, HH=h, time=1:364)) %>%
      arrange(time) %>%
      replace_na(list(incid=0))
    
    if(nrow(data2)==1) {
      if(R0 %>% is.null()){next}
      R0 <- R0 %>% rbind(c(data2$time+7, data2$time+45, rep(0, ncol(R0)-4), j, h))
      next
    }
    
    t_range <- min(data2$time):max(data2$time)
    t_start <- t_range + 7
    t_start <- t_start[t_start<(max(data2$time)-7)]
    t_end <- t_range + 45
    t_end <- t_end[t_end<=max(data2$time)]
    
    if (length(t_end) < length(t_start)) {
      t_end <- c(t_end, rep(max(data$time), length(t_start) - length(t_end)))
    } 
    
    if (length(t_start) == 0 & length(t_end) == 0) {
      R0 <- R0 %>% rbind(c(data2$time+7, data2$time+45, rep(0, ncol(R0)-4), j, h))
      next
    }

    res <- estimate_R(full_data2$incid,
                      method="parametric_si", 
                      config = list(t_start=t_start,
                                    t_end=t_end,
                                    mean_si = 28,
                                    std_si= 2))
    
    R0 <- R0 %>% rbind(res$R %>% mutate(i=j, HH=h))
  }
}

write_csv(R0, here::here("results/separate_models/cumulative_inc_10/r0/r0_hh_75h25c.csv"))
R0_summ <- R0 %>% group_by(i, t_start) %>% summarise(R0=mean(`Mean(R)`,na.rm=T)) %>% 
  group_by(t_start) %>% summarise(R0=mean(R0, na.rm=T))
R0_summ %>% ggplot(aes(t_start, R0)) + geom_line()
inf1 %>% group_by(i, HH) %>% arrange(HH, time) %>% select()

R0 <- NULL
for (j in unique(t$i)) {
  data <- t %>% filter(i==j) 
  
  full_data <- data %>%
    full_join(expand.grid(i=j, time=1:364)) %>%
    arrange(time) %>%
    replace_na(list(incid=0))
    
  if(nrow(data)==1) {
    R0 <- R0 %>% rbind(c(data$time+7, data$time+45, rep(0, ncol(R0)-3), j))
    next
  }
    
  t_range <- min(data$time):max(data$time)
  t_start <- t_range + 7
  t_start <- t_start[t_start<(max(data$time)-7)]
  t_end <- t_range + 45
  t_end <- t_end[t_end<=max(data$time)]
    
  if (length(t_end) < length(t_start)) {
    t_end <- c(t_end, rep(max(data$time), length(t_start) - length(t_end)))
  } 
    
  if (length(t_start) == 0 & length(t_end) == 0) {
    R0 <- R0 %>% rbind(c(data$time+7, data$time+45, rep(0, ncol(R0)-3), j))
    next
  }
    
  res <- estimate_R(full_data$incid,
                    method="parametric_si", 
                    config = list(t_start=t_start,
                                  t_end=t_end,
                                  mean_si = 28,
                                  std_si= 2))
    
  R0 <- R0 %>% rbind(res$R %>% mutate(i=j))
}
R0_summ_overall <- R0 %>% group_by(i, t_start) %>% summarise(R0=mean(`Mean(R)`,na.rm=T)) %>% 
  group_by(t_start) %>% summarise(R0=mean(R0, na.rm=T))

R0_summ_hh <- read_csv("r0/r0_hh_75h25c.csv") %>% group_by(i, t_start) %>% summarise(R0=mean(`Mean(R)`,na.rm=T)) %>% 
  group_by(t_start) %>% summarise(R0=mean(R0, na.rm=T))

R0_summ <- R0 %>% group_by(i, t_start) %>% summarise(R0=mean(`Mean(R)`,na.rm=T)) %>% 
  group_by(t_start) %>% summarise(R0=mean(R0, na.rm=T))

R0_summ_hh %>% ggplot(aes(t_start, R0)) + geom_line(aes(color="75% household transmission")) + 
  geom_line(data = R0_summ, aes(t_start, R0, color="Environmental"))

