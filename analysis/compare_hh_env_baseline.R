rm(list=ls())
gc()
library(tidyverse)

setwd(here::here("results/separate_models/cumulative_inc_5/"))
#hh1 <- read_csv("simulated_data/hh_risk1.csv")
hh2 <- read_csv("simulated_data/hh_risk2.csv")
hh5 <- read_csv("simulated_data/hh_risk5.csv")

env <- read_csv("simulated_data/environmental.csv")

months <- rep(1:13, each=35)[1:364]
days_months <- data.frame(day=1:364, month=months)

prep_results <- function(data){
  data %>% group_by(i, month) %>% 
    summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
    group_by(month) %>% summarise_all(mean) %>% select(!i)
}

#hh1_summ <- hh1 %>% prep_results()
hh2_summ <- hh2 %>% prep_results()
hh5_summ <- hh5 %>% prep_results()
env_summ <- env %>% prep_results()

inf_type <- env_summ %>% 
#  full_join(hh1_summ, by="month") %>% 
  full_join(hh2_summ, by="month") %>% 
  full_join(hh5_summ, by="month")

p <- inf_type %>% ggplot(aes(month)) + 
  geom_line(aes(y=has_hh.x, color="Environmental transmission only")) + 
#  geom_line(aes(y=has_hh.y, color="Household relative risk = 1")) + 
  geom_line(aes(y=has_hh.y, color="Household relative risk = 2")) + 
  geom_line(aes(y=has_hh, color="Household relative risk = 5")) + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Proportion of new cases with any household\ninfection in the last 7-45 days") +
  scale_color_brewer(palette="Dark2", direction=-1, 
                     breaks=c("Environmental transmission only",
                              "Household relative risk = 2","Household relative risk = 5")) +
  theme(legend.title = element_blank()) +
  labs(title="Comparison of models at 5% cumulative incidence")
p %>% ggsave(filename = "figures/comparison.jpg")


mean((hh1 %>% group_by(i) %>% summarise(HH=length(unique(HH))))$HH)/200
mean((hh2 %>% group_by(i) %>% summarise(HH=length(unique(HH))))$HH)/200
mean((hh5 %>% group_by(i) %>% summarise(HH=length(unique(HH))))$HH)/200
mean((env %>% group_by(i) %>% summarise(HH=length(unique(HH))))$HH)/200

hh1_test <- hh1 %>% select(!month) %>% left_join(days_months, by=c("time"="day"))
hh1_test %>% filter(!Type %>% is.na()) %>% prep_results() %>% ggplot(aes(month, has_hh)) + geom_line() + 
  scale_x_continuous("Time (months - 35 days each)") + 
  scale_y_continuous("Proportion of new cases with any household\ninfection in the last 7-45 days")
hh1_test %>% group_by(i, month) %>% 
  summarise(cases=n(), has_hh = mean(has_hh,na.rm=T)) %>% 
  group_by(month) %>% summarise(cases=sum(cases)/100, has_hh=mean(has_hh)) 
