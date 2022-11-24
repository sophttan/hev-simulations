rm(list=ls())
gc()
library(tidyverse)

setwd("results/cumulative_inc_30/")
hh1 <- read_csv("simulated_data/hh_risk1.csv")
hh2 <- read_csv("simulated_data/hh_risk2.csv")
hh5 <- read_csv("simulated_data/hh_risk5.csv")

env <- read_csv("simulated_data/environmental.csv")

prep_results <- function(data){
  data %>% group_by(i, month) %>% 
    summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
    group_by(month) %>% summarise_all(mean) %>% select(!i)
}

hh1_summ <- hh1 %>% prep_results()
hh2_summ <- hh2 %>% prep_results()
hh5_summ <- hh5 %>% prep_results()
env_summ <- env %>% prep_results()

inf_type <- env_summ %>% 
  full_join(hh1_summ, by="month") %>% 
  full_join(hh2_summ, by="month") %>% 
  full_join(hh5_summ, by="month")

p <- inf_type %>% ggplot(aes(month)) + 
  geom_line(aes(y=has_hh.x, color="Environmental transmission only")) + 
  geom_line(aes(y=has_hh.y, color="Household relative risk = 1")) + 
  geom_line(aes(y=has_hh.x.x, color="Household relative risk = 2")) + 
  geom_line(aes(y=has_hh.y.y, color="Household relative risk = 5")) + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Proportion of new cases with any household\ninfection in the last 7-45 days") +
  scale_color_brewer(palette="Dark2", direction=-1, 
                     breaks=c("Environmental transmission only","Household relative risk = 1",
                              "Household relative risk = 2","Household relative risk = 5")) +
  theme(legend.title = element_blank()) +
  labs(title="Comparison of models at 30% cumulative incidence")
p %>% ggsave(filename = "figures/comparison.jpg")

