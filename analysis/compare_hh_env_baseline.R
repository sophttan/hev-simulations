rm(list=ls())
gc()
library(tidyverse)

hh2 <- read_csv("results/hh_risk2.csv")
hh4 <- read_csv("results/hh_risk4.csv")
hh10 <- read_csv("results/hh_risk10.csv")

env <- read_csv("results/environmental.csv")

prep_results <- function(data){
  data %>% group_by(i, month) %>% 
    summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
    group_by(month) %>% summarise_all(mean) %>% select(!i)
}

hh2_summ <- hh2 %>% prep_results()
hh4_summ <- hh4 %>% prep_results()
hh10_summ <- hh10 %>% prep_results()
env_summ <- env %>% prep_results()

inf_type <- env_summ %>% full_join(hh2_summ, by="month") %>% full_join(hh4_summ, by="month") %>% full_join(hh10_summ, by="month")

p <- inf_type %>% ggplot(aes(month)) + 
  geom_line(aes(y=has_hh.x, color="Environmental transmission only")) + 
  geom_line(aes(y=has_hh.y, color="Household relative risk = 2")) + 
  geom_line(aes(y=has_hh.x.x, color="Household relative risk = 4")) + 
  geom_line(aes(y=has_hh.y.y, color="Household relative risk = 10")) + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Proportion of new cases with any household\ninfection in the last 30 days") +
  scale_color_brewer(palette="Dark2", direction=-1, 
                     breaks=c("Environmental transmission only","Household relative risk = 2",
                              "Household relative risk = 4","Household relative risk = 10")) +
  theme(legend.title = element_blank())
p %>% ggsave(filename = "results/hh_environ_comparison/prevalence20_comparison.jpg")

