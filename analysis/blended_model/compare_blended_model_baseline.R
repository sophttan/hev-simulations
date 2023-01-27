rm(list=ls())
gc()
library(tidyverse)

setwd(here::here("results/blended_model/cumulative_inc_30/"))
hh31 <- read_csv("simulated_data/75p25e_hhrisk5.csv")
hh11 <- read_csv("simulated_data/50p50e_hhrisk5.csv")
hh13 <- read_csv("simulated_data/25p75e_hhrisk5.csv")

env <- read_csv(here::here("results/separate_models/cumulative_inc_30/simulated_data/environmental.csv"))
hh <- read_csv(here::here("results/separate_models/cumulative_inc_30/simulated_data/hh_risk5.csv"))

# months <- rep(1:13, each=35)[1:364]
# days_months <- data.frame(day=1:364, month=months)
# hh1 %>% left_join(days_months, by=c("time"="day"))

prep_results <- function(data){
  data %>% group_by(i, month) %>% 
    summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
    group_by(month) %>% summarise_all(mean) %>% select(!i)
}

hh31_summ <- hh31 %>% prep_results()
hh11_summ <- hh11 %>% prep_results()
hh13_summ <- hh13 %>% prep_results()

env_summ <- env %>% prep_results()
hh_summ <- hh %>% prep_results()

inf_type <- hh31_summ %>% 
  full_join(hh11_summ, by="month") %>% 
  full_join(hh13_summ, by="month") 

p <- inf_type %>% ggplot(aes(month)) + 
  geom_line(aes(y=has_hh.x, color="75:25 ratio")) + 
  geom_line(aes(y=has_hh.y, color="50:50 ratio")) + 
  geom_line(aes(y=has_hh, color="25:75 ratio")) + 
  geom_line(data = env_summ, aes(x=month, y=has_hh, color="Environmental only")) + 
  geom_line(data = hh_summ, aes(x=month, y=has_hh, color="Person-person only")) + 
  scale_x_continuous("Time (months)", breaks=1:13, expand=c(0,0)) +
  scale_y_continuous("Proportion of new cases with any household\ninfection in the last 7-45 days",
                     expand=c(0.01, 0), limits = c(0,1)) +
  scale_color_brewer(palette="Dark2", direction=-1, 
                     name="Ratio of transmission\nfrom person-person contact\nand environmental source",
                     breaks=c("Person-person only", "75:25 ratio","50:50 ratio","25:75 ratio", "Environmental only")) +
  theme(panel.grid.minor = element_blank()) +
  labs(title="Comparison of models at 30% cumulative incidence and\nhousehold relative risk of 5")
p
p %>% ggsave(filename = "figures/comparison_hhrisk5.jpg")
