rm(list=ls())
gc()
library(tidyverse)

setwd(here::here("results/separate_models/cumulative_inc_5/"))

env <- read_csv("simulated_data/environmental.csv")
env2 <- read_csv("simulated_data/environmental_fixedhhrisk.csv")

months <- rep(1:13, each=35)[1:364]
days_months <- data.frame(day=1:364, month=months)

prep_results <- function(data){
  data %>% group_by(i, month) %>% 
    summarise(has_hh = mean(has_hh,na.rm=T)) %>% 
    group_by(month) %>% summarise_all(mean) %>% select(!i)
}

env_summ <- env %>% prep_results()
env2_summ <- env2 %>% prep_results()

inf_type <- env_summ %>% 
  full_join(env2_summ, by="month")

p <- inf_type %>% ggplot(aes(month)) + 
  geom_line(aes(y=has_hh.x, color="Environmental transmission\nto individuals\n")) + 
  geom_line(aes(y=has_hh.y, color="Environmental transmission\nto households")) + 
  scale_x_continuous("Time (months)") +
  scale_y_continuous("Proportion of new cases with any household\ninfection in the last 7-45 days", limits=c(0, 0.15)) +
  scale_color_brewer(palette="Dark2", direction=-1, 
                     breaks=c("Environmental transmission\nto individuals\n",
                              "Environmental transmission\nto households")) +
  theme(legend.title = element_blank()) +
  labs(title="Comparison of models at 5% cumulative incidence")
p %>% ggsave(filename = "figures/comparison_env.jpg")


