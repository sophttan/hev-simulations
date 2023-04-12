setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

inc_5 = read.csv("data_inc_5.csv") %>% mutate(i_percent = 5)
inc_10 = read.csv("data_inc_10.csv") %>% mutate(i_percent = 10)
inc_30 = read.csv("data_inc_30.csv") %>% mutate(i_percent = 30)
inc = rbind(inc_5, inc_10, inc_30)
inc$i_percent <- as.factor(inc$i_percent)

# Create frame for storing results.
# ID: ID of individual.
# SIZE: size of individual's household.
# HH: ID of individual's household.
# TYPE: the kind of infection: household (H), community (C), or both (B).
# TIME: when the individual became infectious.
# S_num: number of susceptible people in individual's household when their 
#        infectious period begins.
# I_num: number of people in household that this individual infected over 
#        their infectious period.
# i: ID of the simulation

prep_results <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(i_percent, i, month, TYPE) %>%  
    summarise(incidence = n()/1000) %>%
    group_by(i_percent, month, TYPE) %>% summarise_all(mean) %>% select(!i) %>%
    group_by(i_percent, TYPE) %>% mutate(cum_sum = cumsum(incidence))
  prepped
}

prep_df = prep_results(inc)

plot_results <- function(data){
  p = ggplot(data, aes(x = month, y = cum_sum, colour = TYPE)) +
    geom_line()
  p + 
    ggtitle("cumulative incidence over time by final cumulative incidence (%)") +
    facet_grid(cols = vars(i_percent)) +
    ylab("cumulative incidence") +
    labs(linetype="cumulative incidence (%)", colour="type of infection") +
    scale_color_discrete(name = "type of infection", labels = c("initial", "two causes of infection", "environmental", "person-to-person")) +
    theme_minimal()
}

plot_results(prep_df)
