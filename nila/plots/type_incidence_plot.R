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
    group_by(i_percent, month, TYPE) %>% summarise_all(mean) %>% select(!i)
  prepped
}

prep_inc = prep_results(inc)

plot_results <- function(data){
  p = ggplot(data = data) + geom_line(aes(x = month, y = incidence, colour = TYPE))
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Incidence over time (months)") +
    ylab("Incidence") + xlab("Months") +
    facet_grid(cols = vars(i_percent), labeller = as_labeller(cum_labels)) +
    labs(linetype="cumulative incidence (%)", colour="type of infection") +
    scale_color_discrete(name = "Type of Infection", labels = c("initial", "two causes of infection", "environmental", "person-to-person")) +
    theme_minimal()
}

plot_results(prep_inc)
