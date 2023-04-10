library(tidyverse)

df = read.csv("data_inc_10.csv")
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
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(i, month, TYPE) %>%  
    summarise(incidence = n()/1000) %>%
    group_by(month, TYPE) %>% summarise_all(mean) %>% select(!i) %>%
    group_by(TYPE) %>% mutate(cum_sum = cumsum(incidence))
  prepped
}

prep_df = prep_results(df)

plot_results <- function(data){
  p = ggplot(data, aes(x = month, y = cum_sum, colour = TYPE)) +
    geom_line()
  p + 
    ylab("cumulative incidence") +
    scale_color_discrete(name = "type of infection", labels = c("initial", "two causes of infection", "environmental", "person-to-person"))
}

plot_results(prep_df)
