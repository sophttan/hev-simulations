setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

source("load_hh.R")

inc = load_hh(c(
  "data_inc_5.csv",
  "data_inc_10.csv",
  "data_inc_30.csv"
))

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
  gg_plot = p + 
    ggtitle("cumulative incidence over time by final cumulative incidence (%)") +
    facet_grid(cols = vars(i_percent)) +
    ylab("cumulative incidence") +
    labs(linetype="cumulative incidence (%)", colour="type of infection") +
    scale_color_discrete(name = "type of infection", labels = c("initial", "two causes of infection", "environmental", "person-to-person")) +
    theme_minimal()
  gg_plot
}

plot_results(prep_df)
