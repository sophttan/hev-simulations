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
