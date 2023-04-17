setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

source("load_hh.R")

inc = load_hh(c(
  "data_inc_5.csv",
  "data_inc_10.csv",
  "data_inc_30.csv"
))

prep_results <- function(data){
  prepped_wide = data %>% mutate(month = TIME %/% 31) %>%
    mutate(has_hh = TYPE == "H", has_env = TYPE == "C") %>%
    group_by(i_percent, i, month) %>% summarise(has_hh = mean(has_hh), has_env = mean(has_env)) %>%
    group_by(i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped_long <- gather(prepped_wide, type, proportion, has_hh:has_env, factor_key=TRUE)
}

prep_df = prep_results(inc)

plot_results <- function(data){
  p = ggplot(data, aes(x = month, y = proportion, colour = type)) +
    geom_line()
  p +
    ggtitle("proportion of infection type by cumulative incidence (%)") +
    facet_grid(cols = vars(i_percent)) +
    scale_color_discrete(name = "Type of transmission", label = c("person-to-person", "environmental")) +
    theme_minimal()
}

plot_results(prep_df)
