library(tidyverse)
library(here)

source(here("nila", "plots", "loading-files", "hh.R"))

prep_type_cum <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(i_percent, i, month, TYPE) %>%  
    summarise(incidence = n()/1000) %>%
    group_by(i_percent, month, TYPE) %>% summarise_all(mean) %>% select(!i) %>%
    group_by(i_percent, TYPE) %>% mutate(cum_sum = cumsum(incidence))
  prepped
}

plot_type_cum <- function(data){
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
  
  ggsave("type_cum_inc_plot.png",
         path = here("nila", "plots", "images"),
         width = 2000, height = 1500, units = "px")
}

type_cum_inc <- function(data){
  prepped = prep_type_cum(data)
  plot_type_cum(prepped)
}

type_cum_inc(inc)
