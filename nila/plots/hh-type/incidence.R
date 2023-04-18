library(tidyverse)
library(here)

source(here("nila", "plots", "loading-files", "hh.R"))

prep_type_inc <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(i_percent, i, month, TYPE) %>%  
    summarise(incidence = n()/1000) %>%
    group_by(i_percent, month, TYPE) %>% summarise_all(mean) %>% select(!i)
  prepped
}

plot_type_inc <- function(data){
  p = ggplot(data = data) + geom_line(aes(x = month, y = incidence, colour = TYPE))
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  gg_plot = p +
    ggtitle("Incidence over time (months)") +
    ylab("Incidence") + xlab("Months") +
    facet_grid(cols = vars(i_percent), labeller = as_labeller(cum_labels)) +
    labs(linetype="cumulative incidence (%)", colour="type of infection") +
    scale_color_discrete(name = "Type of Infection", labels = c("initial", "two causes of infection", "environmental", "person-to-person")) +
    theme_minimal()
  
  ggsave("type_incidence_plot.png",
         path = here("nila", "plots", "images"),
         width = 2000, height = 1500, units = "px")
}

type_incidence <- function(data){
  prepped = prep_type_inc(data)
  plot_type_inc(prepped)
}

type_incidence(inc)
