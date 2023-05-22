library(tidyverse)
library(here)

source(here("nila", "plots", "loading-files", "hh.R"))

prep_type_prop <- function(data){
  prepped_wide = data %>% mutate(month = TIME %/% 31) %>%
    mutate(has_hh = TYPE == "H", has_env = TYPE == "C") %>%
    group_by(i_percent, i, month) %>% summarise(has_hh = mean(has_hh), has_env = mean(has_env)) %>%
    group_by(i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped_long <- gather(prepped_wide, type, proportion, has_hh:has_env, factor_key=TRUE)
}

plot_type_prop <- function(data){
  p = ggplot(data, aes(x = month, y = proportion, colour = type)) +
    geom_line()
  p +
    ggtitle("proportion of infection type by cumulative incidence (%)") +
    facet_grid(cols = vars(i_percent)) +
    scale_color_discrete(name = "Type of transmission", label = c("person-to-person", "environmental")) +
    theme_minimal()
  
  ggsave("type_proportion_plot.png",
         path = here("nila", "plots", "images"),
         width = 2000, height = 1500, units = "px")
}

type_proportion <- function(data){
  prepped = prep_type_prop(data)
  plot_type_prop(prepped)
}

type_proportion(inc)
