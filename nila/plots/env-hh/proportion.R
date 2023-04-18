library(tidyverse)
library(here)

source(here("nila", "plots", "loading-files", "env_hh.R"))

# getting proportion of person-to-person transmission for person-to-person
prep_inc_prop <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>%
    mutate(has_hh = TYPE == "H") %>%
    group_by(enviro, i_percent, i, month) %>% summarise(has_hh = mean(has_hh)) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# getting proportion of person-to-person transmission for evironmental transmission 
prep_env_prop <- function(data){
  prepped= data %>% mutate(month = time %/% 31) %>%
    group_by(enviro, i_percent, i, month) %>% summarise(has_hh = mean(has_hh)) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# combining together functions to get one single dataset
prepping_prop <- function(hh, env){
  prep_inc = prep_inc_prop(hh)
  prep_env = prep_env_prop(env)
  prepped <- rbind(prep_inc, prep_env)
  prepped
}

# plotting the final dataset
env_hh_prop_gg <- function(data){
  p = ggplot(data, aes(x = month, y = has_hh, colour = enviro)) +
    geom_line()
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Proportion of household cases over time (months)") +
    ylab("New household cases in last 7-45 days (%)") +
    facet_grid(cols = vars(i_percent),  labeller = as_labeller(cum_labels)) +
    scale_color_discrete(name = "Type of transmission") +
    theme_minimal()
  
  ggsave("env_hh_prop_plot.png",
         path = here("nila", "plots", "images"),
         width = 2000, height = 1500, units = "px")
}

# combined function to do this all together
env_hh_prop_plot <- function(hh, env){
  prep = prepping_prop(hh, env)
  env_hh_prop_gg(prep)
}

env_hh_prop_plot(hh, env)