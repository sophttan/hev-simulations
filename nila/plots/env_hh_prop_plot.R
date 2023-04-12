setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

source("load_env_hh.R")

hh_relpath = c("data_inc_5.csv", "data_inc_10.csv", "data_inc_30.csv")
env_relpath =
  c("../results/separate_models/cumulative_inc_5/simulated_data/environmental.csv",
    "../results/separate_models/cumulative_inc_10/simulated_data/environmental.csv",
    "../results/separate_models/cumulative_inc_30/simulated_data/environmental.csv")

hh_env = load_env_inc(hh_relpath, env_relpath)
hh = hh_env[[1]]
env = hh_env[[2]]

# getting proportion of person-to-person transmission for person-to-person
prep_inc_results <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>%
    mutate(has_hh = TYPE == "H") %>%
    group_by(enviro, i_percent, i, month) %>% summarise(has_hh = mean(has_hh)) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# getting proportion of person-to-person transmission for evironmental transmission 
prep_env_results <- function(data){
  prepped= data %>% mutate(month = time %/% 31) %>%
    group_by(enviro, i_percent, i, month) %>% summarise(has_hh = mean(has_hh)) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# combining together functions to get one single dataset
prepping <- function(hh, env){
  prep_inc = prep_inc_results(inc)
  prep_env = prep_env_results(env)
  prepped <- rbind(prep_inc, prep_env)
  prepped
}

# plotting the final dataset
plot_results <- function(data){
  p = ggplot(data, aes(x = month, y = has_hh, colour = enviro)) +
    geom_line()
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Proportion of household cases over time (months)") +
    ylab("New household cases in last 7-45 days (%)") +
    facet_grid(cols = vars(i_percent),  labeller = as_labeller(cum_labels)) +
    scale_color_discrete(name = "Type of transmission") +
    theme_minimal()
}

# combined function to do this all together
env_hh_prop_plot <- function(hh, env){
  prep = prepping(inc, env)
  plotly::ggplotly(plot_results(prep))
}

env_hh_prop_plot(hh, env)
