setwd("C:/Users/water/OneDrive/Documents/ucsf/hev-simulations/nila")

library(tidyverse)

source("load_env_hh.R")

hh_relpath = c("data_inc_5.csv", "data_inc_10.csv", "data_inc_30.csv")
env_relpath =
  c("../results/separate_models/cumulative_inc_5/simulated_data/environmental.csv",
    "../results/separate_models/cumulative_inc_10/simulated_data/environmental.csv",
    "../results/separate_models/cumulative_inc_30/simulated_data/environmental.csv")

hh_env = load_env_hh(hh_relpath, env_relpath)
hh = hh_env[[1]]
env = hh_env[[2]]

# calculating incidence for person-to-person data
prep_hh_results <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(enviro, i_percent, i, month) %>%  
    summarise(incidence = n()) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# calculating incidence for environmental data
prep_env_results <- function(data){
  prepped = data %>% mutate(month = time %/% 31) %>% group_by(enviro, i_percent, i, month) %>%  
    summarise(incidence = n()) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# combining together functions to get one single dataset
prepping <- function(hh, env){
  prep_hh = prep_hh_results(inc)
  prep_env = prep_env_results(env)
  prepped <- rbind(prep_hh, prep_env)
  prepped
}

# # cumulative incidence for environmental data
# prep_env_cum <- function(data){
#   prepped = data %>% mutate(month = time %/% 31) %>% group_by(enviro, i_percent, i, month) %>%  
#     summarise(incidence = n()/1000) %>%
#     group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i) %>%
#     group_by(enviro, i_percent) %>% mutate(cum_sum = cumsum(incidence))
#   prepped
# }

plot_results <- function(data){
  p = ggplot(data = data) + geom_line(aes(x = month, y = incidence, colour = enviro))
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Incidence over time (months)") +
    ylab("Incidence") + xlab("Months") +
    facet_grid(cols = vars(i_percent), labeller = as_labeller(cum_labels)) +
    scale_color_discrete(name = "Type of transmission", labels = c("Environmental only", "Person-to-person only")) +
    theme_minimal()
}

# combined function to do this all together
env_hh_inc_plot <- function(hh, env){
  prep = prepping(inc, env)
  plotly::ggplotly(plot_results(prep))
}

env_hh_inc_plot(hh, env)
