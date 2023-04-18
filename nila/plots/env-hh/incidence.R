library(tidyverse)
library(here)

source(here("nila", "plots", "loading-files", "env_hh.R"))

# calculating incidence for person-to-person data
prep_hh_inc <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(enviro, i_percent, i, month) %>%  
    summarise(incidence = n()) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# calculating incidence for environmental data
prep_env_inc <- function(data){
  prepped = data %>% mutate(month = time %/% 31) %>% group_by(enviro, i_percent, i, month) %>%  
    summarise(incidence = n()) %>%
    group_by(enviro, i_percent, month) %>% summarise_all(mean) %>% select(!i)
  prepped
}

# combining together functions to get one single dataset
inc_prepping <- function(hh, env){
  prep_hh = prep_hh_inc(hh)
  prep_env = prep_env_inc(env)
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

env_hh_inc_results <- function(data){
  p = ggplot(data = data) + geom_line(aes(x = month, y = incidence, colour = enviro))
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Incidence over time (months)") +
    ylab("Incidence") + xlab("Months") +
    facet_grid(cols = vars(i_percent), labeller = as_labeller(cum_labels)) +
    scale_color_discrete(name = "Type of transmission", labels = c("Environmental only", "Person-to-person only")) +
    theme_minimal()
  
  ggsave("env_hh_inc_plot.png",
         path = here("nila", "plots", "images"),
         width = 2000, height = 1500, units = "px")
}

# combined function to do this all together
env_hh_inc_plot <- function(hh, env){
  prep = inc_prepping(hh, env)
  env_hh_inc_results(prep)
}

env_hh_inc_plot(hh, env)
