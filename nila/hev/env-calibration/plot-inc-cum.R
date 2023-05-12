library(tidyverse)
library(here)

env_relpath = here("nila", "hev", "env-calibration", "data", c("env-results-5.csv", "env-results-10.csv", "env-results-30.csv"))

# loading environmental datasets
env_5 = read.csv(env_relpath[1]) %>% mutate(i_percent = 5)
env_10 = read.csv(env_relpath[2]) %>% mutate(i_percent = 10)
env_30 = read.csv(env_relpath[3]) %>% mutate(i_percent = 30)

# binding together all three environment files
env = rbind(env_5, env_10, env_30)
env$i_percent <- as.factor(env$i_percent)

# calculating incidence + cumulative incidence for environmental data
prep_env <- function(data){
  prepped = data %>% mutate(month = TIME %/% 31) %>% group_by(i_percent, i, month) %>%  
    summarise(incidence = n()) %>% mutate(prop_inc = incidence/1000) %>%
    group_by(i_percent, month) %>% summarise_all(mean) %>% select(!i) %>%
    group_by(i_percent) %>% mutate(cum_sum = cumsum(incidence))
  prepped
}

env_inc_results <- function(data){
  p = ggplot(data = data) + geom_line(aes(x = month, y = incidence))
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Incidence over time (months)") +
    ylab("Incidence") + xlab("Months") +
    facet_grid(cols = vars(i_percent), labeller = as_labeller(cum_labels)) +
    theme_minimal()
  
  ggsave("inc-plot.png",
         path = here("nila", "hev", "env-calibration", "plots"),
         width = 2000, height = 1500, units = "px")
}

env_cum_results <- function(data){
  p = ggplot(data = data) + geom_line(aes(x = month, y = cum_sum))
  
  cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")
  
  p +
    ggtitle("Cumulative incidence over time (months)") +
    ylab("Cumulative incidence") + xlab("Months") +
    facet_grid(cols = vars(i_percent), labeller = as_labeller(cum_labels)) +
    theme_minimal()
  
  ggsave("cum-plot.png",
         path = here("nila", "hev", "env-calibration", "plots"),
         width = 2000, height = 1500, units = "px")
}

# combined function to do this all together
env_inc_plot <- function(env){
  prep = prep_env(env)
  env_inc_results(prep)
}

env_cum_plot <- function(env){
  prep = prep_env(env)
  env_cum_results(prep)
}

env_inc_plot(env)
env_cum_plot(env)