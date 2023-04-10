library(tidyverse)

df = read.csv("data_inc_10.csv")

prep_results <- function(data){
  prepped_wide = data %>% mutate(month = TIME %/% 31) %>%
    mutate(has_hh = TYPE == "H", has_env = TYPE == "C") %>%
    group_by(i, month) %>% summarise(has_hh = mean(has_hh), has_env = mean(has_env)) %>%
    group_by(month) %>% summarise_all(mean) %>% select(!i)
  prepped_long <- gather(prepped_wide, type, proportion, has_hh:has_env, factor_key=TRUE)
}

prep_df = prep_results(df)

plot_results <- function(data){
  p = ggplot(data, aes(x = month, y = proportion, colour = type)) +
    geom_line()
  p + scale_color_discrete(name = "type of infection", label = c("person-to-person", "environmental"))
}

plot_results(prep_df)
