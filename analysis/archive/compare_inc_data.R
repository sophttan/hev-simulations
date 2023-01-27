hh <- read_csv("/Users/sophiatan/Dropbox/ID Modeling/Transmission Methods Modeling/inc_data_inc30_rr2_hh.csv")
env <- read_csv("/Users/sophiatan/Dropbox/ID Modeling/Transmission Methods Modeling/inc_data_inc30_env.csv")
both <- read_csv("/Users/sophiatan/Dropbox/ID Modeling/Transmission Methods Modeling/inc_data_inc30_rr2_blend.csv")

summ_hh <- hh %>% mutate(week=floor(time/7)) %>% 
  group_by(week) %>% summarise(count=n()) %>% full_join(expand.grid(week=0:12)) %>% replace_na(list(count=0))
summ_env <- env %>% mutate(week=floor(time/7)) %>% 
  group_by(week) %>% summarise(count=n()) %>% full_join(expand.grid(week=0:12)) %>% replace_na(list(count=0))
summ_blend <- both %>% mutate(week=floor(time/7)) %>% 
  group_by(week) %>% summarise(count=n()) %>% full_join(expand.grid(week=0:12)) %>% replace_na(list(count=0))

summ_hh %>% ggplot(aes(week, count)) + geom_line(aes(color="Person-person")) + 
  geom_line(data=summ_env, aes(week, count, color="Environment")) + 
  geom_line(data=summ_blend, aes(week, count, color="Both")) + 
  scale_x_continuous("Time (weeks)", breaks=seq(0,13,1), expand=c(0.01,0)) +
  scale_y_continuous("Number of infections", limits=c(0,30), expand=c(0.01,0)) + 
  scale_color_brewer(palette="Dark2", direction=-1, breaks=c("Person-person", "Environment", "Both")) +
  labs(title="Weekly infections over 90 days (5-6% cumulative incidence)") + 
  theme(legend.title = element_blank())

hh <- hh %>% group_by(HH) %>% mutate(day_limits = list(time)) %>% 
  ungroup() %>% rowwise() %>% 
  mutate(has_hh = any((time - unlist(day_limits)) < 45 & (time - unlist(day_limits)) > 7)) %>% 
  select(!day_limits)
env <- env %>% group_by(HH) %>% mutate(day_limits = list(time)) %>% 
  ungroup() %>% rowwise() %>% 
  mutate(has_hh = any((time - unlist(day_limits)) < 45 & (time - unlist(day_limits)) > 7)) %>% 
  select(!day_limits)
both <- both %>% group_by(HH) %>% mutate(day_limits = list(time)) %>% 
  ungroup() %>% rowwise() %>% 
  mutate(has_hh = any((time - unlist(day_limits)) < 45 & (time - unlist(day_limits)) > 7)) %>% 
  select(!day_limits)

mean(hh$has_hh)
mean(env$has_hh)
mean(both$has_hh)
