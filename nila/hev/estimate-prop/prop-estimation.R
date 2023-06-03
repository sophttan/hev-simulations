set.seed(1234)

gc()
library(foreach)
library(doParallel)
library(tidyverse)
library(here)
library(latex2exp)

# choosing the number of cores to do parallelisation on
numCores = detectCores()-2
registerDoParallel(numCores)

# getting the big dataset of blended simulations
blended = read.csv(here("nila/hev/estimate-prop/blend-sims/blended.csv"))
source(here::here("nila/hev/estimate-prop/r_hh-r-functions.R"))

#################################
# proportion of household cases #
#################################

# R_hh/R #
##########

source(here::here("nila/hev/estimate-prop/true-prop.R"))
prop_true = true_prop(blended)$avg_prop

r = estimate_prop(blended)
r_hh = r$R_hh / r$R

parameters = c(1:3,1:3,1:3)*0.25

# 3 tables with 25%, 50%, and 70%
incidence = c(5,5,5,10,10,10,30,30,30)
prop = data.frame(incidence, parameters, prop_true, r_hh)

# Interval estimate #
#####################

source(here::here("nila/hev/estimate-prop/interval-functions.R"))

interval_prop = interval_blend(blended, sims)
prop$interval = interval_prop$avg_prop

# adding env + ptp data estimates #
###################################

# ptp = person-to-person
# env = environment

prop_5 = prop %>% slice(1:3)
prop_10 = prop %>% slice(4:6)
prop_30 = prop %>% slice(7:9)
prop_list = list(prop_5, prop_10, prop_30)

# we need to get the environmental and only person-to-person models as well
# note that the true value for environmental and person-to-person are 0 and 1 respectively

# ~environmental~

# loading datasets
env_inc_5 <- read.csv(here("nila/hev/env-calibration/data/env-results-5.csv")) %>% mutate(inc = 5)
env_inc_10 <- read.csv(here("nila/hev/env-calibration/data/env-results-10.csv")) %>% mutate(inc = 10)
env_inc_30 <- read.csv(here("nila/hev/env-calibration/data/env-results-30.csv")) %>% mutate(inc = 30)
env = rbind(env_inc_5, env_inc_10, env_inc_30) %>% rename(SIZE = HHsize, TYPE = Type, ID = No)

est_env = estimate_int_rhh(env, sims, 0, 0)
prop = rbind(prop, est_env)

# ~person-to-person~

hh_relpath = here("nila", "plots", "data", c("inc_5.csv", "inc_10.csv", "inc_30.csv"))
ptp_5 = read.csv(hh_relpath[1]) %>% mutate(inc = 5) %>% filter(i <= 500)
ptp_10 = read.csv(hh_relpath[2]) %>% mutate(inc = 10) %>% filter(i <= 500)
ptp_30 = read.csv(hh_relpath[3]) %>% mutate(inc = 30) %>% filter(i <= 500)
ptp = rbind(ptp_5, ptp_10, ptp_30)
ptp$inc <- as.factor(ptp$inc)

est_ptp = estimate_int_rhh(ptp, sims, 1, true_hh(ptp))
prop = rbind(prop, est_ptp)

############
# Plotting #
############

cum_labels = c(`5` = "cumulative incidence = 5%", `10` = "10%", `30` = "30%")

prop_pretty = prop %>% pivot_longer(prop, cols = 3:5, names_to = "type")
prop_pretty$type = str_replace(prop_pretty$type, "prop_true", "true")

prop_pretty %>%
  ggplot(aes(x = parameters, y = value, colour = type)) +
  geom_line() +
  facet_grid(cols = vars(incidence), labeller = as_labeller(cum_labels)) +
  ylab("proportion of household cases") +
  xlab("proportion of person-to-person") +
  scale_color_discrete(name = "type of estimate", labels = c("interval", TeX("$R_{hh}/R$"), "true")) +
  ggtitle(TeX("estimating via $R_{hh}/R$ with a log normal infection interval")) +
  theme_minimal()

ggsave("estimate-line-plot.png",
       path = here("nila", "hev", "estimate-prop"),
       width = 2500, height = 1500,
       units = "px")