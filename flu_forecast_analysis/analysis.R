source("./simulation_functions.R")
library(dplyr)
library(ggplot2)
library(cmdstanr)

cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/normal_quantiles.stan')

forecast <- 
  read.csv("../../forecast-hub/FluSight-forecast-hub/model-output/UM-DeepOutbreak/2023-11-18-UM-DeepOutbreak.csv")

forc <- forecast %>% 
  filter(output_type == "quantile", horizon == 0, location == "US") %>%
  mutate(prob = as.numeric(output_type_id), quantile = log(value + 1))
  

forecast %>% 
  filter(output_type == "quantile", horizon == 0, location == "US") %>%
  mutate(prob = as.numeric(output_type_id), quantile = log(value + 1)) %>% 
  ggplot() +
  geom_point(aes(x = prob, y = quantile))
  
stan_data <- forecast %>% 
  filter(output_type == "quantile", horizon == 0, location == "US") %>%
  mutate(prob = as.numeric(output_type_id), quantile = log(value + 1)) %>% 
  select(quantile, prob) %>% 
  make_stan_data(size = 40, comps = 3)



fit <- stan_fit_draws(cltmod, stan_data,
                          sampler = "variational",
                          elbo = 300, grad = 10,
                          out_s = 2000)



fqs <- ecdf(fit$draws)(stan_data$Q)








