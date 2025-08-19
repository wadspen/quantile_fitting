setwd("~/quantile_fitting/")
source("./simulation/simulation_functions.R")
library(dplyr)
library(ggplot2)
library(stringr)
library(distr)
library(tidyverse)
library(cmdstanr)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)

qgp_stan <- cmdstan_model(stan_file = 
                            './stan_models/normal_t_mix4_quantiles.stan')

burn <- 1000
sample <- 1000

mod_loc <- "~/forecast-hub/FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

mod <- models[30]
date <- sub_dates[14]
samp_mods <- sample(length(models), 16)
# 37 28 17  4 31  3 12 13 39  7 34 11
h <- 0
all_forecasts <- data.frame()
for (i in samp_mods) {
  forc_file <- list.files(paste0(mod_loc, models[i]), pattern = date)
  forecasts <- read.csv(paste0(mod_loc, models[i], "/", forc_file)) %>% 
    filter(horizon == h)
  
  forecasts$model <- models[i]
  all_forecasts <- rbind(all_forecasts, forecasts)
}


all_forecasts %>% 
  filter(location == "US") %>% 
  mutate(output_type_id = as.numeric(output_type_id)) %>% 
  filter(output_type_id > 0, 
         !(model %in% c("ISU_NiemiLab-NLH", "ISU_NiemiLab-SIR"))) %>% 
  ggplot() +
  geom_point(aes(x = output_type_id, y = log(value + 1))) + 
  scale_x_continuous(breaks=seq(0, 1, .5)) + 
  xlab("Probability") +
  ylab("log(hospitalizations + 1)") +
  facet_wrap(~model, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=18),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none")

all_forecasts %>% 
  filter(location == "12", model == "UGA_flucast-Copycat") %>% 
  mutate(output_type_id = as.numeric(output_type_id)) %>% 
  ggplot() +
  geom_point(aes(x = output_type_id, y = log(value + 1)))

forecast %>%
  ggplot() +
  geom_point(aes(x = output_type_id, y = log(value + 1)))


quantiles <- log(as.numeric(forecast$value) + 1)
probs <- as.numeric(forecast$output_type_id)
dat <- data.frame(quantile = quantiles,
                  prob = probs)

dat <- dat %>% 
  filter(quantile != 0)
# dat <- dat[-c(21:23),]

stan_data <- make_stan_data(dat, size = 1, comps = 5)
stan_samps <- qgp_stan$sample(data = stan_data,
                              iter_warmup = burn,
                              iter_sampling = sample,
                              chains = 1,
                              adapt_delta = .9999,
                              refresh = 100)

draws <- stan_samps$draws(format = "df")


preds <- draws[, str_detect(colnames(draws), "pred_q")]

stuff <- apply(samp_quantiles, MARGIN = 2, FUN = median, na.rm = TRUE)
# stuff <- quantile(draws$dist_samps, probs = probs)
plot(stuff[-1] ~ stan_data$Q[-1]); abline(a = 0 , b = 1)
apply(samp_quantiles, MARGIN = 2, FUN = quantile, probs = c(.025, .975), na.rm = TRUE)

all_pis <- draws[, str_detect(colnames(draws), "pi")]
all_mus <- draws[, str_detect(colnames(draws), "mus")]
all_sigmas <- draws[, str_detect(colnames(draws), "sigmas")]
all_ns <- draws$n

num_samps <- nrow(all_pis)
m <- 1
samp_quantiles <- matrix(NA, nrow = num_samps, ncol = length(probs))
repeat{
  num <- sample(num_samps, 1)
  mus <- unlist(all_mus[num,])
  sigmas <- unlist(all_sigmas[num,])
  pi <- unlist(all_pis[num,])
  pi[which(pi < 0)] <- 0
  if (sum(pi) < 1) {
    pi[which.max(pi)] <- pi[which.max(pi)] + (1 - sum(pi))
  } else if (sum(pi) > 1) {
    pi[which.min(pi)] <- pi[which.min(pi)] + (1 - sum(pi))
  }
  n <- unlist(all_ns[num])
  
  normmix <- UnivarMixingDistribution(Norm(mus[1], sigmas[1]),
                                      Norm(mus[2], sigmas[2]), 
                                      Norm(mus[3], sigmas[3]), 
                                      Norm(mus[4], sigmas[4]),
                                      Norm(mus[5], sigmas[5]), 
                                      mixCoeff = pi)
  samp <- r(normmix)(n)
  samp_quantiles[m,] <- quantile(samp, probs = probs)
  
  m <- m + 1
  if (m > num_samps) {break}
  print(m)
}


quant_bounds <- apply(samp_quantiles, MARGIN = 2, 
                      FUN = quantile, 
                      probs = probs, na.rm = TRUE)

quant_bounds <- data.frame(t(quant_bounds))
colnames(quant_bounds) <- as.character(probs)
quant_bounds$prob <- probs
quant_bounds$quantile <- quantiles

quant_bounds <- quant_bounds %>% 
  mutate(cover98 = between(quantile, `0.01`, `0.99`),
         cover95 = between(quantile, `0.025`, `0.975`), 
         cover90 = between(quantile, `0.05`, `0.95`),
         cover80 = between(quantile, `0.1`, `0.9`),
         cover70 = between(quantile, `0.15`, `0.85`),
         cover60 = between(quantile, `0.2`, `0.8`),
         cover50 = between(quantile, `0.25`, `0.75`),
         cover40 = between(quantile, `0.3`, `0.7`),
         cover30 = between(quantile, `0.35`, `0.65`),
         cover20 = between(quantile, `0.4`, `0.6`),
         cover10 = between(quantile, `0.45`, `0.55`))


quant_bounds %>% 
  filter(between(prob, .05, .95)) %>%
  ggplot() +
  geom_segment(aes(x = prob, y = `0.025`, yend = `0.975`),
               size = 1.5) +
  geom_segment(aes(x = prob, y = `0.25`, yend = `0.75`), 
               colour = "red",
               size = 10) +
  geom_point(aes(x = prob, y = `0.5`), colour = "pink",
             fill = "pink", shape = 24, size = 2) +
  geom_point(aes(x = prob, y = quantile), size = 2,
             colour = "violet") +
  theme_bw()
  

apply(quant_bounds, MARGIN = 2, FUN = mean)







